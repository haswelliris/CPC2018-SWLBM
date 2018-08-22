#include "../header/LbmSlave.h"
#include "../header/Config.h"

#include "slave.h"
#include "dma.h"
#include "simd.h"
#include <stdio.h>
#include <memory.h>
#include <math.h>

#ifndef __thread_local
#define __thread_local const //本地IDE美化
#endif

#ifndef realt
#define realt float
#endif


__thread_local int my_core, core_x, core_y;
__thread_local struct lbm_init_param param;                 // 28 * (4B / 8B) + 4B
__thread_local volatile long local_flag[FLAG_SIZE];          // 16 * 8B
__thread_local volatile long out_flag[FLAG_SIZE];          // 16 * 8B
__thread_local long slave_flag_to_wait;
__thread_local volatile unsigned long get_reply, put_reply; // 2 * 8B
__thread_local volatile unsigned long get_reply_god, put_reply_god; // 2 * 8B
__thread_local volatile unsigned long get_reply_target, put_reply_target; // 2 * 8B
__thread_local volatile intv8 bcast_tmp256;                 // 8B

void* nw_ldm_malloc(int sz)//will not see FBI WARNING during compile
{
    long ret = ldm_malloc(sz);
    return (void*)ret;
}

//conghe profile module (￣△￣；)
__thread_local long pf_stcc_, pf_edcc_;
inline void stpf()
{
    RPCC(pf_stcc_);
}

inline void edpf(int n)
{
    RPCC(pf_edcc_);
    out_flag[n] += pf_edcc_ - pf_stcc_;
    pf_stcc_ = pf_edcc_;
}

//align
void* _nice_ptr_head(char* ptr, int alignment, int right_shft)
{
    long ptrl = (long)ptr;
    int rem = ptrl % alignment;
    if(rem > 0)
    {
        return ptr + (alignment - rem + right_shft);
    }
    else
    {
        return ptr + right_shft;
    }
}

inline static void wait_all_async_write_data()
{
    while(put_reply_god != put_reply_target);
}

inline static void wait_all_async_get_data()
{
    while(get_reply_god != get_reply_target);
}

void cpe0_bcast_256(int my_core, volatile intv8* bcast_data)
{
    bcast_tmp256 = *bcast_data;

    if (my_core / 8 == 0)
    {
        __builtin_sw64_putC(bcast_tmp256, 0x0000000F);
    } else {
        bcast_tmp256 = __builtin_sw64_getxC(bcast_tmp256);
    }

    if (my_core % 8 == 0)
    {
        asm volatile ("nop":::"memory");
        __builtin_sw64_putR(bcast_tmp256, 0x0000000F);
        *bcast_data = bcast_tmp256;
    } else {
        asm volatile ("nop":::"memory");
        *bcast_data = __builtin_sw64_getxR(bcast_tmp256);
    }
}

void standard_wait_flag(int my_core)
{
    if(my_core == 0)
    {
        //printf("SPE waiting for %ld\n", slave_flag_to_wait);

        while(1)
        {
            get_reply = 0;
            athread_get(
                PE_MODE, param.host_flag, (void*)&local_flag[0],
                sizeof(long) * FLAG_SIZE, (void*)(&get_reply), 0, 0, 0
            );
            while(get_reply != 1);
            asm volatile ("#nop":::"memory");
            if(local_flag[0] >= slave_flag_to_wait)
                break;
        }

        slave_flag_to_wait++;
    }

    int i;
    for(i = 0; i < FLAG_SIZE / 4; i++)
    {
        cpe0_bcast_256(my_core, (intv8*)&local_flag[i * 4]);
    }
}

void standard_write_flag(int my_core)
{
    out_flag[0] = out_flag[0] + 1;

    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);//this is necessary, thanks to Luo Yangze
    // athread_dma_barrier();
    if(my_core == 0)
    {
        put_reply = 0;
        athread_put(
            PE_MODE, (void*)&out_flag[0], param.slave_flag,
            sizeof(long) * FLAG_SIZE, (void*)(&put_reply), 0, 0
        );
        while(put_reply != 1);

    }
    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
}

void seq_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size)
{
    put_reply = 0;
    athread_put(
        PE_MODE, cpe_data_ptr, mpe_data_ptr,
        size, (void*)(&put_reply), 0, 0
    );
    while (put_reply != 1);
}

void seq_read_data(void* cpe_data_ptr, void* mpe_data_ptr, int size)
{
    get_reply = 0;
    athread_get(
                PE_MODE, mpe_data_ptr, cpe_data_ptr,
                size, (void*)(&get_reply), 0, 0, 0
            );
    while (get_reply != 1);
}

//需要异步读取数据
void async_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size)
{
    athread_put(
        PE_MODE, cpe_data_ptr, mpe_data_ptr,
        size, (void*)(&put_reply_god), 0, 0
    );
    put_reply_target++;
}

void async_read_data(void* cpe_data_ptr, void* mpe_data_ptr, int size)
{
    athread_get(
        PE_MODE, mpe_data_ptr, cpe_data_ptr,
        size, (void*)(&get_reply_god), 0, 0, 0
    );
    get_reply_target++;
}

//std area

//const data
#define S_FLUID    0
#define S_NOSLIP   1
#define S_VELOCITY 2
#define S_BOUNCE   3
#define S_PRESSURE 4

//data area
__thread_local realt s_w[19] __attribute__((aligned(32))) = {
        (1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };
__thread_local int s_e_x[19] __attribute__((aligned(32))) = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0};
__thread_local int s_e_y[19] __attribute__((aligned(32))) = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0};
__thread_local int s_e_z[19] __attribute__((aligned(32)))  = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0};
__thread_local realt s_e_T[19][4] __attribute__((aligned(32))) = {
	{0.,1.,0.,0}, {0.,-1.,0.,0.}, {1.,0.,0.,0.}, {-1.,0.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,-1.,0.},
	{1.,1.,0.,0.}, {-1.,1.,0.,0.}, {1.,-1.,0.,0.}, {-1.,-1.,0.,0.}, {0.,1.,1.,0.}, {0.,1.,-1.,0.},
	{0.,-1.,1.,0.}, {0.,-1.,-1.,0.}, {1.,0.,1.,0.}, {1.,0.,-1.,0.}, {-1.,0.,1.,0.}, {-1.,0.,-1.,0.}, {0.,0.,0.,0.}
};
__thread_local int s_dfInv[19] __attribute__((aligned(32))) = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};


void s_stream_insane_step(realt s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], realt s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE])
{
    int inv, l;

    int local_x, local_y;

    int sis;
#define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(local_x = 0; local_x < 2; local_x++)
    {
        for(local_y = 0; local_y < 2; local_y++)
        {
            for(sis = 0; sis < INSANE_SIZE; sis++)
            {
                if(s_insane_flags[local_x][local_y][sis] == S_FLUID) {
                    for(l = 0; l < 19; l++) {
                        inv = s_dfInv[l];
                        s_insane_current_line[local_x][local_y][sis * 19 + l] = s_insane_other_line[local_x + s_e_x[inv] + 1][local_y + s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                    }
                }
                if(s_insane_flags[local_x][local_y][sis] == S_BOUNCE) {

                    for(l = 0; l < 19; l++) {
                        inv = s_dfInv[l];
                        if(s_insane_walls[local_x][local_y][sis * 19 + l]) {
                            s_insane_current_line[local_x][local_y][sis * 19 + l] = s_insane_other_line[local_x + 1][local_y + 1][other_step_index(sis, inv)];
                        } else
                        {
                            s_insane_current_line[local_x][local_y][sis * 19 + l] = s_insane_other_line[local_x + s_e_x[inv] + 1][local_y + s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                        }
                    }
                }
            }
        }
    }
}

void s_collide_insane_step(realt s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], realt s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE])
{

    int local_x, local_y;

    int sis, l;

    realt rho, u_x, u_y, u_z;
    realt fi;
    realt feq[19], Qo, omegaNew, Sij[12];

#define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))

    for(local_x = 0; local_x < 2; local_x++)
    {
        for(local_y = 0; local_y < 2; local_y++)
        {
            for(sis = 0; sis < INSANE_SIZE; sis++)
            {
                if(s_insane_flags[local_x][local_y][sis] == S_FLUID || s_insane_flags[local_x][local_y][sis] == S_BOUNCE)
                {
                    floatv4 U = simd_set_floatv4(0., 0., 0., 0.), E;
                    #pragma unroll[19]
                    for (l = 0; l < 19; l++) {
                        // fi = a_current_line[l];
                        fi = s_insane_current_line[local_x][local_y][sis * 19 + l];
                        E = simd_set_floatv4(1., s_e_T[l][0], s_e_T[l][1], s_e_T[l][2]);
                        U = fi * E + U;
                    }
                    simd_store(U, &(Sij[0]));
                    rho = Sij[0], u_x = Sij[1], u_y = Sij[2], u_z = Sij[3];

                    u_x /= rho;
                    u_y /= rho;
                    u_z /= rho;

                    for (l = 0; l < 19; l++) {
                        const realt tmp = (s_e_x[l] * u_x + s_e_y[l] * u_y + s_e_z[l] * u_z);
                        feq[l] = s_w[l] * rho * (1.0 -
                                               (3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
                                               (3.0 *     tmp) +
                                               (9.0/2.0 * tmp * tmp));
                    }

                    realt S;
                    floatv4 vs1, vs2, vs3, vs4, vs5 ,ve, vk1, vk2, vk3, vk4, vk5;
                    int x1, y1, k1;
                    
                    vs1 = simd_set_floatv4(0., 0., 0., 0.);
                    vs2 = simd_set_floatv4(0., 0., 0., 0.);
                    vs3 = simd_set_floatv4(0., 0., 0., 0.);
                    #pragma unroll[19]
                    for(k1 = 0; k1 < 19; k1++) {
                        float tmp = (s_insane_current_line[local_x][local_y][sis * 19 + k1] - feq[k1]);
                        // float tmp = (a_current_line[k1] - feq[k1]);
                        simd_load(ve, &(s_e_T[k1][0]));
                        vk1 = tmp * s_e_T[k1][0];
                        vk2 = tmp * s_e_T[k1][1];
                        vk3 = tmp * s_e_T[k1][2];
                        
                        vs1 = vk1 * ve + vs1;
                        vs2 = vk2 * ve + vs2;
                        vs3 = vk3 * ve + vs3;
                    }
                    simd_store(vs1, &(Sij[0]));
                    simd_store(vs2, &(Sij[4]));
                    simd_store(vs3, &(Sij[8]));

                    //用simd乘加会丢精度，也不会变快
                    Qo=0;
                    #pragma unroll[12]
                    for(x1 = 0; x1 < 12; ++x1)
                    {
                        Qo += Sij[x1] * Sij[x1];
                    }

                    param.nu = (2.0 / param.omega - 1.0) / 6.0;
                    S = (-param.nu + sqrt(param.nu * param.nu + 18 * param.CSmago * param.CSmago * sqrt(Qo))) / 6.0 / param.CSmago / param.CSmago;
                    omegaNew = 1.0 / (3.0 * (param.nu + param.CSmago * param.CSmago * S) + 0.5);

                    
                    
                    for (l = 0; l < 19; l++) {
                        s_insane_current_line[local_x][local_y][sis * 19 + l] =
                                (1.0 - omegaNew) * s_insane_current_line[local_x][local_y][sis * 19 + l] +
                                omegaNew * feq[l];
                    }
                }
            }
        }
    }
}

void s_stream_std(realt s_std_other_line[3][3][(STEP_SIZE + 2) * 19], realt s_std_current_line[STEP_SIZE * 19], int s_std_walls[STEP_SIZE * 19], int s_std_flags[STEP_SIZE])
{
    int inv, l;

    int sis;
    #define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
        if(s_std_flags[sis] == S_FLUID) {
            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                s_std_current_line[sis * 19 + l] = s_std_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
            }
        }
        if(s_std_flags[sis] == S_BOUNCE) {

            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                if(s_std_walls[sis * 19 + l]) {
                    s_std_current_line[sis * 19 + l] = s_std_other_line[1][1][other_step_index(sis, inv)];
                } else 
                {
                    s_std_current_line[sis * 19 + l] = s_std_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                }
            }
        }
    }
}

void s_collide_std(realt s_std_other_line[3][3][(STEP_SIZE + 2) * 19], realt s_std_current_line[STEP_SIZE * 19], int s_std_walls[STEP_SIZE * 19], int s_std_flags[STEP_SIZE])
{
    int l;

    realt rho, u_x, u_y, u_z;
    realt fi;
    realt feq[19], Qo, omegaNew, Sij[12];

    int sis;
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
        if(s_std_flags[sis] == S_FLUID || s_std_flags[sis] == S_BOUNCE)
        {
            floatv4 U = simd_set_floatv4(0., 0., 0., 0.), E;
            #pragma unroll[19]
            for (l = 0; l < 19; l++) {
                // fi = a_current_line[l];
                fi = s_std_current_line[sis * 19 + l];
                E = simd_set_floatv4(1., s_e_T[l][0], s_e_T[l][1], s_e_T[l][2]);
                U = fi * E + U;
            }
            simd_store(U, &(Sij[0]));
            rho = Sij[0], u_x = Sij[1], u_y = Sij[2], u_z = Sij[3];

            u_x /= rho;
            u_y /= rho;
            u_z /= rho;

            for (l = 0; l < 19; l++) {
                const realt tmp = (s_e_x[l] * u_x + s_e_y[l] * u_y + s_e_z[l] * u_z);
                feq[l] = s_w[l] * rho * (1.0 -
                                       (3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
                                       (3.0 *     tmp) +
                                       (9.0/2.0 * tmp * tmp));
            }

            realt S;
			floatv4 vs1, vs2, vs3, vs4, vs5 ,ve, vk1, vk2, vk3, vk4, vk5;
            int x1, y1, k1;
			
			vs1 = simd_set_floatv4(0., 0., 0., 0.);
			vs2 = simd_set_floatv4(0., 0., 0., 0.);
			vs3 = simd_set_floatv4(0., 0., 0., 0.);
            #pragma unroll[19]
			for(k1 = 0; k1 < 19; k1++) {
				float tmp = (s_std_current_line[sis * 19 + k1] - feq[k1]);
				// float tmp = (a_current_line[k1] - feq[k1]);
				simd_load(ve, &(s_e_T[k1][0]));
				vk1 = tmp * s_e_T[k1][0];
				vk2 = tmp * s_e_T[k1][1];
				vk3 = tmp * s_e_T[k1][2];
				
				vs1 = vk1 * ve + vs1;
				vs2 = vk2 * ve + vs2;
				vs3 = vk3 * ve + vs3;
            }
			simd_store(vs1, &(Sij[0]));
			simd_store(vs2, &(Sij[4]));
			simd_store(vs3, &(Sij[8]));

            //用simd乘加会丢精度，也不会变快
            Qo=0;
            #pragma unroll[12]
			for(x1 = 0; x1 < 12; ++x1)
			{
				Qo += Sij[x1] * Sij[x1];
			}


            param.nu = (2.0 / param.omega - 1.0) / 6.0;
            S = (-param.nu + sqrt(param.nu * param.nu + 18 * param.CSmago * param.CSmago * sqrt(Qo))) / 6.0 / param.CSmago / param.CSmago;
            omegaNew = 1.0 / (3.0 * (param.nu + param.CSmago * param.CSmago * S) + 0.5);


            for (l = 0; l < 19; l++) {
               s_std_current_line[sis * 19 + l] =
                       (1.0 - omegaNew) * s_std_current_line[sis * 19 + l] +
                       omegaNew * feq[l];
            }
        }
    }
}

void stream_and_collide_xy_insane(int nz,
                           realt* current_line,
                           int* walls,
                           int* flags,
                           realt* other_line,
                           realt s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], realt s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE])
{
    int k, l;

    int x_node_stride = (param.y_sec + 2) * param.Z * 19;
    int y_node_stride = param.Z * 19;
    int x_wall_stride = (param.y_sec) * param.Z * 19;
    int y_wall_stride = param.Z * 19;
    int x_flag_stride = (param.y_sec + 2) * param.Z;
    int y_flag_stride = param.Z;


    int total_step = nz / INSANE_SIZE;

    int cur_step;
    for(cur_step = 0; cur_step < total_step; cur_step++)
    {
        k = cur_step * INSANE_SIZE;
        int z;

        stpf();
        int ex, ey;
        for(ex = 0; ex < 2; ex++)
        {
            for(ey = 0; ey < 2; ey++)
            {
                async_read_data(s_insane_flags[ex][ey], flags + (k + ex * x_flag_stride + ey * y_flag_stride), sizeof(int) * INSANE_SIZE);
                async_read_data(s_insane_walls[ex][ey], walls + (k * 19 + ex * x_wall_stride + ey * y_wall_stride), sizeof(int) * INSANE_SIZE * 19);
                async_read_data(s_insane_current_line[ex][ey], current_line + (k * 19 + ex * x_node_stride + ey * y_node_stride), sizeof(realt) * INSANE_SIZE * 19);
            }
        }
        edpf(CPF1_INSANE_READ1);

        stpf();
        k = cur_step * INSANE_SIZE;
        int z_st = k - 1, z_ed = k + INSANE_SIZE;
        if(k == 0)
        {
            z_st = 0;
        }
        if(z_ed == nz)
        {
            z_ed--;
        }

        for(ex = -1; ex <= 2; ex++)
        {
            for(ey = -1; ey <= 2; ey++)
            {
                async_read_data(s_insane_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19,
                        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19],
                        sizeof(realt) * (z_ed - z_st + 1) * 19);
            }
        }
        wait_all_async_get_data();
        edpf(CPF1_INSANE_READ2);


        stpf();
        s_stream_insane_step(s_insane_other_line, s_insane_current_line, s_insane_walls, s_insane_flags);
        edpf(CPF1_INSANE_STREAM);

        stpf();
        s_collide_insane_step(s_insane_other_line, s_insane_current_line, s_insane_walls, s_insane_flags);
        edpf(CPF1_INSANE_COLLIDE);


        stpf();
        for(ex = 0; ex < 2; ex++)
        {
            for(ey = 0; ey < 2; ey++)
            {
                async_write_data(current_line + (k * 19 + ex * x_node_stride + ey * y_node_stride), s_insane_current_line[ex][ey], sizeof(realt) * INSANE_SIZE * 19);
            }
        }
        wait_all_async_write_data();
        edpf(CPF1_INSANE_WRITE);

    }

}

__thread_local int dbg_print_cnt;

__thread_local int dbg_x, dbg_y;

void s_stream_and_collide_xy(int nz,
                           realt* current_line,
                           int* walls,
                           int* flags,
                           realt* other_line,
                           realt s_std_other_line[3][3][(STEP_SIZE + 2) * 19], realt s_std_current_line[STEP_SIZE * 19], int s_std_walls[STEP_SIZE * 19], int s_std_flags[STEP_SIZE])
{
    int k, l;

    int x_node_stride = (param.y_sec + 2) * param.Z * 19;
    int y_node_stride = param.Z * 19;

    int total_step = param.Z / STEP_SIZE;


    int cur_step;
    for(cur_step = 0; cur_step < total_step; cur_step++)
    {
        k = cur_step * STEP_SIZE;
        int z;

        stpf();
        async_read_data(s_std_flags, flags + k, sizeof(int) * STEP_SIZE);
        async_read_data(s_std_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        async_read_data(s_std_current_line, current_line + k * 19, sizeof(realt) * STEP_SIZE * 19);
        edpf(CPF1_STD_READ1);

        // memcpy(s_std_flags, flags + k, sizeof(int) * STEP_SIZE);
        // memcpy(s_std_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        // memcpy(s_std_current_line, current_line + k * 19, sizeof(int) * STEP_SIZE * 19);

        k = cur_step * STEP_SIZE;
        int z_st = k - 1, z_ed = k + STEP_SIZE;
        if(k == 0)
        {
            z_st = 0;
        }
        if(z_ed == nz)
        {
            z_ed--;
        }

        stpf();
        int ex, ey;
        for(ex = -1; ex <= 1; ex++)
        {
            for(ey = -1; ey <= 1; ey++)
            {
                async_read_data(s_std_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                       &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                       sizeof(realt) * (z_ed - z_st + 1) * 19);
                // memcpy(s_std_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                //        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                //        sizeof(realt) * (z_ed - z_st + 1) * 19);
            }
        }
        wait_all_async_get_data();
        edpf(CPF1_STD_READ2);

        stpf();
        s_stream_std(s_std_other_line, s_std_current_line, s_std_walls, s_std_flags);
        edpf(CPF1_STD_STREAM);
        
        stpf();
        s_collide_std(s_std_other_line, s_std_current_line, s_std_walls, s_std_flags);
        edpf(CPF1_STD_COLLIDE);

        stpf();
        seq_write_data(current_line + k * 19, s_std_current_line, sizeof(realt) * STEP_SIZE * 19);
        edpf(CPF1_STD_WRITE);
        // memcpy(current_line + k * 19, s_std_current_line, sizeof(realt) * STEP_SIZE * 19);
    }

}

#define s_nodes_idx(d2, d3, d4, d5) ((((d2) * (param.y_sec + 2) + (d3)) * param.Z + (d4)) * 19 + (d5))
#define s_flags_idx(d2, d3, d4) (((d2) * (param.y_sec + 2) + (d3)) * param.Z + (d4))
#define s_walls_idx(d2, d3, d4, d5) ((((d2) * param.y_sec + (d3)) * param.Z + (d4)) * 19 + (d5))

void insane_lbm()
{
    realt s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19];//in
    realt s_insane_current_line[2][2][INSANE_SIZE * 19];
    int s_insane_walls[2][2][INSANE_SIZE * 19];
    int s_insane_flags[2][2][INSANE_SIZE];
    // realt s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], realt s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE]

    int  insane_point_cnt = local_flag[INSANE_POINT_CNT];
    int* insane_xs_ptr = (int*)local_flag[INSANE_XS_PTR];
    int* insane_ys_ptr = (int*)local_flag[INSANE_YS_PTR];
    realt* current_head = (realt*)local_flag[STD_CURRENT_HEAD];
    realt* other_head = (realt*)local_flag[STD_OTHER_HEAD];

    int insane_my_start, insane_my_part, insane_max_part;
    insane_my_start = BLOCK_LOW(my_core, 64, insane_point_cnt);
    insane_my_part = BLOCK_SIZE(my_core, 64, insane_point_cnt);
    insane_max_part = BLOCK_SIZE(0, 64, insane_point_cnt);

    int* s_insane_xs = insane_xs_ptr;
    int* s_insane_ys = insane_ys_ptr;

    int i, j;
    for(i = insane_my_start; i < insane_my_start + insane_my_part; i++)
    {
        int x = s_insane_xs[i];//GLD!!!
        int y = s_insane_ys[i];//GLD!!!

        realt* point_current = current_head + s_nodes_idx(x, y, 0, 0);
        realt* point_other = other_head + s_nodes_idx(x, y, 0, 0);
        int* point_walls = param.walls + s_walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = param.flags + s_flags_idx(x, y, 0);

        stream_and_collide_xy_insane(param.Z, point_current, point_walls,
                                     point_flags, point_other,
                                     s_insane_other_line, s_insane_current_line, s_insane_walls, s_insane_flags);
    }
}

void std_lbm()
{
    realt s_std_other_line[3][3][(STEP_SIZE + 2) * 19];//in
    realt s_std_current_line[STEP_SIZE * 19];
    int s_std_walls[STEP_SIZE * 19];
    int s_std_flags[STEP_SIZE];
    //realt s_std_other_line[3][3][(STEP_SIZE + 2) * 19], realt s_std_current_line[STEP_SIZE * 19], int s_std_walls[STEP_SIZE * 19], int s_std_flags[STEP_SIZE]

    dbg_print_cnt = 0;
    int point_cnt = local_flag[STD_POINT_CNT];
    int* xs_ptr = (int*)local_flag[STD_XS_PTR];
    int* ys_ptr = (int*)local_flag[STD_YS_PTR];
    realt* current_head = (realt*)local_flag[STD_CURRENT_HEAD];
    realt* other_head = (realt*)local_flag[STD_OTHER_HEAD];

    int my_start, my_part, max_part;
    my_start = BLOCK_LOW(my_core, 64, point_cnt);
    my_part = BLOCK_SIZE(my_core, 64, point_cnt);
    max_part = BLOCK_SIZE(0, 64, point_cnt);

    int* s_xs = xs_ptr;
    int* s_ys = ys_ptr;

    int i, j;
    for(i = my_start; i < my_start + my_part; i++)
    {
        int x = s_xs[i];//GLD!!!
        int y = s_ys[i];//GLD!!!
        dbg_x = x;
        dbg_y = y;

        //var for this point
        realt* point_current = current_head + s_nodes_idx(x, y, 0, 0);
        realt* point_other = other_head + s_nodes_idx(x, y, 0, 0);
        int* point_walls = param.walls + s_walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = param.flags + s_flags_idx(x, y, 0);

        s_stream_and_collide_xy(param.Z, point_current, point_walls, 
                              point_flags, point_other,
                              s_std_other_line, s_std_current_line, s_std_walls, s_std_flags);
    }

    // standard_write_flag(my_core);//removed!
}

void cpe_athread_daemon(void *_param)
{
    my_core = athread_get_id(-1);

    int i, j, k;

    param = *((struct lbm_init_param*) _param);

    core_x = my_core % 8;
    core_y = my_core / 8;

    if(my_core == 0 && param.my_id == 0)
    {
        printf("@@@ MESSAGE FROM CPE @@@\n");
        printf("iter = %d, x_sec = %d, y_sec = %d, Z = %d\n", param.iter, param.x_sec, param.y_sec, param.Z);
        printf("val: nu = %f, omega = %f, CSmago = %f\n", param.nu, param.omega, param.CSmago);
        printf("ptrs: %x, %x, %x, %x\n",
               param.host_flag, param.slave_flag, param.walls, param.flags);
        printf("sizeof(int) = %d, sizeof(realt) = %d\n", sizeof(int), sizeof(realt));
        printf("@@@ MESSAGE FROM CPE @@@\n");
    }

    for (i = 0; i < FLAG_SIZE; i++)
    {
        local_flag[i] = 0; out_flag[i] = 0;
    }
    slave_flag_to_wait = 1;

    get_reply_god = 0, put_reply_god = 0;//the Jesus do not like waiting for someone
    get_reply_target = 0, put_reply_target = 0;

    //initial arr

    if(my_core == 0 && param.my_id == 0)
    {
        printf("CPE ENTER LOOP\n");
    }
    while(1)
    {
        stpf();
        standard_wait_flag(my_core);
        edpf(CPF1_WAIT);

        if (local_flag[KERNEL_ACTION] == EXIT_FLAG)
        {
            standard_write_flag(my_core);
            standard_write_flag(my_core);
            break;
        }
        if (local_flag[KERNEL_ACTION] == HANDSHAKE_FLAG)
        {
            if(my_core == 0 && param.my_id == 5)
            {
                printf("@@@ CPE HANDSHAKE_FLAG @@@\n");
                fflush(stdout);
            }
            standard_write_flag(my_core);
        }
        else if (local_flag[KERNEL_ACTION] == STD_LBM_FLAG)
        {
            if(my_core == 0 && param.my_id == 5)
            {
                printf("@@@ CPE STD_LBM_FLAG @@@\n");
                fflush(stdout);
            }
            std_lbm();
            insane_lbm();
            standard_write_flag(my_core);
        }
        // else if(local_flag[KERNEL_ACTION] == INSANE_LBM_FLAG)
        // {
        //     insane_lbm();
        // }
        else
        {
            if(my_core == 0 && param.my_id == 0)
            {
                printf("unknown area!!!!!!!!!\n");
            }
        }
    }
}

