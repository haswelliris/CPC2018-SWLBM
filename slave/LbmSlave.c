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
__thread_local realt s_step_other_line[3][3][(STEP_SIZE + 2) * 19] __attribute__((aligned(32)));//in
__thread_local realt s_step_current_line[STEP_SIZE * 19] __attribute__((aligned(32)));
__thread_local int s_step_walls[STEP_SIZE * 19] __attribute__((aligned(32)));
__thread_local int s_step_flags[STEP_SIZE] __attribute__((aligned(32)));
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
__thread_local realt  Sij[12] __attribute__((aligned(32)));
__thread_local realt  feq[20] __attribute__((aligned(32)));
__thread_local floatv4  simd_feq[19] __attribute__((aligned(32)));
__thread_local floatv4  simd_Sij[9] __attribute__((aligned(32)));
// __thread_local realt  a_current_line[20] __attribute__((aligned(32)));

void s_stream_step()
{
    int inv, l;
	int sis;
    #define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
	//stream
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
		
		// for(l = 0; l < 19; l++) {
			// a_current_line[l] = s_step_current_line[sis * 19 + l];
		// }
        if(s_step_flags[sis] == S_FLUID) {
            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
				// a_current_line[l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                s_step_current_line[sis * 19 + l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
            }
        }
        if(s_step_flags[sis] == S_BOUNCE) {

            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                if(s_step_walls[sis * 19 + l]) {
					// a_current_line[l] = s_step_other_line[1][1][other_step_index(sis, inv)];
                    s_step_current_line[sis * 19 + l] = s_step_other_line[1][1][other_step_index(sis, inv)];
                } else 
                {
					// a_current_line[l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                    s_step_current_line[sis * 19 + l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                }
            }
        }
	}
}

//半SIMD化
void s_collide_step(int st, int ed)
{
	realt rho, u_x, u_y, u_z;
    realt fi;
    realt Qo, omegaNew;
	
	int sis, l;
	//coli
	for(sis = st; sis < ed; sis++)
    {
        if(s_step_flags[sis] == S_FLUID || s_step_flags[sis] == S_BOUNCE)
        {
            //rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
			floatv4 U = simd_set_floatv4(0., 0., 0., 0.), E;
            for (l = 0; l < 19; l++) {
				// fi = a_current_line[l];
				fi = s_step_current_line[sis * 19 + l];
				E = simd_set_floatv4(1., s_e_T[l][0], s_e_T[l][1], s_e_T[l][2]);
				U = fi * E + U;
                //rho += fi;
                //u_x += s_e_x[l] * fi;
                //u_y += s_e_y[l] * fi;
                //u_z += s_e_z[l] * fi;
            }
			simd_store(U, &(Sij[0]));
			rho = Sij[0];
            u_x = Sij[1];
            u_y = Sij[2];
            u_z = Sij[3];
			
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

            Qo=0;
            realt S;
			floatv4 vs1, vs2, vs3, vs4, vs5 ,ve, vk1, vk2, vk3, vk4, vk5;
            int x1, y1, k1;
			
			//for(x1 = 0; x1 < 12; ++x1)
			//{
			//	Sij[x1] = 0;
			//}
			vs1 = simd_set_floatv4(0., 0., 0., 0.);
			vs2 = simd_set_floatv4(0., 0., 0., 0.);
			vs3 = simd_set_floatv4(0., 0., 0., 0.);
			for(k1 = 0; k1 < 19; k1++) {
				float tmp = (s_step_current_line[sis * 19 + k1] - feq[k1]);
				// float tmp = (a_current_line[k1] - feq[k1]);
				simd_load(ve, &(s_e_T[k1][0]));
				vk1 = tmp * s_e_T[k1][0];
				vk2 = tmp * s_e_T[k1][1];
				vk3 = tmp * s_e_T[k1][2];
				
				vs1 = vk1 * ve + vs1;
				vs2 = vk2 * ve + vs2;
				vs3 = vk3 * ve + vs3;
				
				//for(x1 = 0; x1 < 3; x1++) {
				//	simd_load(vsij, &(Sij[x1 * 4]));
				//	vsij = (tmp * s_e_T[k1][x1]) * ve + vsij;
				//	simd_store(vsij, &(Sij[x1 * 4]));
					//Sij[x1 * 3 + 0] += s_e_T[k1][x1] * s_e_T[k1][0] * tmp;
					//Sij[x1 * 3 + 1] += s_e_T[k1][x1] * s_e_T[k1][1] * tmp;
					//Sij[x1 * 3 + 2] += s_e_T[k1][x1] * s_e_T[k1][2] * tmp;
                //}
            }
			simd_store(vs1, &(Sij[0]));
			simd_store(vs2, &(Sij[4]));
			simd_store(vs3, &(Sij[8]));
			for(x1 = 0; x1 < 12; ++x1)
			{
				Qo += Sij[x1] * Sij[x1];
			}

            param.nu = (2.0 / param.omega - 1.0) / 6.0;
            S = (-param.nu + sqrt(param.nu * param.nu + 18 * param.CSmago * param.CSmago * sqrt(Qo))) / 6.0 / param.CSmago / param.CSmago;
            omegaNew = 1.0 / (3.0 * (param.nu + param.CSmago * param.CSmago * S) + 0.5);


			//精度丢失还只快1s，算了算了
			// a_current_line[19] = 0;
			
			// simd_load(vs1, &(a_current_line[0]));
			// simd_load(vs2, &(a_current_line[4]));
			// simd_load(vs3, &(a_current_line[8]));
			// simd_load(vs4, &(a_current_line[12]));
			// simd_load(vs5, &(a_current_line[16]));
			
			// S = (1.0 - omegaNew);
			
			// simd_load(vk1, &(feq[0]));
			// simd_load(vk2, &(feq[4]));
			// simd_load(vk3, &(feq[8]));
			// simd_load(vk4, &(feq[12]));
			// simd_load(vk5, &(feq[16]));
			
			
			// vs1 = S * vs1;
			// vs2 = S * vs2;
			// vs3 = S * vs3;
			// vs4 = S * vs4;
			// vs5 = S * vs5;
			
			// vs1 = omegaNew * vk1 + vs1;
			// vs2 = omegaNew * vk2 + vs2;
			// vs3 = omegaNew * vk3 + vs3;
			// vs4 = omegaNew * vk4 + vs4;
			// vs5 = omegaNew * vk5 + vs5;
			
			// simd_store(vs1, &(a_current_line[0]));
			// simd_store(vs2, &(a_current_line[4]));
			// simd_store(vs3, &(a_current_line[8]));
			// simd_store(vs4, &(a_current_line[12]));
			// simd_store(vs5, &(a_current_line[16]));
			
			
            for (l = 0; l < 19; l++) {
              s_step_current_line[sis * 19 + l] =
                      (1.0 - omegaNew) * s_step_current_line[sis * 19 + l] +
                      omegaNew * feq[l];
            }
			//for (l = 0; l < 19; l++) {
            //   a_current_line[l] =
            //           (1.0 - omegaNew) * a_current_line[l] +
            //           omegaNew * feq[l];
            //}
			// for (l = 0; l < 19; l++) {
				// s_step_current_line[sis * 19 + l] = a_current_line[l];
			// }
        }
    }
}

//TODO !!!!精度丢失严重，需要重写或者debug
void s_simd_collide_step()
{
	int STEP_SIZE_MAIN = (STEP_SIZE / 4) * 4;
	
	floatv4 rho, u_x, u_y, u_z;
	floatv4 fi;
    floatv4 Qo, omegaNew;
	int sis, l;
	for(sis = 0; sis < STEP_SIZE_MAIN; sis+=4)
    {
		if(1)
		//if(s_step_flags[sis] == S_FLUID || s_step_flags[sis] == S_BOUNCE)
        {
            rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
            for (l = 0; l < 19; l++) {
				fi = simd_set_floatv4(
					s_step_current_line[sis * 19 + l],
					s_step_current_line[(sis + 1) * 19 + l],
					s_step_current_line[(sis + 2) * 19 + l],
					s_step_current_line[(sis + 3) * 19 + l]
				);
                rho = rho + fi;
                u_x = s_e_x[l] * fi + u_x;
                u_y = s_e_y[l] * fi + u_y;
                u_z = s_e_z[l] * fi + u_z;
            }

            u_x = u_x / rho;
            u_y = u_y / rho;
            u_z = u_z / rho;

            for (l = 0; l < 19; l++) {
                floatv4 tmp = (s_e_x[l] * u_x + s_e_y[l] * u_y + s_e_z[l] * u_z);
                simd_feq[l] = s_w[l] * rho * (1.0 -
                                       (3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
                                       (3.0 *     tmp) +
                                       (9.0/2.0 * tmp * tmp));
            }

            Qo=0;
            floatv4 S;
            int x1, y1, k1;
			for(x1 = 0; x1 < 9; x1++) {
				simd_Sij[x1] = 0;
			}
			for(k1 = 0; k1 < 19; k1++) {
				floatv4 tmp = simd_set_floatv4(
					s_step_current_line[sis * 19 + k1],
					s_step_current_line[(sis + 1) * 19 + k1],
					s_step_current_line[(sis + 2) * 19 + k1],
					s_step_current_line[(sis + 3) * 19 + k1]
				);
				for(x1 = 0; x1 < 3; x1++) {
					for(y1 = 0; y1 < 3; y1++) {
                       simd_Sij[x1 * 3 + y1] = s_e_T[k1][x1] * s_e_T[k1][y1] * (tmp - simd_feq[k1]) + simd_Sij[x1 * 3 + y1];
                    }
				}
            }
			for(x1 = 0; x1 < 9; x1++) {
				Qo = simd_Sij[x1] * simd_Sij[x1] + Qo;
			}
			
			Qo = __builtin_sw64_sgqrts(Qo);
			Qo = param.nu * param.nu + 18 * param.CSmago * param.CSmago * Qo;
			Qo = __builtin_sw64_sgqrts(Qo);
            param.nu = (2.0 / param.omega - 1.0) / 6.0;
            S = (-param.nu + Qo) / 6.0 / param.CSmago / param.CSmago;
            omegaNew = 1.0 / (3.0 * (param.nu + param.CSmago * param.CSmago * S) + 0.5);



			
            for (l = 0; l < 19; l++) {
				floatv4 tmp = simd_set_floatv4(
					s_step_current_line[sis * 19 + l],
					s_step_current_line[(sis + 1) * 19 + l],
					s_step_current_line[(sis + 2) * 19 + l],
					s_step_current_line[(sis + 3) * 19 + l]
				);
				simd_feq[l] = (1.0 - omegaNew) * tmp + omegaNew * simd_feq[l];
            }
			
			for (l = 0; l < 19; l++) {
				simd_store(simd_feq[l], &(feq[0]));
				for(x1 = 0; x1 < 4; ++x1)
				{
					if(s_step_flags[sis + x1] == S_FLUID || s_step_flags[sis + x1] == S_BOUNCE)
					{
						s_step_current_line[(sis + x1) * 19 + l] = feq[x1];
					}
				}
			}
        }
	}
	
	s_collide_step(STEP_SIZE_MAIN, STEP_SIZE);
}

__thread_local int dbg_print_cnt;

__thread_local int dbg_x, dbg_y;

void s_stream_and_collide_xy(int nz,
                           realt* current_line,
                           int* walls,
                           int* flags,
                           realt* other_line)
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
        seq_read_data(s_step_flags, flags + k, sizeof(int) * STEP_SIZE);
        seq_read_data(s_step_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        seq_read_data(s_step_current_line, current_line + k * 19, sizeof(realt) * STEP_SIZE * 19);
        edpf(CPF1_READ1);

        // memcpy(s_step_flags, flags + k, sizeof(int) * STEP_SIZE);
        // memcpy(s_step_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        // memcpy(s_step_current_line, current_line + k * 19, sizeof(int) * STEP_SIZE * 19);

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
                seq_read_data(s_step_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                       &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                       sizeof(realt) * (z_ed - z_st + 1) * 19);
                // memcpy(s_step_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                //        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                //        sizeof(realt) * (z_ed - z_st + 1) * 19);
            }
        }
        edpf(CPF1_READ2);

        stpf();
        s_stream_step();
        edpf(CPF1_STREAM);
        
        stpf();
        s_collide_step(0, STEP_SIZE);
        edpf(CPF1_COLLIDE);

        stpf();
        seq_write_data(current_line + k * 19, s_step_current_line, sizeof(realt) * STEP_SIZE * 19);
        edpf(CPF1_WRITE);
        // memcpy(current_line + k * 19, s_step_current_line, sizeof(realt) * STEP_SIZE * 19);
    }

}

void std_lbm()
{
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

    int* fk_xs = xs_ptr;
    int* fk_ys = ys_ptr;

#define s_nodes_idx(d2, d3, d4, d5) ((((d2) * (param.y_sec + 2) + (d3)) * param.Z + (d4)) * 19 + (d5))
#define s_flags_idx(d2, d3, d4) (((d2) * (param.y_sec + 2) + (d3)) * param.Z + (d4))
#define s_walls_idx(d2, d3, d4, d5) ((((d2) * param.y_sec + (d3)) * param.Z + (d4)) * 19 + (d5))

    int i, j;
    for(i = my_start; i < my_start + my_part; i++)
    {
        int x = fk_xs[i];//GLD!!!
        int y = fk_ys[i];//GLD!!!
        dbg_x = x;
        dbg_y = y;

        //var for this point
        realt* point_current = current_head + s_nodes_idx(x, y, 0, 0);
        realt* point_other = other_head + s_nodes_idx(x, y, 0, 0);
        int* point_walls = param.walls + s_walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = param.flags + s_flags_idx(x, y, 0);

        s_stream_and_collide_xy(param.Z, point_current, point_walls, 
                              point_flags, point_other);
    }

    standard_write_flag(my_core);
}

void insane_lbm()
{

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
