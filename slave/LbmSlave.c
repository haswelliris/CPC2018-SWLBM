#include "../header/LbmSlave.h"

#include <slave.h>

#include "dma.h"
#include "simd.h"
#include <stdio.h>
#include <memory.h>

#include "../header/Array.h"

#ifndef __thread_local
#define __thread_local const //本地IDE美化(__thread_local 是神威从核编译的关键字, 相当于cuda的片上内存)
#endif

__thread_local int slave_id, core_x, core_y; // 从核id，所在行，所在列
// 我觉得是7 * 4B + 4B
__thread_local struct lbm_param param;                 // 28 * (4B / 8B) + 4B
__thread_local volatile long local_flag[FLAG_SIZE];          // 16 * 8B
__thread_local volatile long out_flag[FLAG_SIZE];          // 16 * 8B
__thread_local long slave_flag_to_wait;
__thread_local volatile unsigned long get_reply, put_reply; // 2 * 8B
__thread_local volatile unsigned long get_reply_god, put_reply_god; // 2 * 8B
__thread_local volatile unsigned long get_reply_target, put_reply_target; // 2 * 8B
__thread_local volatile intv8 bcast_tmp256;                 // 8 *4B

void cpe0_bcast_256(int slave_id, volatile intv8 *bcast_data);
void standard_wait_flag(int slave_id);
void standard_write_flag(int slave_id);
void std_lbm();
void insane_lbm();

// 守护进程函数，主要逻辑步骤都在这儿
void cpe_athread_daemon(void *_param)
{
    int i, j, k;
    // 从核id
    param = *((struct lbm_param *)_param);
    slave_id = athread_get_id(-1);
    core_x = slave_id % 8;
    core_y = slave_id / 8;

    if (slave_id == 0 && param.master_id == 0)
    {
        printf("\n@@@ MESSAGE FROM CPE 0 @@@\n");
        printf("iter = %d, x_sec = %d, y_sec = %d, z_sec = %d\n", param.iter, param.x_sec, param.y_sec, param.Z);
        printf("val: nu = %f, omega = %f, CSmago = %f\n", param.nu, param.omega, param.CSmago);
        printf("ptrs: %x, %x, %x, %x\n",
               param.host_flag, param.slave_flag, param.walls, param.flags);
        printf("sizeof(int) = %d, sizeof(Real) = %d\n", sizeof(int), sizeof(Real));
        printf("\n");
    }
    // 
    for (i = 0; i < FLAG_SIZE; i++)
    {
        local_flag[i] = 0;
        out_flag[i] = 0;
    }
    slave_flag_to_wait = 1;

    get_reply_god = 0, put_reply_god = 0; //the Jesus do not like waiting for someone
    get_reply_target = 0, put_reply_target = 0;

    if (slave_id == 0 && param.master_id == 0)
    {
        printf("\n@@@ MESSAGE FROM CPE 0 @@@\n");
        printf("CPE ENTER LOOP\n");
        printf("\n");
    }

    while (1)
    {
        standard_wait_flag(slave_id);

        if (local_flag[KERNEL_ACTION] == EXIT_FLAG)
        {
            standard_write_flag(slave_id);
            standard_write_flag(slave_id);
            break;
        }
        // 每一个从核函数，加在下方
                if (local_flag[KERNEL_ACTION] == HANDSHAKE_FLAG)
        {
            standard_write_flag(slave_id);
        }
        else if (local_flag[KERNEL_ACTION] == STD_LBM_FLAG)
        {            if(slave_id == 0 && param.master_id == 5)
            {
                printf("@@@ CPE STD_LBM_FLAG @@@\n");
            }
            std_lbm();
        }
        // else if (local_flag[KERNEL_ACTION] == INSANE_LBM_FLAG)
        // {
        //     insane_lbm();
        // }
        else
        {
            // todo
        }
    }
}

void cpe0_bcast_256(int slave_id, volatile intv8* bcast_data)
{
    bcast_tmp256 = *bcast_data;

    if (slave_id / 8 == 0)
    {
        __builtin_sw64_putC(bcast_tmp256, 0x0000000F);
    } else {
        bcast_tmp256 = __builtin_sw64_getxC(bcast_tmp256);
    }

    if (slave_id % 8 == 0)
    {
        asm volatile ("nop":::"memory");
        __builtin_sw64_putR(bcast_tmp256, 0x0000000F);
        *bcast_data = bcast_tmp256;
    } else {
        asm volatile ("nop":::"memory");
        *bcast_data = __builtin_sw64_getxR(bcast_tmp256);
    }
}

void standard_wait_flag(int slave_id)
{
    if(slave_id == 0)
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
            //内存屏障
            asm volatile ("nop":::"memory");
            // 这个不是很懂
            if(local_flag[0] >= slave_flag_to_wait)
                break;
        }

        slave_flag_to_wait++;
    }
    // 广播16B*8
    int i;
    for(i = 0; i < FLAG_SIZE / 4; i++)
    {
        cpe0_bcast_256(slave_id, (intv8*)&local_flag[i * 4]);
    }
}

void standard_write_flag(int slave_id)
{
    out_flag[0] = out_flag[0] + 1;

    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);//this is necessary, thanks to Luo Yangze
    // athread_dma_barrier();
    if(slave_id == 0)
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


//const data
#define S_FLUID    0
#define S_NOSLIP   1
#define S_VELOCITY 2
#define S_BOUNCE   3
#define S_PRESSURE 4

//data area
__thread_local Real s_step_other_line[3][3][(STEP_SIZE + 2) * 19];//in
__thread_local Real s_step_current_line[STEP_SIZE * 19];
__thread_local int s_step_walls[STEP_SIZE * 19];
__thread_local int s_step_flags[STEP_SIZE];
__thread_local const Real s_w[19] = {
        (1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };
__thread_local const int s_e_x[19] = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0};
__thread_local const int s_e_y[19] = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0};
__thread_local const int s_e_z[19] = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0};
__thread_local const int s_dfInv[19] = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};


void s_stream_and_collide_step()
{
    int inv, l;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew;

    int sis;
    #define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
        if(s_step_flags[sis] == S_FLUID) {
            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                s_step_current_line[sis * 19 + l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
            }
        }
        if(s_step_flags[sis] == S_BOUNCE) {

            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                if(s_step_walls[sis * 19 + l]) {
                    s_step_current_line[sis * 19 + l] = s_step_other_line[1][1][other_step_index(sis, inv)];
                } else 
                {
                    s_step_current_line[sis * 19 + l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                }
            }
        }

        if(s_step_flags[sis] == S_FLUID || s_step_flags[sis] == S_BOUNCE)
        {
            rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
            for (l = 0; l < 19; l++) {
               fi = s_step_current_line[sis * 19 + l];
                rho += fi;
                u_x += s_e_x[l] * fi;
                u_y += s_e_y[l] * fi;
                u_z += s_e_z[l] * fi;
            }

            u_x /= rho;
            u_y /= rho;
            u_z /= rho;

            for (l = 0; l < 19; l++) {
                const Real tmp = (s_e_x[l] * u_x + s_e_y[l] * u_y + s_e_z[l] * u_z);
                feq[l] = s_w[l] * rho * (1.0 -
                                       (3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
                                       (3.0 *     tmp) +
                                       (9.0/2.0 * tmp * tmp));
            }

            Qo=0;
            Real Sij[3][3],S;
            Real e[19][3];
            int x1, y1, k1;
            for(k1 = 0; k1 < 19; k1++) {
                e[k1][0] = s_e_x[k1];
                e[k1][1] = s_e_y[k1];
                e[k1][2] = s_e_z[k1];
            }

            for(x1 = 0; x1 < 3; x1++) {
                for(y1 = 0; y1 < 3; y1++) {
                    Sij[x1][y1] = 0;
                    for(k1 = 0; k1 < 19; k1++) {
                       Sij[x1][y1] += e[k1][x1] * e[k1][y1] * (s_step_current_line[sis * 19 + k1] - feq[k1]);
                    }
                    Qo += Sij[x1][y1] * Sij[x1][y1];
                }
            }

            param.nu = (2.0 / param.omega - 1.0) / 6.0;
            S = (-param.nu + sqrt(param.nu * param.nu + 18 * param.CSmago * param.CSmago * sqrt(Qo))) / 6.0 / param.CSmago / param.CSmago;
            omegaNew = 1.0 / (3.0 * (param.nu + param.CSmago * param.CSmago * S) + 0.5);




            for (l = 0; l < 19; l++) {
               s_step_current_line[sis * 19 + l] =
                       (1.0 - omegaNew) * s_step_current_line[sis * 19 + l] +
                       omegaNew * feq[l];
            }
        }
    }
}

__thread_local int dbg_print_cnt;

__thread_local int dbg_x, dbg_y;

void s_stream_and_collide_xy(int nz,
                           Real* current_line,
                           int* walls,
                           int* flags,
                           Real* other_line)
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

        seq_read_data(s_step_flags, flags + k, sizeof(int) * STEP_SIZE);
        seq_read_data(s_step_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        seq_read_data(s_step_current_line, current_line + k * 19, sizeof(Real) * STEP_SIZE * 19);

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

        int ex, ey;
        for(ex = -1; ex <= 1; ex++)
        {
            for(ey = -1; ey <= 1; ey++)
            {
                seq_read_data(s_step_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                       &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                       sizeof(Real) * (z_ed - z_st + 1) * 19);
                // memcpy(s_step_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                //        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                //        sizeof(Real) * (z_ed - z_st + 1) * 19);
            }
        }

        s_stream_and_collide_step();

        seq_write_data(current_line + k * 19, s_step_current_line, sizeof(Real) * STEP_SIZE * 19);
        // memcpy(current_line + k * 19, s_step_current_line, sizeof(Real) * STEP_SIZE * 19);
    }

}

void std_lbm()
{
    dbg_print_cnt = 0;
    int point_cnt = local_flag[STD_POINT_CNT];
    int* xs_ptr = (int*)local_flag[STD_XS_PTR];
    int* ys_ptr = (int*)local_flag[STD_YS_PTR];
    Real* current_head = (Real*)local_flag[STD_CURRENT_HEAD];
    Real* other_head = (Real*)local_flag[STD_OTHER_HEAD];

    int my_start, my_part, max_part;
    my_start = BLOCK_LOW(slave_id, 64, point_cnt);
    my_part = BLOCK_SIZE(slave_id, 64, point_cnt);
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
        Real* point_current = current_head + s_nodes_idx(x, y, 0, 0);
        Real* point_other = other_head + s_nodes_idx(x, y, 0, 0);
        int* point_walls = param.walls + s_walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = param.flags + s_flags_idx(x, y, 0);

        s_stream_and_collide_xy(param.Z, point_current, point_walls, 
                              point_flags, point_other);
    }

    standard_write_flag(slave_id);
}

void insane_lbm()
{

}


// 其他暂时未用到的函数，分开放，避免文件过长，也方便多人编辑
// 全部完成后复制到下方
#include "SlaveKernalFunc.c"
#include "SlaveWorks.c"
