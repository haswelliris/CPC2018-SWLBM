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
__thread_local volatile intv8 bcast_tmp256;                 // 8B

void cpe0_bcast_256(int slave_id, volatile intv8 *bcast_data);
void standard_wait_flag(int slave_id);
void standard_write_flag(int slave_id);
void seq_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int slave_id, int write_flag);
void async_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int slave_id, int write_flag);
void std_lbm();
void insane_lbm();

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
        printf("@@@ MESSAGE FROM CPE 0 @@@\n");
        printf("iter = %d, x_sec = %d, y_sec = %d, z_sec = %d\n", param.iter, param.x_sec, param.y_sec, param.Z);
        printf("ptrs: %x, %x\n", param.host_flag, param.slave_flag);
        printf("@@@ MESSAGE FROM CPE 0 @@@\n");
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
        printf("@@@ MESSAGE FROM CPE 0 @@@\n");
        printf("CPE ENTER LOOP\n");
        printf("@@@ MESSAGE FROM CPE 0 @@@\n");
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
        else if (local_flag[KERNEL_ACTION] == STD_LBM_FLAG)
        {
            std_lbm();
        }
        else if (local_flag[KERNEL_ACTION] == INSANE_LBM_FLAG)
        {
            insane_lbm();
        }
        else
        {
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
            asm volatile ("#nop":::"memory");
            if(local_flag[0] >= slave_flag_to_wait)
                break;
        }

        slave_flag_to_wait++;
    }

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

void seq_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int slave_id, int write_flag)
{
//    int wcid;
//    for (wcid = 63; wcid >= 0; wcid--)
//    {
//        if (slave_id == wcid && write_flag != 0)
//        {
            put_reply = 0;
            athread_put(
                PE_MODE, cpe_data_ptr, mpe_data_ptr,
                size, (void*)(&put_reply), 0, 0
            );
            while (put_reply != 1);
//        }
//        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
//    }
}

//需要异步读取数据

void async_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int slave_id, int write_flag)
{
    int wcid;
    //for (wcid = 63; wcid >= 0; wcid--)
    //{
    //    if (slave_id == wcid && write_flag != 0)
    //    {
            //put_reply = 0;
            athread_put(
                PE_MODE, cpe_data_ptr, mpe_data_ptr,
                size, (void*)(&put_reply_god), 0, 0
            );
            put_reply_target++;
            //while (put_reply != 1);
    //    }
    //    //athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    //}
}

//std area

void std_lbm()
{
    //get the x & y cnt

}

void insane_lbm()
{

}
