#ifndef __LBM_SLAVE_H__
#define __LBM_SLAVE_H__

#include <slave.h>

#include "dma.h"
#include "simd.h"
#include <stdio.h>
#include <memory.h>


const int slave_id, core_x, core_y;  // 从核id，所在行，所在列
// 我觉得是7 * 4B + 4B
const struct lbm_init_param param;   // 28 * (4B / 8B) + 4B
//  cache line 128 Bytes. ( 32B * 4 )
const volatile long local_flag[FLAG_SIZE];        // 16 * 8B or 32 * 4B
const volatile long out_flag[FLAG_SIZE];          // 16 * 8B or 32 * 4B
const long slave_flag_to_wait;
const volatile unsigned long get_reply, put_reply; // 2 * 8B
const volatile unsigned long get_reply_god, put_reply_god; // 2 * 8B
const volatile unsigned long get_reply_target, put_reply_target; // 2 * 8B
const volatile intv8 bcast_tmp256;                 // 8B

struct lbm_init_param
{
    // 主核进程的mpi rank
    int master_id;
    // 用于控制程序步骤的flags
    long* host_flag, *slave_flag;
    // 迭代次数
    int iter;
    // 分割之后的数据大小
    int x_sec, y_sec, Z;
};

void cpe_athread_daemon(void *_param);

#define EXIT_FLAG        12450
#define STD_LBM_FLAG     2333
#define INSANE_LBM_FLAG  114514

//only prepared for
#define FLAG_SIZE        32
#define MPI_RANK         1
#define KERNEL_ACTION    2 // work and exit

//std lbm
//#define GROUP_SIZE       3
//#define REMAIN_POINT     4
//#define IN_PTR           5
//#define OUT_PTR          6
//#define IN_STRIDE        7
//#define OUT_STRIDE       8
//#define REQUIRE_IN       9
//#define REQUIRE_OUT      10

//step 40, 3 for rem, 39 for step rem, and extra 6 for safety
#define SLAVE_SAFE_PAD 48
#define STEP_SIZE 50
#define CPE_TOTAL_SYNC 0x0000FFFF

void cpe0_bcast_256(int slave_id, volatile intv8* bcast_data);
void standard_wait_flag(int slave_id);
void standard_write_flag(int slave_id);
void seq_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int slave_id, int write_flag);
void async_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int slave_id, int write_flag);
void std_lbm();
void insane_lbm();

#endif
