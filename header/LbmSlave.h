#ifndef __LBM_SLAVE_H__
#define __LBM_SLAVE_H__

typedef struct lbm_param
{
    int master_id;

    long *host_flag;
    long *slave_flag;

    int iter;
    int x_sec, y_sec, Z;
} lbm_param;

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

#endif
