#ifndef __LBM_SLAVE_H__
#define __LBM_SLAVE_H__

#ifndef Real
#define Real float
#endif

typedef struct lbm_param
{
    int master_id;
    // master传入数据起始地址，slave接收地址
    long* host_flag, *slave_flag;
    // 迭代次数
    int iter;
    int x_sec, y_sec, Z;
    int* walls;
    int* flags;

    float nu, omega, CSmago;
} lbm_param;

void cpe_athread_daemon(void *_param);

#define EXIT_FLAG        12450
#define HANDSHAKE_FLAG   19260817
#define STD_LBM_FLAG     2333
#define INSANE_LBM_FLAG  114514

//only prepared for
#define FLAG_SIZE        32
#define MPI_RANK         1
#define KERNEL_ACTION    2 // work and exit

//STD_LBM
#define STD_POINT_CNT        3
#define STD_XS_PTR           4
#define STD_YS_PTR           5
#define STD_CURRENT_HEAD     6
#define STD_OTHER_HEAD       7

//step 40, 3 for rem, 39 for step rem, and extra 6 for safety
#define SLAVE_SAFE_PAD 48
#define STEP_SIZE 50
#define CPE_TOTAL_SYNC 0x0000FFFF
#define CPE_TOTAL_SYNC 0x0000FFFF

#define BLOCK_SIZE(block_id, total_blocks, n) (((n) / (total_blocks)) + (((n) % (total_blocks) > (block_id)) ? 1 : 0))
#define BLOCK_LOW(block_id, total_blocks, n) (((n) / (total_blocks)) * (block_id) + (((n) % (total_blocks) > (block_id)) ? (block_id) : (n) % (total_blocks)))


#endif
