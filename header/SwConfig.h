#ifndef __SW_CONFIG_H__
#define __SW_CONFIG_H__

#define PARAM_SIZE                          32
    #define PARAM_MPI_RANK                  1
    #define PARAM_DEVICE_ACTION             2
        #define DEVICE_ACTION_EXIT                  999
        #define DEVICE_ACTION_RUN                   666
    // ---- CUSTOM PARAM
    #define STD_POINT_CNT        3
    #define STD_XS_PTR           4
    #define STD_YS_PTR           5
    #define STD_CURRENT_HEAD     6
    #define STD_OTHER_HEAD       7
    #define INSANE_POINT_CNT     8
    #define INSANE_XS_PTR        9
    #define INSANE_YS_PTR        10
    // ----

#define CORE_SYNC_64 0x0000FFFF
#define CORE_SYNC_8  0x0000000F

#define DEVICE_SAFE_PAD 48

struct InitParam
{
    int my_id;
    long* host_flag, *slave_flag;

    // ---- CUSTOM VAR
    int iter;
    int x_sec, y_sec, Z;

    int* walls;
    int* flags;

    float nu, omega, CSmago;
    // ----
};

// ---- CUSTOM CONFIG
#define STEP_SIZE 50
#define INSANE_SIZE 20

#define BLOCK_SIZE(ID, ID_SZ, N) (((N) / (ID_SZ)) + (((N) % (ID_SZ) > (ID)) ? 1 : 0))
#define BLOCK_HEAD(ID, ID_SZ, N) (((N) / (ID_SZ)) * (ID) + (((N) % (ID_SZ) > (ID)) ? (ID) : (N) % (ID_SZ)))
// ----

#endif
