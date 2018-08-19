#ifndef __LBM_SLAVE_H__
#define __LBM_SLAVE_H__

#include "../header/Array.h"

typedef struct lbm_param
{
    int myrank;

    long *host_flag;
    long *slave_flag;

    int STEPS;
    int x_sec, y_sec, Z;
} lbm_param;

#endif
