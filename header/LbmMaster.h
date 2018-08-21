#ifndef __LBM_MASTER_H__
#define __LBM_MASTER_H__
#include "../header/LbmSlave.h"
#include "../header/Config.h"
#include "../header/Argument.h"

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "mpi.h"
#include "athread.h"
#include <assert.h>
#include <memory.h>

//static inline unsigned long rpcc()
//{
//    unsigned long time;
//    asm("rtc %0": "=r" (time) : );
//    return time;
//}

static inline double ccts(unsigned long cc)
{
    return cc * 1.0 / CCPS;
}


extern void SLAVE_FUN(cpe_athread_daemon)();

void lbm_data_init(int _X, int _Y, int _Z, int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec,
                   realt *****_ori_nodes, int ***_ori_flags, int ****_ori_walls, int _STEPS,
                   Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                   Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                   Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                   Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send);
void main_iter(int *_s);

#endif
