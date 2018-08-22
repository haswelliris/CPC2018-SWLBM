#ifndef __LBM_MASTER_H__
#define __LBM_MASTER_H__
#include "../header/Argument.h"
#include "../header/SwHost.h"

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "mpi.h"
#include <assert.h>
#include <memory.h>

void lbm_data_init(int _X, int _Y, int _Z, int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec,
                   Real *****_ori_nodes, int ***_ori_flags, int ****_ori_walls, int _STEPS,
                   Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                   Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                   Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                   Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send);
void main_iter(int *_s);

#endif
