#ifndef __LBM_MASTER_H__
#define __LBM_MASTER_H__

#include "../header/Argument.h"
#include "../header/Array.h"
#include "../header/LbmSlave.h"
#include "../header/Config.h"
#include "athread.h"

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "mpi.h"
#include <assert.h>
#include <memory.h>

#define CCPS 1450000000.0

extern void SLAVE_FUN(cpe_athread_daemon)();

int my_rank, comm_sz;

MPI_Comm *mycomm; 
int *dims, *coords;

int X, Y, Z, Xst, Xed, Yst, Yed, x_sec, y_sec;
int STEPS;

Real *****ori_nodes;
int ***ori_flags, ****ori_walls;

Real ***right_recv;
Real ***left_recv;
Real ***down_recv;
Real ***up_recv;
Real ***right_send;
Real ***left_send;
Real ***down_send;
Real ***up_send;

Real **lu_recv;
Real **ld_recv;
Real **ru_recv;
Real **rd_recv;
Real **lu_send;

Real **ld_send;
Real **ru_send;
Real **rd_send;

void wait_slave_flag();
void terminate_athread_daemon();
void athread_handshake();

volatile long host_flag[FLAG_SIZE];
volatile long slave_flag[FLAG_SIZE];
volatile long flag_to_wait;

// put data into global field
void data_init( int _myrank, int _comm_sz, MPI_Comm *_mycomm, 
		        int *_dims, int *_coords,
                int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec, 
                Real *****_nodes, int ***_flags, int ****_walls,
                Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send);

//void main_iter();

void lbm_data_init(int _X, int _Y, int _Z, int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec,
                   Real *****_ori_nodes, int ***_ori_flags, int ****_ori_walls, int _STEPS,
                   Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                   Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                   Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                   Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send);
void main_iter(int *_s);


static inline double ccts(unsigned long cc)
{
    return cc * 1.0 / CCPS;
}



#endif
