#ifndef __SW_HOST_H__
#define __SW_HOST_H__
#include "../header/SwConfig.h"

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "athread.h"

extern volatile long host_flag[PARAM_SIZE];
extern volatile long slave_flag[PARAM_SIZE];
extern volatile long local_cc[PARAM_SIZE];
extern volatile long flag_to_wait;
extern unsigned long mpe_cc_cur_[100];
extern unsigned long mpe_cc_total_[100];
extern int myrank, comm_sz;
extern struct InitParam host_param;

void wait_slave_flag();

void terminate_athread_daemon();

#endif
