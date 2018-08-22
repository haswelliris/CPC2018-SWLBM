#include "../header/SwHost.h"

#include <stdlib.h>
#include "mpi.h"
#include "athread.h"

volatile long host_flag[PARAM_SIZE];
volatile long slave_flag[PARAM_SIZE];
volatile long local_cc[PARAM_SIZE];
volatile long flag_to_wait;
unsigned long mpe_cc_cur_[100];
unsigned long mpe_cc_total_[100];
int myrank, comm_sz;
struct InitParam host_param;

void wait_slave_flag()
{
    while(slave_flag[0] < flag_to_wait);
    flag_to_wait++;
}

void terminate_athread_daemon()
{
    host_flag[PARAM_MPI_RANK] = myrank;
    host_flag[PARAM_DEVICE_ACTION] = DEVICE_ACTION_EXIT;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();

    athread_join();
}
