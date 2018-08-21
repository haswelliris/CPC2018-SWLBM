#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int _my_rank, _comm_sz;

//hch intel timer
#define MAX_HCH_TIMER_LEN 256
#define HCH_TIMER_CNT 8

#if defined (__i386__)
#define HCH_CC_TYPE unsigned long long
#define CCPS 2900000000
static __inline__ unsigned long long GetCycleCount(void)
{
        unsigned long long int x;
        __asm__ volatile("rdtsc":"=A"(x));
        return x;
}
#elif defined (__x86_64__)
#define HCH_CC_TYPE unsigned long long
#define CCPS 2900000000
static __inline__ unsigned long long GetCycleCount(void)
{
        unsigned hi,lo;
        __asm__ volatile("rdtsc":"=a"(lo),"=d"(hi));
        return ((unsigned long long)lo)|(((unsigned long long)hi)<<32);
}
#elif (defined SW2) || (defined SW5)
#define HCH_CC_TYPE unsigned long
#define CCPS 1450000000
static __inline__ unsigned long GetCycleCount(void)
{
    unsigned long time;
    asm("rtc %0": "=r" (time) : );
    return time;
}
#endif

HCH_CC_TYPE hch_cc[MAX_HCH_TIMER_LEN];
HCH_CC_TYPE hch_tmp_cc[MAX_HCH_TIMER_LEN];

void hch_timer_init_()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &_my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_comm_sz);

    int i;
    for(i = 0; i < MAX_HCH_TIMER_LEN; i++)
        hch_cc[i] = 0;
}

//hch intel timer

void hch_timer_finalize_()
{
    if(_my_rank == 0)
    {
        // printf("hch_timer_finalize_\n");
        FILE * fp;
        fp = fopen ("hch_timer_profile.csv", "w");
        fprintf(fp, "my_rank, athread init, send init, send, wait, update, stream, collide, athread join");
        fprintf(fp, "\n");
        fflush(fp);
        fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int i, j, k;
    for(i = 0; i < _comm_sz; i++)
    {
        if(_my_rank == i)
        {
            FILE * fp;

            fp = fopen ("hch_timer_profile.csv", "a");

            // printf("%d,", _my_rank);
            fprintf(fp, "%d,", _my_rank);

            #ifndef HCH_TIMER_CNT
            for(j = 0; j < MAX_HCH_TIMER_LEN; j++)
            {
                if(hch_cc[j] == 0)
                    break;
                // printf("%.4f,", hch_cc[j] * 1.0 / CCPS);
                fprintf(fp, "%.4f,", hch_cc[j] * 1.0 / CCPS);
            }
            #else
            for(j = 0; j < HCH_TIMER_CNT; j++)
            {
                // printf("%.4f,", hch_cc[j] * 1.0 / CCPS);
                fprintf(fp, "%.4f,", hch_cc[j] * 1.0 / CCPS);
            }
            #endif

            fprintf(fp, "\n");
            // printf("\n");

            fflush(fp);
            fclose(fp);

            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void hch_timer_start(int num)
{
    //printf("_my_rank = %d, num = %d\n", _my_rank, *num);
    hch_tmp_cc[num] = (HCH_CC_TYPE)GetCycleCount();
}

void hch_timer_stop(int num)
{
    //printf("_my_rank = %d, num = %d\n", _my_rank, *num);

    hch_cc[num] += (HCH_CC_TYPE)GetCycleCount() - hch_tmp_cc[num];
}

void hch_timer_start_(int* num)
{
    //printf("_my_rank = %d, num = %d\n", _my_rank, *num);
    hch_tmp_cc[*num] = (HCH_CC_TYPE)GetCycleCount();
}

void hch_timer_stop_(int* num)
{
    //printf("_my_rank = %d, num = %d\n", _my_rank, *num);

    hch_cc[*num] += (HCH_CC_TYPE)GetCycleCount() - hch_tmp_cc[*num];
}