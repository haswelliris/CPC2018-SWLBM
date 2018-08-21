#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int my_rank, comm_sz;

//hch intel timer
#define MAX_HCH_TIMER_LEN 256
#define HCH_TIMER_CNT 6

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
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    int i;
    for(i = 0; i < MAX_HCH_TIMER_LEN; i++)
        hch_cc[i] = 0;
}

//hch intel timer

void hch_timer_manual(long* arr, int cnt, const char* fname, const char* header)
{
    if(my_rank == 0)
    {
        FILE * fp;
        fp = fopen (fname, "w");
        fprintf(fp, header);
        fprintf(fp, "\n");
        fflush(fp);
        fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int i, j, k;
    for(i = 0; i < comm_sz; i++)
    {
        if(my_rank == i)
        {
            FILE * fp;

            fp = fopen (fname, "a");

            fprintf(fp, "%d,", my_rank);

            for(j = 0; j < cnt; j++)
            {
                fprintf(fp, "%.4f,", arr[j] * 1.0 / CCPS);
            }

            fprintf(fp, "\n");

            fflush(fp);
            fclose(fp);

            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void hch_timer_finalize_()
{
    if(my_rank == 0)
    {
        printf("hch_timer_finalize_\n");
        FILE * fp;
        fp = fopen ("hch_timer_profile.csv", "w");
        fprintf(fp, "my_rank, MPF1_SEND_INIT, MPF1_COMM, MPF1_WAIT_MPI, MPF1_UPDATE, MPF1_WAIT_INNER, MPF1_HALO");
        fprintf(fp, "\n");
        fflush(fp);
        fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int i, j, k;
    for(i = 0; i < comm_sz; i++)
    {
        if(my_rank == i)
        {
            FILE * fp;

            fp = fopen ("hch_timer_profile.csv", "a");

            // printf("%d,", my_rank);
            fprintf(fp, "%d,", my_rank);

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
    //printf("my_rank = %d, num = %d\n", my_rank, *num);
    hch_tmp_cc[num] = (HCH_CC_TYPE)GetCycleCount();
}

void hch_timer_stop(int num)
{
    //printf("my_rank = %d, num = %d\n", my_rank, *num);

    hch_cc[num] += (HCH_CC_TYPE)GetCycleCount() - hch_tmp_cc[num];
}

void hch_timer_start_(int* num)
{
    //printf("my_rank = %d, num = %d\n", my_rank, *num);
    hch_tmp_cc[*num] = (HCH_CC_TYPE)GetCycleCount();
}

void hch_timer_stop_(int* num)
{
    //printf("my_rank = %d, num = %d\n", my_rank, *num);

    hch_cc[*num] += (HCH_CC_TYPE)GetCycleCount() - hch_tmp_cc[*num];
}

