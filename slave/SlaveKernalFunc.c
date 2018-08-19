#include "../header/LbmSlave.h"

// 暂时没用到，以后可能会用的东西。。。先放这里
#ifndef __thread_local
#define __thread_local const //本地IDE美化(__thread_local 是神威从核编译的关键字, 相当于cuda的片上内存)
#endif

int BLOCK_SIZE(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) + ((n % total_blocks > block_id) ? 1 : 0);
}

int BLOCK_LOW(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) * block_id + ((n % total_blocks > block_id) ? block_id : n % total_blocks);
}

void* nw_ldm_malloc(int sz)//will not see FBI WARNING during compile
{
    long ret = ldm_malloc(sz);
    return (void*)ret;
}

//conghe profile module (￣△￣；)
__thread_local long pf_stcc_, pf_edcc_;
inline void stpf()
{
    RPCC(pf_stcc_);
}

inline void edpf(int n)
{
    RPCC(pf_edcc_);
    out_flag[n] += pf_edcc_ - pf_stcc_;
    pf_stcc_ = pf_edcc_;
}

//align
void* _nice_ptr_head(char* ptr, int alignment, int right_shft)
{
    long ptrl = (long)ptr;
    int rem = ptrl % alignment;
    if(rem > 0)
    {
        return ptr + (alignment - rem + right_shft);
    }
    else
    {
        return ptr + right_shft;
    }
}

inline static void wait_all_async_write_data()
{
    while(put_reply_god != put_reply_target);
}

inline static void wait_all_async_get_data()
{
    while(get_reply_god != get_reply_target);
}
