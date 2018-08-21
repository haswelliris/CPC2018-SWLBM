#include "../header/LbmSlave.h"
#include "../header/Config.h"

#include "slave.h"
#include "dma.h"
#include "simd.h"
#include <stdio.h>
#include <memory.h>
#include <math.h>

#ifndef __thread_local
#define __thread_local const //本地IDE美化
#endif

#ifndef realt
#define realt float
#endif


__thread_local int my_core, core_x, core_y;
__thread_local struct lbm_init_param param;                 // 28 * (4B / 8B) + 4B
__thread_local volatile long local_flag[FLAG_SIZE];          // 16 * 8B
__thread_local volatile long out_flag[FLAG_SIZE];          // 16 * 8B
__thread_local long slave_flag_to_wait;
__thread_local volatile unsigned long get_reply, put_reply; // 2 * 8B
__thread_local volatile unsigned long get_reply_god, put_reply_god; // 2 * 8B
__thread_local volatile unsigned long get_reply_target, put_reply_target; // 2 * 8B
__thread_local volatile intv8 bcast_tmp256;                 // 8B

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

void cpe0_bcast_256(int my_core, volatile intv8* bcast_data)
{
    bcast_tmp256 = *bcast_data;

    if (my_core / 8 == 0)
    {
        __builtin_sw64_putC(bcast_tmp256, 0x0000000F);
    } else {
        bcast_tmp256 = __builtin_sw64_getxC(bcast_tmp256);
    }

    if (my_core % 8 == 0)
    {
        asm volatile ("nop":::"memory");
        __builtin_sw64_putR(bcast_tmp256, 0x0000000F);
        *bcast_data = bcast_tmp256;
    } else {
        asm volatile ("nop":::"memory");
        *bcast_data = __builtin_sw64_getxR(bcast_tmp256);
    }
}

void standard_wait_flag(int my_core)
{
    if(my_core == 0)
    {
        //printf("SPE waiting for %ld\n", slave_flag_to_wait);

        while(1)
        {
            get_reply = 0;
            athread_get(
                PE_MODE, param.host_flag, (void*)&local_flag[0],
                sizeof(long) * FLAG_SIZE, (void*)(&get_reply), 0, 0, 0
            );
            while(get_reply != 1);
            asm volatile ("#nop":::"memory");
            if(local_flag[0] >= slave_flag_to_wait)
                break;
        }

        slave_flag_to_wait++;
    }

    int i;
    for(i = 0; i < FLAG_SIZE / 4; i++)
    {
        cpe0_bcast_256(my_core, (intv8*)&local_flag[i * 4]);
    }
}

void standard_write_flag(int my_core)
{
    out_flag[0] = out_flag[0] + 1;

    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);//this is necessary, thanks to Luo Yangze
    // athread_dma_barrier();
    if(my_core == 0)
    {
        put_reply = 0;
        athread_put(
            PE_MODE, (void*)&out_flag[0], param.slave_flag,
            sizeof(long) * FLAG_SIZE, (void*)(&put_reply), 0, 0
        );
        while(put_reply != 1);

    }
    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
}

void seq_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size)
{
    put_reply = 0;
    athread_put(
        PE_MODE, cpe_data_ptr, mpe_data_ptr,
        size, (void*)(&put_reply), 0, 0
    );
    while (put_reply != 1);
}

void seq_read_data(void* cpe_data_ptr, void* mpe_data_ptr, int size)
{
    get_reply = 0;
    athread_get(
                PE_MODE, mpe_data_ptr, cpe_data_ptr,
                size, (void*)(&get_reply), 0, 0, 0
            );
    while (get_reply != 1);
}

//需要异步读取数据
void async_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size)
{
    athread_put(
        PE_MODE, cpe_data_ptr, mpe_data_ptr,
        size, (void*)(&put_reply_god), 0, 0
    );
    put_reply_target++;
}

void async_read_data(void* cpe_data_ptr, void* mpe_data_ptr, int size)
{
    athread_get(
        PE_MODE, mpe_data_ptr, cpe_data_ptr,
        size, (void*)(&get_reply_god), 0, 0, 0
    );
    get_reply_target++;
}

//std area

//const data
#define S_FLUID    0
#define S_NOSLIP   1
#define S_VELOCITY 2
#define S_BOUNCE   3
#define S_PRESSURE 4

//data area
__thread_local realt s_step_other_line[3][3][(STEP_SIZE + 2) * 19] __attribute__((aligned(32)));//in
__thread_local realt s_step_current_line[STEP_SIZE * 19] __attribute__((aligned(32)));
__thread_local int s_step_walls[STEP_SIZE * 19] __attribute__((aligned(32)));
__thread_local int s_step_flags[STEP_SIZE] __attribute__((aligned(32)));
__thread_local realt s_w[19] __attribute__((aligned(32))) = {
        (1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };
__thread_local int s_e_x[19] __attribute__((aligned(32))) = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0};
__thread_local int s_e_y[19] __attribute__((aligned(32))) = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0};
__thread_local int s_e_z[19] __attribute__((aligned(32)))  = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0};
__thread_local realt s_e_T[19][4] __attribute__((aligned(32))) = {
	{0.,1.,0.,0}, {0.,-1.,0.,0.}, {1.,0.,0.,0.}, {-1.,0.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,-1.,0.},
	{1.,1.,0.,0.}, {-1.,1.,0.,0.}, {1.,-1.,0.,0.}, {-1.,-1.,0.,0.}, {0.,1.,1.,0.}, {0.,1.,-1.,0.},
	{0.,-1.,1.,0.}, {0.,-1.,-1.,0.}, {1.,0.,1.,0.}, {1.,0.,-1.,0.}, {-1.,0.,1.,0.}, {-1.,0.,-1.,0.}, {0.,0.,0.,0.}
};
__thread_local int s_dfInv[19] __attribute__((aligned(32))) = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};
// __thread_local realt  a_current_line[20] __attribute__((aligned(32)));

void s_stream_step()
{
    int inv, l;
	int sis;
    #define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
	//stream
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
		
		// for(l = 0; l < 19; l++) {
			// a_current_line[l] = s_step_current_line[sis * 19 + l];
		// }
        if(s_step_flags[sis] == S_FLUID) {
            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
				// a_current_line[l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                s_step_current_line[sis * 19 + l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
            }
        }
        if(s_step_flags[sis] == S_BOUNCE) {

            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                if(s_step_walls[sis * 19 + l]) {
					// a_current_line[l] = s_step_other_line[1][1][other_step_index(sis, inv)];
                    s_step_current_line[sis * 19 + l] = s_step_other_line[1][1][other_step_index(sis, inv)];
                } else 
                {
					// a_current_line[l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                    s_step_current_line[sis * 19 + l] = s_step_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                }
            }
        }
	}
}

__thread_local realt  Sij[12] __attribute__((aligned(32)));
__thread_local realt  feq[20] __attribute__((aligned(32)));
__thread_local realt  res[20] __attribute__((aligned(32)));
__thread_local    realt _prho, _pu_x, _pu_y, _pu_z;
__thread_local    realt _pfi;
__thread_local    realt _pQo, _pomegaNew;
__thread_local    int _psis, _pl;
__thread_local    realt _ptmp;
__thread_local    realt _pS;
__thread_local	  floatv4 _pvs1, _pvs2, _pvs3,_pve, _pvk1, _pvk2, _pvk3, _pvk4, _pvk5;
__thread_local    int _px1, _py1, _pk1;
__thread_local    floatv4 _pU, _pE;

s_collide_step_a() {

}

s_collide_step_b() {
    _pU = simd_set_floatv4(0., 0., 0., 0.), _pE;
    #pragma unroll[19]
    for (_pl = 0; _pl < 19; _pl++) {
		// fi = a_current_line[l];
		_pfi = s_step_current_line[_psis * 19 + _pl];
		_pE = simd_set_floatv4(1., s_e_T[_pl][0], s_e_T[_pl][1], s_e_T[_pl][2]);
		_pU = _pfi * _pE + _pU;
    }
	simd_store(_pU, &(Sij[0]));
	_prho = Sij[0], _pu_x = Sij[1], _pu_y = Sij[2], _pu_z = Sij[3];
}

s_collide_step_c() {
    _pu_x /= _prho;
    _pu_y /= _prho;
    _pu_z /= _prho;
}

s_collide_step_d() {
    //用simd会丢精度, 但是会快(2/100)s/次
    #pragma unroll[19]
    for (_pl = 0; _pl < 19; _pl++) {
        _ptmp = (s_e_x[_pl] * _pu_x + s_e_y[_pl] * _pu_y + s_e_z[_pl] * _pu_z);
        feq[_pl] = s_w[_pl] * _prho * (1.0 -
                               (3.0/2.0 * (_pu_x * _pu_x + _pu_y * _pu_y + _pu_z * _pu_z)) +
                               (3.0 *     _ptmp) +
                               (9.0/2.0 * _ptmp * _ptmp));
    }
}

s_collide_step_e() {
    _pvs1 = simd_set_floatv4(0., 0., 0., 0.);
	_pvs2 = simd_set_floatv4(0., 0., 0., 0.);
	_pvs3 = simd_set_floatv4(0., 0., 0., 0.);
    #pragma unroll[19]
	for(_pk1 = 0; _pk1 < 19; _pk1++) {
		_ptmp = (s_step_current_line[_psis * 19 + _pk1] - feq[_pk1]);
		// float tmp = (a_current_line[k1] - feq[k1]);
		simd_load(_pve, &(s_e_T[_pk1][0]));
		_pvk1 = _ptmp * s_e_T[_pk1][0];
		_pvk2 = _ptmp * s_e_T[_pk1][1];
		_pvk3 = _ptmp * s_e_T[_pk1][2];
				
		_pvs1 = _pvk1 * _pve + _pvs1;
		_pvs2 = _pvk2 * _pve + _pvs2;
		_pvs3 = _pvk3 * _pve + _pvs3;
    }
	simd_store(_pvs1, &(Sij[0]));
	simd_store(_pvs2, &(Sij[4]));
	simd_store(_pvs3, &(Sij[8]));
}

s_collide_step_f() {
    //用simd乘加会丢精度，也不会变快
    _pQo=0;
    #pragma unroll[12]
	for(_px1 = 0; _px1 < 12; ++_px1)
	{
		_pQo += Sij[_px1] * Sij[_px1];
	}
}

s_collide_step_g() {
    param.nu = (2.0 / param.omega - 1.0) / 6.0;
    _pS = (-param.nu + sqrt(param.nu * param.nu + 18 * param.CSmago * param.CSmago * sqrt(_pQo))) / 6.0 / param.CSmago / param.CSmago;
    _pomegaNew = 1.0 / (3.0 * (param.nu + param.CSmago * param.CSmago * _pS) + 0.5);
}

s_collide_step_h() {
    //没有必要上simd
    #pragma unroll[19]
    for (_pl = 0; _pl < 19; _pl++) {
        res[_pl] =
              (1.0 - _pomegaNew) * s_step_current_line[_psis * 19 + _pl] +
              _pomegaNew * feq[_pl];
    }
}

s_collide_step_i() {
    //没有必要上simd
    #pragma unroll[19]
    for (_pl = 0; _pl < 19; _pl++) {
        s_step_current_line[_psis * 19 + _pl] = res[_pl];
    }
}

//半SIMD化
void s_collide_step(int st, int ed)
{
	//coli
	for(_psis = st; _psis < ed; _psis++)
    {
        if(s_step_flags[_psis] == S_FLUID || s_step_flags[_psis] == S_BOUNCE)
        {
            s_collide_step_a();
            s_collide_step_b();
            s_collide_step_c();
            s_collide_step_d();
            s_collide_step_e();
            s_collide_step_f();
            s_collide_step_g();
            s_collide_step_h();
            s_collide_step_i();
        }
    }
}

//TODO !!!!修了一些bug，但是精度丢失还是很严重，变快5s
__thread_local float Qo4[4] __attribute__((aligned(32)));
__thread_local float S4[4] __attribute__((aligned(32)));
__thread_local float omegaNew4[4] __attribute__((aligned(32)));
__thread_local float line4[4] __attribute__((aligned(32)));
__thread_local float u_x4[4] __attribute__((aligned(32)));
__thread_local float u_y4[4] __attribute__((aligned(32)));
__thread_local float u_z4[4] __attribute__((aligned(32)));
__thread_local float rho4[4] __attribute__((aligned(32)));
__thread_local float tmp4[4] __attribute__((aligned(32)));

#define STEP_SIZE_MAIN ((STEP_SIZE / 4) * 4)

__thread_local    floatv4  _vpfeqv[19];
__thread_local    floatv4  _vpSijv[9];
__thread_local    floatv4  _vpcurrentv[19];
__thread_local
__thread_local    floatv4 _vprho, _vpu_x, _vpu_y, _vpu_z;
__thread_local	  floatv4 _vpfi;
__thread_local    floatv4 _vpQo, _vpomegaNew;
__thread_local	  int _vpsis, _vpl, _vpii;
__thread_local    floatv4 _vptmp;
__thread_local    floatv4 _vpS;
__thread_local    int _vpx1, _vpy1, _vpk1;

void s_simd_collide_step_a()
{
    for (_vpl = 0; _vpl < 19; _vpl++) {
        for(_vpii = 0; _vpii < 4; ++_vpii)
        {
            if(s_step_flags[_vpsis + _vpii] == S_FLUID || s_step_flags[_vpsis + _vpii] == S_BOUNCE)
			{
				line4[_vpii] = s_step_current_line[(_vpsis + _vpii) * 19 + _vpl];
			}
            else
            {
                line4[_vpii] = 1.0;
            }
        }
		simd_load(_vpcurrentv[_vpl], &(line4[0]));
    }
}

void s_simd_collide_step_b()
{
    _vprho = simd_set_floatv4(0., 0., 0., 0.);
    _vpu_x = simd_set_floatv4(0., 0., 0., 0.);
    _vpu_y = simd_set_floatv4(0., 0., 0., 0.);
    _vpu_z = simd_set_floatv4(0., 0., 0., 0.);
    for (_vpl = 0; _vpl < 19; _vpl++) {
		_vpfi = _vpcurrentv[_vpl];
        _vprho = _vprho + _vpfi;
        _vpu_x = s_e_x[_vpl] * _vpfi + _vpu_x;
        _vpu_y = s_e_y[_vpl] * _vpfi + _vpu_y;
        _vpu_z = s_e_z[_vpl] * _vpfi + _vpu_z;
    }
}

void s_simd_collide_step_c()
{
    _vpu_x = _vpu_x / _vprho;
    _vpu_y = _vpu_y / _vprho;
    _vpu_z = _vpu_z / _vprho;
}

void s_simd_collide_step_d()
{
    //扩展成普通计算精度更爆炸了，为什么？
    for (_vpl = 0; _vpl < 19; _vpl++) {
        _vptmp = (s_e_x[_vpl] * _vpu_x + s_e_y[_vpl] * _vpu_y + s_e_z[_vpl] * _vpu_z);
		_vpfeqv[_vpl] = s_w[_vpl] * _prho * (1.0 -
			(3.0/2.0 * (_vpu_x * _vpu_x + _vpu_y * _vpu_y + _vpu_z * _vpu_z)) +
			(3.0 *     _vptmp) +
			(9.0/2.0 * _vptmp * _vptmp));
    }
}

void s_simd_collide_step_e()
{
    for(_vpx1 = 0; _vpx1 < 9; _vpx1++) {
		_vpSijv[_vpx1] = simd_set_floatv4(0., 0., 0., 0.);
	}
	for(_vpk1 = 0; _vpk1 < 19;_vpk1++) {
		for(_vpx1 = 0; _vpx1 < 3; _vpx1++) for(_vpy1 = 0; _vpy1 < 3; _vpy1++) {
            _vpSijv[_vpx1 * 3 + _vpy1] = \
                s_e_T[_vpk1][_vpx1] * s_e_T[_vpk1][_vpy1] * (_vpcurrentv[_vpk1] - _vpfeqv[_vpk1]) + \
                _vpSijv[_vpx1 * 3 + _vpy1];
		}
    }
}

void s_simd_collide_step_f()
{
    _vpQo = simd_set_floatv4(0., 0., 0., 0.);
	for(_vpx1 = 0; _vpx1 < 9; _vpx1++) {
		_vpQo = _vpSijv[_vpx1] * _vpSijv[_vpx1] + _vpQo;
	}
}

void s_simd_collide_step_g()
{
    simd_store(_vpQo, &(Qo4[0]));
    param.nu = (2.0 / param.omega - 1.0) / 6.0;
    for(_vpii = 0; _vpii < 4; _vpii++)
    {
        S4[_vpii] = (-param.nu + sqrt(param.nu * param.nu + 18 * param.CSmago * param.CSmago * sqrt(Qo4[_vpii]))) / 6.0 / param.CSmago / param.CSmago;
		omegaNew4[_vpii] = 1.0 / (3.0 * (param.nu + param.CSmago * param.CSmago * S4[_vpii]) + 0.5);
    }
    simd_load(_vpomegaNew, &(omegaNew4[0]));
}

void s_simd_collide_step_h()
{
    for (_vpl = 0; _vpl < 19; _vpl++) {
		_vpcurrentv[_vpl] = (1.0 - _vpomegaNew) * _vpcurrentv[_vpl] + _vpomegaNew * _vpfeqv[_vpl];
    }
}

void s_simd_collide_step_i()
{
    for (_vpl = 0; _vpl < 19; _vpl++) {
		simd_store(_vpcurrentv[_vpl], &(line4[0]));
		for(_vpx1 = 0; _vpx1 < 4; ++_vpx1)
		{
			if(s_step_flags[_vpsis + _vpx1] == S_FLUID || s_step_flags[_vpsis + _vpx1] == S_BOUNCE)
			{
				s_step_current_line[(_vpsis + _vpx1) * 19 + _vpl] = line4[_vpx1];
			}
		}
	}
}

void s_simd_collide_step()
{
	
	for(_vpsis = 0; _vpsis < STEP_SIZE_MAIN; _vpsis+=4)
    {
        s_simd_collide_step_a();
        s_simd_collide_step_b();
        s_simd_collide_step_c();
        s_simd_collide_step_d();
        s_simd_collide_step_e();
        s_simd_collide_step_f();
        s_simd_collide_step_g();
        s_simd_collide_step_h();
        s_simd_collide_step_i();
	}
	
	s_collide_step(STEP_SIZE_MAIN, STEP_SIZE);
}

__thread_local int dbg_print_cnt;

__thread_local int dbg_x, dbg_y;

void s_stream_and_collide_xy(int nz,
                           realt* current_line,
                           int* walls,
                           int* flags,
                           realt* other_line)
{
    int k, l;

    int x_node_stride = (param.y_sec + 2) * param.Z * 19;
    int y_node_stride = param.Z * 19;

    int total_step = param.Z / STEP_SIZE;


    int cur_step;
    for(cur_step = 0; cur_step < total_step; cur_step++)
    {
        k = cur_step * STEP_SIZE;
        int z;

        stpf();
        seq_read_data(s_step_flags, flags + k, sizeof(int) * STEP_SIZE);
        seq_read_data(s_step_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        seq_read_data(s_step_current_line, current_line + k * 19, sizeof(realt) * STEP_SIZE * 19);
        edpf(CPF1_READ1);

        // memcpy(s_step_flags, flags + k, sizeof(int) * STEP_SIZE);
        // memcpy(s_step_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        // memcpy(s_step_current_line, current_line + k * 19, sizeof(int) * STEP_SIZE * 19);

        k = cur_step * STEP_SIZE;
        int z_st = k - 1, z_ed = k + STEP_SIZE;
        if(k == 0)
        {
            z_st = 0;
        }
        if(z_ed == nz)
        {
            z_ed--;
        }

        stpf();
        int ex, ey;
        for(ex = -1; ex <= 1; ex++)
        {
            for(ey = -1; ey <= 1; ey++)
            {
                seq_read_data(s_step_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                       &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                       sizeof(realt) * (z_ed - z_st + 1) * 19);
                // memcpy(s_step_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                //        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                //        sizeof(realt) * (z_ed - z_st + 1) * 19);
            }
        }
        edpf(CPF1_READ2);

        stpf();
        s_stream_step();
        edpf(CPF1_STREAM);
        
        stpf();
        //s_simd_collide_step();
        s_collide_step(0, STEP_SIZE);
        edpf(CPF1_COLLIDE);

        stpf();
        seq_write_data(current_line + k * 19, s_step_current_line, sizeof(realt) * STEP_SIZE * 19);
        edpf(CPF1_WRITE);
        // memcpy(current_line + k * 19, s_step_current_line, sizeof(realt) * STEP_SIZE * 19);
    }

}

void std_lbm()
{
    dbg_print_cnt = 0;
    int point_cnt = local_flag[STD_POINT_CNT];
    int* xs_ptr = (int*)local_flag[STD_XS_PTR];
    int* ys_ptr = (int*)local_flag[STD_YS_PTR];
    realt* current_head = (realt*)local_flag[STD_CURRENT_HEAD];
    realt* other_head = (realt*)local_flag[STD_OTHER_HEAD];

    int my_start, my_part, max_part;
    my_start = BLOCK_LOW(my_core, 64, point_cnt);
    my_part = BLOCK_SIZE(my_core, 64, point_cnt);
    max_part = BLOCK_SIZE(0, 64, point_cnt);

    int* fk_xs = xs_ptr;
    int* fk_ys = ys_ptr;

#define s_nodes_idx(d2, d3, d4, d5) ((((d2) * (param.y_sec + 2) + (d3)) * param.Z + (d4)) * 19 + (d5))
#define s_flags_idx(d2, d3, d4) (((d2) * (param.y_sec + 2) + (d3)) * param.Z + (d4))
#define s_walls_idx(d2, d3, d4, d5) ((((d2) * param.y_sec + (d3)) * param.Z + (d4)) * 19 + (d5))

    int i, j;
    for(i = my_start; i < my_start + my_part; i++)
    {
        int x = fk_xs[i];//GLD!!!
        int y = fk_ys[i];//GLD!!!
        dbg_x = x;
        dbg_y = y;

        //var for this point
        realt* point_current = current_head + s_nodes_idx(x, y, 0, 0);
        realt* point_other = other_head + s_nodes_idx(x, y, 0, 0);
        int* point_walls = param.walls + s_walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = param.flags + s_flags_idx(x, y, 0);

        s_stream_and_collide_xy(param.Z, point_current, point_walls, 
                              point_flags, point_other);
    }

    standard_write_flag(my_core);
}

void insane_lbm()
{

}


void cpe_athread_daemon(void *_param)
{
    my_core = athread_get_id(-1);

    int i, j, k;

    param = *((struct lbm_init_param*) _param);

    core_x = my_core % 8;
    core_y = my_core / 8;

    if(my_core == 0 && param.my_id == 0)
    {
        printf("@@@ MESSAGE FROM CPE @@@\n");
        printf("iter = %d, x_sec = %d, y_sec = %d, Z = %d\n", param.iter, param.x_sec, param.y_sec, param.Z);
        printf("val: nu = %f, omega = %f, CSmago = %f\n", param.nu, param.omega, param.CSmago);
        printf("ptrs: %x, %x, %x, %x\n",
               param.host_flag, param.slave_flag, param.walls, param.flags);
        printf("sizeof(int) = %d, sizeof(realt) = %d\n", sizeof(int), sizeof(realt));
        printf("@@@ MESSAGE FROM CPE @@@\n");
    }

    for (i = 0; i < FLAG_SIZE; i++)
    {
        local_flag[i] = 0; out_flag[i] = 0;
    }
    slave_flag_to_wait = 1;

    get_reply_god = 0, put_reply_god = 0;//the Jesus do not like waiting for someone
    get_reply_target = 0, put_reply_target = 0;

    //initial arr

    if(my_core == 0 && param.my_id == 0)
    {
        printf("CPE ENTER LOOP\n");
    }
    while(1)
    {
        stpf();
        standard_wait_flag(my_core);
        edpf(CPF1_WAIT);

        if (local_flag[KERNEL_ACTION] == EXIT_FLAG)
        {
            standard_write_flag(my_core);
            standard_write_flag(my_core);
            break;
        }
        if (local_flag[KERNEL_ACTION] == HANDSHAKE_FLAG)
        {
            if(my_core == 0 && param.my_id == 5)
            {
                printf("@@@ CPE HANDSHAKE_FLAG @@@\n");
                fflush(stdout);
            }
            standard_write_flag(my_core);
        }
        else if (local_flag[KERNEL_ACTION] == STD_LBM_FLAG)
        {
            if(my_core == 0 && param.my_id == 5)
            {
                printf("@@@ CPE STD_LBM_FLAG @@@\n");
                fflush(stdout);
            }
            std_lbm();
        }
        // else if(local_flag[KERNEL_ACTION] == INSANE_LBM_FLAG)
        // {
        //     insane_lbm();
        // }
        else
        {
            if(my_core == 0 && param.my_id == 0)
            {
                printf("unknown area!!!!!!!!!\n");
            }
        }
    }
}
