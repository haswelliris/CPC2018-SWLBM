#include <stdio.h>
#include <math.h>
#include "slave.h"
#include "collide_slave.h"

__thread_local const int e[19][3] = {
	{ 0, 1, 0},
	{ 0,-1, 0},
	{ 1, 0, 0},
	{-1, 0, 0},
	{ 0, 0, 1},
	{ 0, 0,-1},
	{ 1, 1, 0},
	{-1, 1, 0},
	{ 1,-1, 0},
	{-1,-1, 0},
	{ 0, 1, 1},
	{ 0, 1,-1},
	{ 0,-1, 1},
	{ 0,-1,-1},
	{ 1, 0, 1},
	{ 1, 0,-1},
	{-1, 0, 1},
	{-1, 0,-1},
	{0, 0, 0,}
};

__thread_local const float ww[19] = {
		(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };

__thread_local struct athread_collide_parameter parameter;
__thread_local int my_id;
__thread_local int Xst, Xed, Yst, Yed, nz, current;
__thread_local float athread_nu, athread_omega, athread_CSmago;
__thread_local int flags[500];
__thread_local float nodes[500][19];
__thread_local int i, j, k, l;
__thread_local float rho, u_x, u_y, u_z;
__thread_local float fi;
__thread_local float feq[19], Qo, omegaNew;

__thread_local float tmp;
__thread_local float Sij[3][3],S;
__thread_local int x_1, y_1, k_1;

__thread_local volatile unsigned long get_reply,put_reply;
__thread_local volatile int total_mask = 0xFFFFFFFF, lpc, trump_i, trump_j;

void athread_collide(void* paras)
{	
	parameter = *(struct athread_collide_parameter*)paras;
	my_id = athread_get_id(-1);

	Xst = parameter.Xst;
	Xed = parameter.Xed;
	Yst = parameter.Yst;
	Yed = parameter.Yed;
	nz = parameter.nz;
	current = parameter.current;
	athread_nu = parameter.nu;
	athread_omega = parameter.omega;
	athread_CSmago = parameter.CSmago;

	for(i = Xst; i < Xed; i++) {
		if ((i-Xst)%64 != my_id)
			continue;
		for(j = Yst; j < Yed; j++) {
			get_reply = 0;
			athread_get(PE_MODE,&parameter.nodes[current][i-Xst+1][j-Yst+1][0][0],&nodes[0][0],500*19*sizeof(float),&get_reply,0,0,0);
			athread_get(PE_MODE,&parameter.flags[i-Xst+1][j-Yst+1][0],flags,500*sizeof(int),&get_reply,0,0,0);
			while(get_reply != 2);
			for(k = 0; k < nz; k++) {
				if(flags[k] == 0 || flags[k] == 3) {
					rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
					for (l = 0; l < 19; l++) {
						fi = nodes[k][l];
						rho += fi;
						u_x += e[l][0] * fi;
						u_y += e[l][1] * fi;
						u_z += e[l][2] * fi;
					}

					u_x /= rho;
                    u_y /= rho;
                    u_z /= rho;

					for (l = 0; l < 19; l++) {
						tmp = (e[l][0] * u_x + e[l][1] * u_y + e[l][2] * u_z);
						feq[l] = ww[l] * rho * (1.0 -
							(3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
							(3.0 *     tmp) +
							(9.0/2.0 * tmp * tmp));
					}

					Qo=0;
					for(x_1 = 0; x_1 < 3; x_1++) {
						for(y_1 = 0; y_1 < 3; y_1++) {
							Sij[x_1][y_1] = 0;
							for(k_1 = 0; k_1 < 19; k_1++) {
								Sij[x_1][y_1] += e[k_1][x_1] * e[k_1][y_1] * (nodes[k][k_1] - feq[k_1]);
							}
							Qo += Sij[x_1][y_1] * Sij[x_1][y_1];
						}
					}

					athread_nu = (2.0 / athread_omega - 1.0) / 6.0;
					S = (-athread_nu + sqrt(athread_nu * athread_nu + 18 * athread_CSmago * athread_CSmago * sqrt(Qo))) / 6.0 / athread_CSmago / athread_CSmago;
					omegaNew = 1.0 / (3.0 * (athread_nu + athread_CSmago * athread_CSmago * S) + 0.5);

					for (l = 0; l < 19; l++) {
						nodes[k][l] =
							(1.0 - omegaNew) * nodes[k][l] +
							omegaNew * feq[l];
					}
				}
			}
			put_reply = 0;
			athread_put(PE_MODE,&nodes[0][0],&parameter.nodes[current][i-Xst+1][j-Yst+1][0][0],500*19*sizeof(float),&put_reply,0,0);
			while(put_reply!=1);
		}
	}
	// 有没有这个athread=0最后再写入一下结果都一样
	// athread_syn(ARRAY_SCOPE, total_mask);
	// if (my_id == 0) {
	// 	get_reply = 0;
	// 	athread_get(PE_MODE,&parameter.flags[1][1][0],flags,500*sizeof(int),&get_reply,0,0,0);
	// 	while(get_reply != 0);

	// 	put_reply = 0;
	// 	athread_put(PE_MODE,flags,&parameter.flags[1][1][0],500*sizeof(int),&put_reply,0,0);
	// 	while(put_reply!=1);

	// }
	// printf("*");
}