#include <stdio.h>
#include <math.h>
#include "slave.h"
#include "SlaveController.h"

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

__thread_local int my_id, my2drank;
__thread_local int Xst, Xed, Yst, Yed, nz, current;
__thread_local float athread_nu, athread_omega, athread_CSmago;
__thread_local int ***flags;
__thread_local int tmp_flags[500];
__thread_local float *****nodes;
__thread_local float tmp_nodes[500][19];
__thread_local int i, j, k, l;
__thread_local float rho, u_x, u_y, u_z;
__thread_local float fi;
__thread_local float feq[19], Qo, omegaNew;

__thread_local float tmp;
__thread_local float Sij[3][3],S;
__thread_local int x_1, y_1, k_1;

__thread_local volatile unsigned long get_reply,put_reply;
__thread_local volatile int total_mask = 0xFFFFFFFF, lpc, trump_i, trump_j;

void slaveInit(void* paras) {
	struct athread_init_parameter parameter = *(struct athread_init_parameter*)paras;
	my_id = athread_get_id(-1);
	my2drank = parameter.my2drank;

	Xst = parameter.Xst;
	Xed = parameter.Xed;
	Yst = parameter.Yst;
	Yed = parameter.Yed;
	nz = parameter.nz;
	flags = parameter.flags;
	nodes = parameter.nodes;
	athread_omega = parameter.omega;
	athread_CSmago = parameter.CSmago;
}

void slaveController(void* paras) {
	// slaveInit(paras);
	
	// while(true) {
	// 	while (noSignal());
	// 	if (signal == SLAVE_CALCULATE) {
	// 		slaveStream();
	// 		slaveCollide();
	// 	}
	// 	else if (signal == SLAVE_EXIT) {
	// 		return;
	// 	}
	// 	else {
	// 		printf("slaveController signal wrong\n");
	// 	}
	// }
}

void SlaveCollide(void* paras)
{
	current = *(int*)paras;
	for(i = Xst+1; i < Xed-1; i++) {
		for(j = Yst+1; j < Yed-1; j++) {
			if (((i-Xst-1)*(Yed-Yst-2)+(j-Yst-1))%64 != my_id)
				continue;

			get_reply = 0;
			athread_get(PE_MODE,&nodes[current][i-Xst+1][j-Yst+1][0][0],&tmp_nodes[0][0],500*19*sizeof(float),&get_reply,0,0,0);
			athread_get(PE_MODE,&flags[i-Xst+1][j-Yst+1][0],tmp_flags,500*sizeof(int),&get_reply,0,0,0);
			while(get_reply != 2);
			for(k = 0; k < nz; k++) {
				if(tmp_flags[k] == 0 || tmp_flags[k] == 3) {
					rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
					for (l = 0; l < 19; l++) {
						fi = tmp_nodes[k][l];
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
								Sij[x_1][y_1] += e[k_1][x_1] * e[k_1][y_1] * (tmp_nodes[k][k_1] - feq[k_1]);
							}
							Qo += Sij[x_1][y_1] * Sij[x_1][y_1];
						}
					}

					athread_nu = (2.0 / athread_omega - 1.0) / 6.0;
					S = (-athread_nu + sqrt(athread_nu * athread_nu + 18 * athread_CSmago * athread_CSmago * sqrt(Qo))) / 6.0 / athread_CSmago / athread_CSmago;
					omegaNew = 1.0 / (3.0 * (athread_nu + athread_CSmago * athread_CSmago * S) + 0.5);

					for (l = 0; l < 19; l++) {
						tmp_nodes[k][l] =
							(1.0 - omegaNew) * tmp_nodes[k][l] +
							omegaNew * feq[l];
					}
				}
			}
			put_reply = 0;
			athread_put(PE_MODE,&tmp_nodes[0][0],&nodes[current][i-Xst+1][j-Yst+1][0][0],500*19*sizeof(float),&put_reply,0,0);
			while(put_reply!=1);
		}
	}
	// 有没有这个athread=0最后再写入一下结果都一样
	// athread_syn(ARRAY_SCOPE, total_mask);
	// if (my_id == 0) {
	// 	get_reply = 0;
	// 	athread_get(PE_MODE,&tmp_flags[1][1][0],tmp_flags,500*sizeof(int),&get_reply,0,0,0);
	// 	while(get_reply != 0);

	// 	put_reply = 0;
	// 	athread_put(PE_MODE,tmp_flags,&tmp_flags[1][1][0],500*sizeof(int),&put_reply,0,0);
	// 	while(put_reply!=1);

	// }
	// printf("*");
}