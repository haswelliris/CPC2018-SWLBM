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
	{ 0, 0, 0}
};

__thread_local const float ww[19] = {
		(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };

__thread_local const int ddfInv[19] = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};

__thread_local const int tl[9][3] = {
	{-1,-1, 8},
	{12,13, 1},
	{-1,-1, 9},
	{14,15, 2},
	{ 4, 5,18},
	{16,17, 3},
	{-1,-1, 6},
	{10,11, 0},
	{-1,-1, 7}
};

__thread_local const int tk[3] = {-1, 1, 0};

__thread_local int my_id, my2drank;
__thread_local int Xst, Xed, Yst, Yed, nz, current, other;
__thread_local float athread_nu, athread_omega, athread_CSmago;
__thread_local int ***flags;

__thread_local float *****nodes;
__thread_local int tmp_walls[64][19];
__thread_local int i, j, k, kk, l, r, xi, yi;
__thread_local float rho, u_x, u_y, u_z;
__thread_local float fi;
__thread_local float feq[19], Qo, omegaNew;

__thread_local float tmp;
__thread_local float Sij[3][3],S;

__thread_local int x_1, y_1, k_1;

__thread_local int now, next;
__thread_local float result_nodes[2][21][19], buffer_nodes[2][3][3][22][19];
__thread_local int tmp_flags[2][21];
__thread_local volatile unsigned long result_nodes_put_reply[2], buffer_nodes_get_reply[2];
__thread_local const int xy_around[9][2] = {
	{-1, 1},
	{ 0, 1},
	{ 1, 1},
	{-1, 0},
	{ 0, 0},
	{ 1, 0},
	{-1,-1},
	{ 0,-1},
	{ 1,-1}
};

__thread_local int inv;

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

inline void load_buffer(int target, int x, int y, int z) {
	buffer_nodes_get_reply[target] = 0;
	athread_get(PE_MODE,&flags[x][y][z],&tmp_flags[target][1],sizeof(int)*20,&buffer_nodes_get_reply[target],0,0,0);
	for (xi = 0; xi < 3; xi++)
		for (yi = 0; yi < 3; yi++) {
			if (z == 0)
				athread_get(PE_MODE,&nodes[other][x+xi-1][y+yi-1][ z ][0],&buffer_nodes[target][xi][yi][1][0],sizeof(float)*21*19,&buffer_nodes_get_reply[target],0,0,0);
			else if (z == 480)
				athread_get(PE_MODE,&nodes[other][x+xi-1][y+yi-1][z-1][0],&buffer_nodes[target][xi][yi][0][0],sizeof(float)*21*19,&buffer_nodes_get_reply[target],0,0,0);
			else
				athread_get(PE_MODE,&nodes[other][x+xi-1][y+yi-1][z-1][0],&buffer_nodes[target][xi][yi][0][0],sizeof(float)*22*19,&buffer_nodes_get_reply[target],0,0,0);
		}
}

inline void stream(int target) {
	for(kk = 1; kk <= 20; kk++) {
		if(tmp_flags[target][kk] & (1 << 20)) {
			for(l = 0; l < 19; l++) {
				inv = ddfInv[l];
				result_nodes[target][kk][l] = buffer_nodes[target][e[inv][0]+1][e[inv][1]+1][kk + e[inv][2]][l];
			}
		}	
		if(tmp_flags[target][kk] & (1 << 19)) {
			for(l = 0; l < 19; l++) {
				inv = ddfInv[l];
				if(tmp_flags[target][kk] & (1 << l)) {
					result_nodes[target][kk][l] = buffer_nodes[target][1][1][kk][inv];
				} else {
					result_nodes[target][kk][l] = buffer_nodes[target][e[inv][0]+1][e[inv][1]+1][kk + e[inv][2]][l];
				}
			}
		}	
	}
}

inline void collide(int target) {
	for(kk = 1; kk <= 20; kk++) {
		if(tmp_flags[target][kk] & (1<<20) || tmp_flags[target][kk] & (1<<19)) {
			rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
			for (l = 0; l < 19; l++) {
				fi = result_nodes[target][kk][l];
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
						Sij[x_1][y_1] += e[k_1][x_1] * e[k_1][y_1] * (result_nodes[target][kk][k_1] - feq[k_1]);
					}
					Qo += Sij[x_1][y_1] * Sij[x_1][y_1];
				}
			}

			athread_nu = (2.0 / athread_omega - 1.0) / 6.0;
			S = (-athread_nu + sqrt(athread_nu * athread_nu + 18 * athread_CSmago * athread_CSmago * sqrt(Qo))) / 6.0 / athread_CSmago / athread_CSmago;
			omegaNew = 1.0 / (3.0 * (athread_nu + athread_CSmago * athread_CSmago * S) + 0.5);

			for (l = 0; l < 19; l++) {
				result_nodes[target][kk][l] = (1.0 - omegaNew) * result_nodes[target][kk][l] + omegaNew * feq[l];
			}
		}
	}
}

inline void store_nodes(int target, int x, int y, int z) {
	result_nodes_put_reply[target] = 0;
	athread_put(PE_MODE,&result_nodes[target][1][0],&nodes[current][x][y][z][0],20*19*sizeof(float),&result_nodes_put_reply[target],0,0);
}

void SlaveCollide(void* paras) {
	current = *(int*)paras;
	other = 1 - current;
	for(i = Xst+1; i < Xed-1; i++) {
		for(j = Yst+1; j < Yed-1; j++) {
			if (((i-Xst-1)*(Yed-Yst-2)+(j-Yst-1))%64 != my_id)
				continue;
			now = 0;
			next = 1;
			result_nodes_put_reply[now] = 1;
			result_nodes_put_reply[next] = 1;
			//load buffer now
			load_buffer(now, i-Xst+1, j-Yst+1, 0);
			for (k = 0; k < nz; k+=20) {
				// start load buffer next
				if (k != 480)
					load_buffer(next, i-Xst+1, j-Yst+1, k+20);
				// end store result now
				// while (result_nodes_put_reply[now] != 1);
				// end load buffer now
				while (buffer_nodes_get_reply[now] != 10);
				// stream nodes now
				stream(now);
				// collide nodes now
				// collide(now);
				// start store nodes now
				store_nodes(now, i-Xst+1, j-Yst+1, k);
				// reverse now next
				now = 1 - now;
				next = 1 - next;
			}
			// while (result_nodes_put_reply[next] != 1);
		}
	}
}
