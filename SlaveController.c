#include <stdio.h>
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

void slaveInit(void* paras) {
	struct athread_init_parameter parameter = *(struct athread_init_parameter*)paras;
	my_id = athread_get_id(-1);
	my2drank = parameter.my2drank;

	Xst = parameter.Xst;
	Xed = parameter.Xed;
	Yst = parameter.Yst;
	Yed = parameter.Yed;
	nz = parameter.nz;
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
