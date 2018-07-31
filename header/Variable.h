#ifndef __LBM_VARIABLE_H__
#define __LBM_VARIABLE_H__

#include <stdio.h>

#include "Config.h"

int X, Y, Z, STEPS;

double t1, t2;
unsigned long r1, r2;

/*------------------------------------------*
 * Declare some constants and globals ...
 * -----------------------------------------*/
int v_log;
Real CSmago;
Real LDC_VELOCITY[3], V0[3] = {0.0};
Real RE;
Real R;
Real L_cylinder;
Real Rmax_suboff;
Real m_per_ft;
Real ox_suboff,oy_suboff,oz_suboff;
Real rho0;
Real nu;
Real dx, dt;
Real pressureScale, refP;
Real omega;     // viscosity of the fluid, 0..2
Real exp_vel;

int tmp_coords[2];
int right_nbr = -1,
    left_nbr  = -1,
    up_nbr    = -1,
    down_nbr  = -1,
    lu_nbr    = -1,
    ru_nbr    = -1,
    ld_nbr    = -1,
    rd_nbr    = -1;

/*----------------------------------------------------*
 *   the weight for the equilibrium distribution
 *   -------------------------------------------------*/
const Real w[19] = {
		(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };

/*-----------------------------------------------------*
 *     convenience variables for handling
 *     the distribution functions in a loop
 *     the 19 lattice vectors
 *-----------------------------------------------------*/
const int e_x[19] = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0};
const int e_y[19] = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0};
const int e_z[19] = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0};

/*-----------------------------------------------------*
 * the index of the inverse for each lattice vector
 * ----------------------------------------------------*/
const int dfInv[19] = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};

typedef struct ImageInfo
{
	int x;
	int y;
	Real u_x;
	Real u_y;
	Real u_z;
	Real MagV;
	Real rho;
} ImageInfo;

static inline unsigned long rpcc()
{
     unsigned long time;
     asm("rtc %0": "=r" (time) : );
     return time;
}

#endif
