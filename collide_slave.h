#ifndef __Collide_Slave_H__
#define __Collide_Slave_H__
// const int e_x[19] = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0};
// const int e_y[19] = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0};
// const int e_z[19] = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0};


struct athread_collide_parameter
{
	float *****nodes;
    int ***flags;
	int Xst;
	int Xed;
	int Yst;
	int Yed;
	int nz;
	int current;
	float nu;
	float omega;
	float CSmago;
};

#endif
