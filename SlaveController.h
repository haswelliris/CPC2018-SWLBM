#ifndef __SLAVECONTROLLER_H__
#define __SLAVECONTROLLER_H__

struct athread_init_parameter
{
	float *****nodes;
    int ***flags;
    int ****walls;
	int Xst;
	int Xed;
	int Yst;
	int Yed;
	int nz;
	float omega;
	float CSmago;
	int my2drank;
};

#endif
