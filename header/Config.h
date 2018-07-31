#ifndef __LBM_CONFIG_H__
#define __LBM_CONFIG_H__

#define NUM_DIMS 2
#define NUM_THREADS 4
#define OUTNAME "result/outlbm"

/*-----------------------------*
 *  different types of cells
 *-----------------------------*/
#define FLUID 0
#define	NOSLIP 1
#define VELOCITY 2
#define BOUNCE 3
#define	PRESSURE 4

typedef float Real;

#endif
