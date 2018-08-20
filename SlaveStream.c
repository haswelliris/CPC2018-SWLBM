#include "Argument.h"

/*------------------------------------
 *      Main Computing Part
 *-----------------------------------*/
void slaveStream(Real *****nodes,
		    int ****walls,
			int ***flags,
			int Xst,
			int Xed,
			int Yst,
			int Yed,
			int nz,
		    int current,
			int other)
{
	int i, j, k, l;
	int inv;

	for(i = Xst+1; i < Xed-1; i++) {
		for(j = Yst+1; j < Yed-1; j++) {
			for(k = 0; k < nz; k++) {
				if(flags[i - Xst + 1][j - Yst + 1][k] & (1 << 20)) {
					for(l = 0; l < 19; l++) {
						inv = dfInv[l];
						nodes[current][i - Xst + 1][j - Yst + 1][k][l] = nodes[other][i - Xst + 1 + e_x[inv]][j - Yst + 1 + e_y[inv]][k + e_z[inv]][l];
					}
				}	
				if(flags[i - Xst + 1][j - Yst + 1][k] & (1 << 19)) {
					for(l = 0; l < 19; l++) {
						inv = dfInv[l];
						if(flags[i - Xst + 1][j - Yst + 1][k] & (1 << l)) {
							nodes[current][i - Xst + 1][j - Yst + 1][k][l] = nodes[other][i - Xst + 1][j - Yst + 1][k][inv];
						} else {
							nodes[current][i - Xst + 1][j - Yst + 1][k][l] = nodes[other][i - Xst + 1 + e_x[inv]][j - Yst + 1 + e_y[inv]][k + e_z[inv]][l];
						}
					}
				}	
			}
		}
	}
}
