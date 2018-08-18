#include "../header/Argument.h"
#include "../header/Array.h"

//#define nodes_idx(d1, d2, d3, d4, d5) (((((d1) * (x_sec + 2) + (d2)) * (y_sec + 2) + (d3)) * nz + (d4)) * 19 + (d5))

/*------------------------------------
 *      Main Computing Part
 *-----------------------------------*/
void stream(Real *****nodes,
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
	int x_sec = Xed-Xst;
	int y_sec = Yed-Yst;

	for(i = Xst; i < Xed; i++) {
		for(j = Yst; j < Yed; j++) {
			for(k = 0; k < nz; k++) {
				if(flags[i - Xst + 1][j - Yst + 1][k] == FLUID) {
					for(l = 0; l < 19; l++) {
						inv = dfInv[l];
						//(&nodes[0][0][0][0][0])[nodes_idx(current, i - Xst + 1, j - Yst + 1, k, l)] = nodes[other][i - Xst + 1 + e_x[inv]][j - Yst + 1 + e_y[inv]][k + e_z[inv]][l];
						Index5D(&nodes[0][0][0][0][0], x_sec + 2, y_sec + 2, nz, 19, current, i - Xst + 1, j - Yst + 1, k, l) = nodes[other][i - Xst + 1 + e_x[inv]][j - Yst + 1 + e_y[inv]][k + e_z[inv]][l];
						nodes[current][i - Xst + 1][j - Yst + 1][k][l] = nodes[other][i - Xst + 1 + e_x[inv]][j - Yst + 1 + e_y[inv]][k + e_z[inv]][l];
					}
				}	
				if(flags[i - Xst + 1][j - Yst + 1][k] == BOUNCE) {
					for(l = 0; l < 19; l++) {
						inv = dfInv[l];
						if(walls[i - Xst][j - Yst][k][l]) {
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
