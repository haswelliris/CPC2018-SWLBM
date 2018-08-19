#include <stdio.h>
#include "slave.h"
#include "collide_slave.h"

__thread_local int my_id;
__thread_local volatile unsigned long get_reply,put_reply;

void athread_collide(void* paras)
{	
	struct athread_collide_parameter parameter = *(struct athread_collide_parameter*)paras;
	my_id = athread_get_id(-1);

	int Xst = parameter.Xst;
	int Xed = parameter.Xed;
	int Yst = parameter.Yst;
	int Yed = parameter.Yed;
	int nz = parameter.nz;
	int current = parameter.current;
	float nu = parameter.nu;
	float omega = parameter.omega;
	float CSmago = parameter.CSmago;
	int flags[500];
	float nodes[500][19];

	int i, j, k, l;
    float rho, u_x, u_y, u_z;
    float fi;
	float feq[19], Qo, omegaNew;


	for(i = Xst; i < Xed; i++) {
		if ((i-Xst)%64 != my_id)
			continue;
		for(j = Yst; j < Yed; j++) {
			get_reply = 0;
			athread_get(PE_MODE,&parameter.nodes[current][i-Xst+1][j-Yst+1][0][0],nodes,500*19*sizeof(float),&get_reply,0,0,0);
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
						const float tmp = (e[l][0] * u_x + e[l][1] * u_y + e[l][2] * u_z);
						feq[l] = ww[l] * rho * (1.0 -
							(3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
							(3.0 *     tmp) +
							(9.0/2.0 * tmp * tmp));
					}

					Qo=0;
					float Sij[3][3],S;
					// float e[19][3];
					int x1, y1, k1;
					// for(k1 = 0; k1 < 19; k1++) {
					// 	e[k1][0] = e_x[k1];
					// 	e[k1][1] = e_y[k1];
					// 	e[k1][2] = e_z[k1];
					// }

					for(x1 = 0; x1 < 3; x1++) {
						for(y1 = 0; y1 < 3; y1++) {
							Sij[x1][y1] = 0;
							for(k1 = 0; k1 < 19; k1++) {
								Sij[x1][y1] += e[k1][x1] * e[k1][y1] * (nodes[k][k1] - feq[k1]);
							}
							Qo += Sij[x1][y1] * Sij[x1][y1];
						}
					}

					nu = (2.0 / omega - 1.0) / 6.0;
					S = (-nu + sqrt(nu * nu + 18 * CSmago * CSmago * sqrt(Qo))) / 6.0 / CSmago / CSmago;
					omegaNew = 1.0 / (3.0 * (nu + CSmago * CSmago * S) + 0.5);

					for (l = 0; l < 19; l++) {
						nodes[k][l] =
							(1.0 - omegaNew) * nodes[k][l] +
							omegaNew * feq[l];
					}
				}
			}
			put_reply = 0;
			athread_put(PE_MODE,nodes,&parameter.nodes[current][i-Xst+1][j-Yst+1][0][0],500*19*sizeof(float),&put_reply,0,0);
			while(put_reply!=1);
		}
		printf(".");
	}
	printf("*");
}