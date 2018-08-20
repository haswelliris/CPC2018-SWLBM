#include <stdio.h>

#define Real long

Real nodes[2][252][127][500][19];

#define encode(i,j,k,l) (i*127*500*19+j*500*19+k*19+l)
#define decode(n) printf("%ld %ld %ld %ld", n/(127*500*19), n%(127*500*19)/(500*19), n%(500*19)/19, n%19)

int Xst = 500;
int Xed = 750;
int Yst = 250;
int Yed = 375;
int nz = 500;

void slaveStream(int* paras) {
	int current, other;
	current = *(int*)paras;
	other = 1 - current;



}

int main() {
	long long int i, j, k, l;
	for (i = 0; i < 252; i++)
		for (j = 0; j < 127; j++)
			for (k = 0; k < 500; k++)
				for (l = 0; l < 19; l++)
					nodes[0][i][j][k][l] = encode(i,j,k,l);
	// decode(nodes[251][126][499][18]);

	int in = 1;
	slaveStream(&in);

	return 0;
}