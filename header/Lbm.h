#ifndef __LBM_LBM_H__
#define __LBM_LBM_H__

#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <string.h>

#include "Config.h"

void TIME_ST();

void TIME_ED();

void setParameter();

Real **array2DF(int row, int col);

void arrayFree2DF(Real **array);

int **array2DI(int row, int col);

void arrayFree2DI(int **array);

Real ***array3DF(int row, int col, int channel);

void arrayFree3DF(Real ***array);

int ***array3DI(int row, int col, int channel);

void arrayFree3DI(int ***array);

int ****array4DI(int row1, int row, int col, int channel);

void arrayFree4DI(int ****array);

Real *****array5DF(int row2, int row1, int row, int col, int channel);

void arrayFree5DF(Real *****array);

void write_result(int X, int Y, int Z, int s, double times, Real **image);

void init_Flag_local(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int ***flags);

void collide(Real *****nodes, int ***flags, int Xst, int Xed, int Yst, int Yed, int nz, int current);

void writeImage_tecplot(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int other, int myrank, Real **local_image, Real *****nodes);

void init_Pointer(int ***flags, Real *****nodes, int ****walls, int Xst, int Xed, int Yst, int Yed, int nz,int x_sec, int y_sec, Real rho,Real* u);

void stream(Real *****nodes, int ****walls, int ***flags, int Xst, int Xed, int Yst, int Yed, int nz, int current, int other);

void SetMPI(MPI_Comm mycomm, int *dims, int *coords);

void INITINPUT(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int x_sec, int y_sec, 
		       int myrank, int size, char *user, int *local_rankinfo, int **rankinfo, int ***flags);

int OUTPUT(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int s, 
		   int myrank, int size, int other, int x_sec, int y_sec, char *user, 
		   Real **local_image, Real **image, int **rankinfo, Real *****nodes);

void bounce_send_init_test (int X, int Y, int Z,
							int Xst, int Xed, int Yst, int Yed, 
							int x_sec, int y_sec, int corrent, int other, 
							Real *****nodes,
						    Real ***temp_left_send,
							Real ***temp_right_send, 
							Real ***temp_up_send, 
							Real ***temp_down_send,
						    Real **temp_ld_send, 
							Real **temp_lu_send, 
							Real **temp_rd_send, 
							Real **temp_ru_send);

void bounce_communicate_Itest(MPI_Comm mycomm, int *dims, int *coords, 
							  int x_sec, int y_sec, int Z, int *l_count,
							  MPI_Status *sta, MPI_Request *req,
							  Real ***temp_left_send, 
							  Real ***temp_right_send, 
							  Real ***temp_up_send, 
							  Real ***temp_down_send, 
							  Real ***temp_left, 
							  Real ***temp_right, 
							  Real ***temp_up, 
							  Real ***temp_down,
							  Real **temp_lu_send, 
							  Real **temp_ld_send, 
							  Real **temp_ru_send, 
							  Real **temp_rd_send, 
							  Real **temp_lu, 
							  Real **temp_ld, 
							  Real **temp_ru, 
							  Real **temp_rd);

void bounce_update_test(int X, int Y, int Z,
						int Xst, int Xed, int Yst, int Yed, 
						int myrank, int x_sec, int y_sec,
						int other,
						Real *****nodes,
					    Real ***temp_left, 
						Real ***temp_right, 
						Real ***temp_up, 
						Real ***temp_down,
					    Real **temp_ld, 
						Real **temp_lu, 
						Real **temp_rd, 
						Real **temp_ru);

#endif
