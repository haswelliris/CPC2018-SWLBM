#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Argument.h"

int main(int argc, char *argv[])                                       
{                                                                      
    int ***flags;                                                      
    int i, j, k, s, l;
    int myrank, my2drank, size;
    int current = 0, other = 1; 
    int periods[NUM_DIMS] = {0, 0};
    int dims[NUM_DIMS] = {0, 0};
    int coords[2];
    Real df;
    unsigned int X_section ;
    unsigned int X_res     ;
    unsigned int Y_section ;
    unsigned int Y_res     ;
    int	Xst ;
    int	Xed ;
    int	Yst ;
    int	Yed ;
    int x_sec ;
    int y_sec ;
    Real *****nodes;
    int ****walls;
    Real ***temp_right;
    Real ***temp_left;
    Real ***temp_down;
    Real ***temp_up;
    Real ***temp_right_send;
    Real ***temp_left_send;
    Real ***temp_down_send;
    Real ***temp_up_send;
    Real **temp_lu;
    Real **temp_ld;
    Real **temp_ru;
    Real **temp_rd;
    Real **temp_lu_send;
    Real **temp_ld_send;
    Real **temp_ru_send;
    Real **temp_rd_send;
    MPI_Comm mycomm;
    MPI_Status sta[16];
    MPI_Request req[16];
    int n = 0;
    int count;
    int local_rankinfo[7];
    int **rankinfo;
    Real **image;
    Real **local_image;
		
	/*------------------------
         * Parameter Set
         * ----------------------*/

    	setParameter();
	
	/*---------------------------*
	 * MPI Init
	 * --------------------------*/
    

    	MPI_Init(&argc, &argv);

    	MPI_Comm_size(MPI_COMM_WORLD, &size);

    	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	MLOG("CPC -- LBM-simulation Benchmark!\n\n");
	 
	/*----------------------------*
         * MPI Mesh Section
         * ---------------------------*/
    
	MPI_Dims_create(size, NUM_DIMS, dims);

    	MPI_Cart_create(MPI_COMM_WORLD, NUM_DIMS, dims, periods, 0, &mycomm);

    	MPI_Comm_rank(mycomm, &my2drank);

    	MPI_Cart_coords(mycomm, my2drank, 2, coords);

	SetMPI(mycomm, dims, coords);
	


    	X_section = X / dims[0];
    	X_res     = X % dims[0];
    	Y_section = Y / dims[1];
    	Y_res     = Y % dims[1];
    
	Xst = (coords[0] < X_res) ? (coords[0] * (X_section + 1))
								  : (coords[0] *  X_section + X_res);
	Xed = (coords[0] < X_res) ? (Xst       + (X_section + 1))
								  : (Xst       +  X_section);
	Yst = (coords[1] < Y_res) ? (coords[1] * (Y_section + 1))
								  : (coords[1] *  Y_section + Y_res);
	Yed = (coords[1] < Y_res) ? (Yst       + (Y_section + 1))
								  : (Yst       +  Y_section);
	
    	x_sec = Xed - Xst;
    	y_sec = Yed - Yst;
	

	/*-----------------------------*
         *-----------------------------* 
         * ----------------------------*/
	
	MLOG("Size                 :  %d x %d x %d\n", X, Y, Z);
        MLOG("Steps                :  %d\n", STEPS);
        MLOG("Number of Process    :  %d\n\n", size);
	
        /*-----------------------------*
         * Space Allocate
         * ----------------------------*/

	

	flags = array3DI(x_sec + 2, y_sec + 2, Z);
	nodes = array5DF(2, x_sec + 2, y_sec + 2, Z, 19);
	walls = array4DI(x_sec, y_sec, Z, 19);
    	temp_lu_send    = array2DF(Z, 19);
    	temp_lu         = array2DF(Z, 19);
    	temp_ld         = array2DF(Z, 19);
    	temp_ru         = array2DF(Z, 19);
    	temp_rd         = array2DF(Z, 19);
    	temp_ld_send    = array2DF(Z, 19);
    	temp_ru_send    = array2DF(Z, 19);
    	temp_rd_send    = array2DF(Z, 19);
    	temp_right      = array3DF(y_sec, Z, 19);
    	temp_left       = array3DF(y_sec, Z, 19);
    	temp_down       = array3DF(x_sec, Z, 19);
    	temp_up         = array3DF(x_sec, Z, 19);
    	temp_right_send = array3DF(y_sec, Z, 19);
    	temp_left_send  = array3DF(y_sec, Z, 19);
    	temp_down_send  = array3DF(x_sec, Z, 19);
    	temp_up_send    = array3DF(x_sec, Z, 19);
    	memset(&walls[0][0][0][0], 0, x_sec * y_sec * Z * 19 * sizeof(int));
	
	/*----------------------------------------*
	 * local rank info
	 * ---------------------------------------*/

	local_image = array2DF(y_sec, x_sec * 5);

	if(my2drank == 0) {
		image = array2DF(Y, X * 5);
		rankinfo = array2DI(size, 7);
	} else {
		image = array2DF(1, 1);
		rankinfo = array2DI(1, 1);
	}

	INITINPUT(X, Y, Z, Xst, Xed, Yst, Yed, x_sec, y_sec, myrank, size, argv[1], local_rankinfo, rankinfo, flags);

	MLOG("Init >> init Pointer!\n");

    	init_Pointer(flags, nodes, walls, Xst, Xed, Yst, Yed, Z, x_sec, y_sec, 1, LDC_VELOCITY);

	MLOG("Step >> Main Steps start!\n\n");
	

	/*----------------------------------------------------*
	 * Main Calculation section
	 * ---------------------------------------------------*/
	TIME_ST();

	masterController(
			flags,
		    myrank,
		    my2drank,
		    current,
		    other, 
		    dims,
		    coords,
		    Xst,
		    Xed,
		    Yst,
		    Yed,
		    x_sec,
		    y_sec,
		    nodes,
		    walls,
		    temp_right,
		    temp_left,
		    temp_down,
		    temp_up,
		    temp_right_send,
		    temp_left_send,
		    temp_down_send,
		    temp_up_send,
		    temp_lu,
		    temp_ld,
		    temp_ru,
		    temp_rd,
		    temp_lu_send,
		    temp_ld_send,
		    temp_ru_send,
		    temp_rd_send,
		    mycomm,
		    sta,
		    req,
		    count
		);
		s = STEPS;

	TIME_ED();
	/*-----------------------------*
 	 * OUTPUT 
 	 *-----------------------------*/

	MLOG("Step >> Main Steps Done!\n\n");

	OUTPUT(X, Y, Z, Xst, Xed, Yst, Yed, s, myrank, size, other, x_sec, y_sec, argv[1], local_image, image, rankinfo, nodes);

	// Check Result
	#define CHECK_RESULT 0
	#define SAVE_RESULT 1
	int dump = CHECK_RESULT;
	FILE *outfile, *infile;
	char fname[100];
	sprintf(fname, "result/lbm_steps-%d_dump_rank-%02d.dat", STEPS, myrank);

	if(dump == SAVE_RESULT) 
	{
		outfile = fopen(fname, "wb");
		if(outfile == NULL)
        {
			MLOG("EMMMM!? cannot write\n");
		}
		else
		{
			MLOG("Writing\n");
			int i0, i1, i2, i3;
            for(i0 = 0; i0 < 2; i0++)
            {
                for(i1 = 1; i1 < x_sec + 1; i1++)
                {
					i2 = 1;
                    for(i3 = 0; i3 < Z; i3++)
                    {
                        fwrite(nodes[i0][i1][i2][i3], sizeof(Real), 1, outfile);
                    }
                }
            }
            fclose(outfile);
		}
	}
	else
	{
		infile = fopen(fname, "rb");
		if(infile == NULL) {
			MLOG("Please dump the data first :3\n");
		}
		else
		{
			MLOG("Checking\n");
			if(1)
			{
            	Real* dumped = malloc(sizeof(Real) * 2 * (x_sec + 2) * 1 * Z * 1);
				fread(dumped, sizeof(Real), 2 * (x_sec + 2) * 1 * Z * 1, infile);

				int i0, i1, i2, i3;
				int iter = 0;
				int error_count = 0;
				for(i0 = 0; i0 < 2; i0++)
				{
					for(i1 = 1; i1 < x_sec + 1; i1++)
					{
						i2 = 1;
						for(i3 = 0; i3 < Z; i3++)
						{
							if(isnan(dumped[iter]) && isnan(nodes[i0][i1][i2][i3][0]))
							{
							}
							else if(isnan(dumped[iter]))
							{
							}
							else if(isnan(nodes[i0][i1][i2][i3][0]))
							{
								OLOG(my2drank, "error data at (%d %d %d %d %d) is nan, at [%d %d]\n", i0, i1, i2, i3, 0, coords[0], coords[1]);
								iter = -1;
								break;
							}
							else if(dumped[iter] - nodes[i0][i1][i2][i3][0] > 1e-6)
							{
								if (error_count == 0)
									OLOG(my2drank, "error data at [%d %d %d %d %d]: %f %f\n", i0, i1, i2, i3, 0, dumped[iter], nodes[i0][i1][i2][i3][0]);
								// iter = -1;
								error_count++;
								break;
							}
							iter++;
						}
						if(iter == -1)
							break;
					}
					if(iter == -1)
						break;
				}
				if (error_count != 0)
					OLOG(my2drank,"error_count = %d\n", error_count);
			}
		}
	}


	arrayFree2DF(image);
	arrayFree2DI(rankinfo);
	arrayFree2DF(local_image);

	arrayFree3DI(flags);
	arrayFree5DF(nodes);
	arrayFree4DI(walls);
	arrayFree3DF(temp_right);
	arrayFree3DF(temp_left);
	arrayFree3DF(temp_down);
	arrayFree3DF(temp_up);
	arrayFree2DF(temp_lu);
	arrayFree2DF(temp_ld);
	arrayFree2DF(temp_ru);
	arrayFree2DF(temp_rd);
	arrayFree3DF(temp_right_send);
	arrayFree3DF(temp_left_send);
	arrayFree3DF(temp_down_send);
	arrayFree3DF(temp_up_send);
	arrayFree2DF(temp_lu_send);
	arrayFree2DF(temp_ld_send);
	arrayFree2DF(temp_ru_send);
	arrayFree2DF(temp_rd_send);

	MLOG("LBM-simulation Done!\n");

    	MPI_Finalize();

	return 0;
}
