#include <athread.h>
#include "Argument.h"
#include "SlaveController.h"

// #define SMALL_STEP 2

extern SLAVE_FUN(slaveController)();
extern SLAVE_FUN(slaveInit)();
extern SLAVE_FUN(SlaveCollide)();

void masterController(
			int ***flags,
		    int myrank,
		    int my2drank,
		    int current,
		    int other, 
		    int *dims,
		    int *coords,
		    int	Xst ,
		    int	Xed ,
		    int	Yst ,
		    int	Yed ,
		    int x_sec ,
		    int y_sec ,
		    Real *****nodes,
		    int ****walls,
		    Real ***temp_right,
		    Real ***temp_left,
		    Real ***temp_down,
		    Real ***temp_up,
		    Real ***temp_right_send,
		    Real ***temp_left_send,
		    Real ***temp_down_send,
		    Real ***temp_up_send,
		    Real **temp_lu,
		    Real **temp_ld,
		    Real **temp_ru,
		    Real **temp_rd,
		    Real **temp_lu_send,
		    Real **temp_ld_send,
		    Real **temp_ru_send,
		    Real **temp_rd_send,
		    MPI_Comm mycomm,
		    MPI_Status *sta,
		    MPI_Request *req,
		    int count
		) {

	#ifdef SMALL_STEP
	STEPS = SMALL_STEP;
	#endif


	int i, j, k, l;
	int ***newFalgs = array3DI(x_sec+2,y_sec+2,Z);
	memset(&newFalgs[0][0][0], 0, sizeof(int)*(x_sec+2)*(y_sec+2)*Z);

	for(i = Xst; i < Xed; i++)
		for(j = Yst; j < Yed; j++)
			for(k = 0; k < Z; k++){
				if (flags[i-Xst+1][j-Yst+1][k] == FLUID)
					newFalgs[i-Xst+1][j-Yst+1][k] |= 1<<20;
				else if(flags[i-Xst+1][j-Yst+1][k] == BOUNCE)
					newFalgs[i-Xst+1][j-Yst+1][k] |= 1<<19;
				for(l = 0; l < 19; l++)
					if(walls[i-Xst][j-Yst][k][l])
						newFalgs[i-Xst+1][j-Yst+1][k] |= 1<<l;
			}

	struct athread_init_parameter parameter;
	parameter.nodes = nodes;
	parameter.flags = newFalgs;
	parameter.Xst = Xst;
	parameter.Xed = Xed;
	parameter.Yst = Yst;
	parameter.Yed = Yed;
	parameter.nz = Z;
	parameter.omega = omega;
	parameter.CSmago = CSmago;

	parameter.my2drank = my2drank;

	athread_init();
	
	// athread_spawn(slaveController, &parameter);

	athread_spawn(slaveInit, &parameter);
	athread_join();

	int s, n = 0;

	hch_timer_init_();

	for (s = 0; s < STEPS; s++) {

		// slaveStream(nodes, walls, newFalgs, Xst, Xed, Yst, Yed, Z, current, other);
		hch_timer_start(0);
		athread_spawn(SlaveCollide, &current);
		hch_timer_stop(0);

		hch_timer_start(1);
        bounce_send_init(X,
			Y,
			Z,
			Xst, 
			Xed, 
			Yst, 
			Yed, 
			x_sec, 
			y_sec, 
			current, 
			other, 
			nodes, 
			temp_left_send, 
			temp_right_send, 
			temp_up_send, 
			temp_down_send, 
			temp_ld_send, 
			temp_lu_send, 
			temp_rd_send, 
			temp_ru_send);
        hch_timer_stop(1);

		hch_timer_start(2);
        bounce_communicate(mycomm, 
			dims, 
			coords, 
			x_sec, 
			y_sec, 
			Z,
			&count,
			sta,
			req,
			temp_left_send, 
			temp_right_send, 
			temp_up_send, 
			temp_down_send, 
			temp_left, 
			temp_right, 
			temp_up, 
			temp_down,
			temp_lu_send, 
			temp_ld_send, 
			temp_ru_send, 
			temp_rd_send, 
			temp_lu, 
			temp_ld, 
			temp_ru, 
			temp_rd);
        hch_timer_stop(2);

		hch_timer_start(3);
		for(i = 0; i < count; i++) {
			MPI_Wait(&req[i], &sta[i]);
		}
        hch_timer_stop(3);

		hch_timer_start(4);
        bounce_update(X,
			Y,
			Z,
			Xst, 
			Xed, 
			Yst, 
			Yed, 
			myrank,
			x_sec, 
			y_sec,
			other,
			nodes,
			temp_left, 
			temp_right, 
			temp_up, 
			temp_down,
			temp_ld, 
			temp_lu, 
			temp_rd, 
			temp_ru);
        hch_timer_stop(4);

		hch_timer_start(5);
		masterStream(nodes, walls, flags, Xst, Xed, Yst, Yed, Z, current, other);
        hch_timer_stop(5);

		hch_timer_start(6);
		masterCollide(nodes, flags, Xst, Xed, Yst, Yed, Z, current);
        hch_timer_stop(6);

		hch_timer_start(7);
		athread_join();
		hch_timer_stop(7);

		other = current;
		current = (current+1)%2;

		if(myrank == 0 && STEPS >= 10 && (s + 1)%(STEPS/10) == 0.0) {
			n += 1;
			MLOG("Step >> [%d/%d] Calculation Completed %d%% \n", s + 1, STEPS, n * 10);
		}
	
	}
	hch_timer_finalize_();
}
