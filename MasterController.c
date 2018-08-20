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

	struct athread_init_parameter parameter;
	parameter.nodes = nodes;
	parameter.flags = flags;
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

	int s, i, n = 0;

	for (s = 0; s < STEPS; s++) {

		slaveStream(nodes, walls, flags, Xst, Xed, Yst, Yed, Z, current, other);
		// slaveCollide(nodes, flags, Xst, Xed, Yst, Yed, Z, current);
		athread_spawn(SlaveCollide, &current);

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

		for(i = 0; i < count; i++) {
			MPI_Wait(&req[i], &sta[i]);
		}

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

		masterStream(nodes, walls, flags, Xst, Xed, Yst, Yed, Z, current, other);
		masterCollide(nodes, flags, Xst, Xed, Yst, Yed, Z, current);
		
		athread_join();

		other = current;
		current = (current+1)%2;

		if(myrank == 0 && STEPS >= 10 && (s + 1)%(STEPS/10) == 0.0) {
			n += 1;
			MLOG("Step >> [%d/%d] Calculation Completed %d%% \n", s + 1, STEPS, n * 10);
		}
	
	}
	// athread_join();
}
