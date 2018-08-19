#include "../header/LbmMaster.h"

void data_init( int _myrank, int _comm_sz, MPI_Comm *_mycomm, 
		        int *_dims, int *_coords,
                int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec, 
                Real *****_nodes, int ***_flags, int ****_walls,
                Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send)
{
    myrank = _myrank;
    comm_sz = _comm_sz;

    mycomm = _mycomm;
    dims = _dims;
    coords = _coords;

    Xst = _Xst;
    Xed = _Xed;
    Yst = _Yst;
    Yed = _Yed;
    x_sec = _x_sec;
    y_sec = _y_sec;

    nodes = _nodes;
    walls = _walls;
    flags = _flags;

    right_recv = _temp_right;
    left_recv = _temp_left;
    down_recv = _temp_down;
    up_recv = _temp_up;
    right_send = _temp_right_send;
    left_send = _temp_left_send;
    down_send = _temp_down_send;
    up_send = _temp_up_send;

    lu_recv = _temp_lu;
    ld_recv = _temp_ld;
    ru_recv = _temp_ru;
    rd_recv = _temp_rd;
    lu_send = _temp_lu_send;
    ld_send = _temp_ld_send;
    ru_send = _temp_ru_send;
    rd_send = _temp_rd_send;

    param.myrank = _myrank;
    param.x_sec = _x_sec;
    param.y_sec = _y_sec;
    param.Z = Z;
    param.STEPS = STEPS;

    param.host_flag = (long*)&host_flag[0];
    param.slave_flag = (long*)&slave_flag[0];
}

void main_iter()
{
    //athread_init();
    //athread_spawn(cpe_athread_daemon, (void*)&param);

    int i, j, k, l, s;

    int current = 0, other = 1;
    MPI_Status sta[16];
    MPI_Request req[16];
    int n = 0;
    int count;

    for (s = 0; s < STEPS; s++) {

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
                         left_send, 
			 right_send, 
			 up_send, 
			 down_send, 
                         ld_send, 
			 lu_send, 
			 rd_send, 
			 ru_send);

        bounce_communicate(*mycomm, 
		           dims, 
			   coords, 
			   x_sec, 
			   y_sec, 
			   Z,
			   &count,
			   sta,
			   req,
			   left_send, 
			   right_send, 
			   up_send, 
			   down_send, 
			   left_recv, 
			   right_recv, 
			   up_recv, 
			   down_recv,
			   lu_send, 
			   ld_send, 
			   ru_send, 
			   rd_send, 
			   lu_recv, 
			   ld_recv, 
			   ru_recv, 
			   rd_recv);

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
                      left_recv, 
		      right_recv, 
		      up_recv, 
		      down_recv,
                      ld_recv, 
		      lu_recv, 
		      rd_recv, 
		      ru_recv);

		stream(nodes, walls, flags, Xst, Xed, Yst, Yed, Z, current, other);

		collide(nodes, flags, Xst, Xed, Yst, Yed, Z, current);

		other = current;
		current = (current+1)%2;

		if(myrank == 0 && STEPS >= 10 && (s + 1)%(STEPS/10) == 0.0) {
			n += 1;
			MLOG("Step >> [%d/%d] Calculation Completed %d%% \n", s + 1, STEPS, n * 10);
		}
	
	}

    //terminate_athread_daemon();
}
