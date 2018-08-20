#include "../header/LbmMaster.h"

//new data

int halo_point_cnt;
int inner_point_cnt;

int *halo_point_xs, *inner_point_xs;
int *halo_point_ys, *inner_point_ys;
int *halo_point_ns, *inner_point_ns;
//float**
Real* new_nodes;
int* new_flags, *new_walls;

lbm_param host_param;

// analyse
extern void hch_timer_start(int num);
extern void hch_timer_stop(int num);

#define nodes_idx(d1, d2, d3, d4, d5) (((((d1) * (x_sec + 2) + (d2)) * (y_sec + 2) + (d3)) * Z + (d4)) * 19 + (d5))
#define flags_idx(d2, d3, d4) (((d2) * (y_sec + 2) + (d3)) * Z + (d4))
#define walls_idx(d2, d3, d4, d5) ((((d2) * y_sec + (d3)) * Z + (d4)) * 19 + (d5))


//faker

Real* fk_step_other_line[3][3];//in

Real* fk_step_current_line;
int* fk_step_walls;
int* fk_step_flags;

int cur_main_iter;

void fk_stream_and_collide_step()
{
    int inv, l;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew;

    int sis;
#define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
        if(fk_step_flags[sis] == FLUID) {
            for(l = 0; l < 19; l++) {
                inv = dfInv[l];
                fk_step_current_line[sis * 19 + l] = fk_step_other_line[e_x[inv] + 1][e_y[inv] + 1][other_step_index(sis + e_z[inv], l)];
            }
        }
        if(fk_step_flags[sis] == BOUNCE) {

            for(l = 0; l < 19; l++) {
                inv = dfInv[l];
                if(fk_step_walls[sis * 19 + l]) {
                    fk_step_current_line[sis * 19 + l] = fk_step_other_line[1][1][other_step_index(sis, inv)];
                } else
                {
                    fk_step_current_line[sis * 19 + l] = fk_step_other_line[e_x[inv] + 1][e_y[inv] + 1][other_step_index(sis + e_z[inv], l)];
                }
            }
        }

        if(fk_step_flags[sis] == FLUID || fk_step_flags[sis] == BOUNCE)
        {
            rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
            for (l = 0; l < 19; l++) {
                fi = fk_step_current_line[sis * 19 + l];
                rho += fi;
                u_x += e_x[l] * fi;
                u_y += e_y[l] * fi;
                u_z += e_z[l] * fi;
            }

            u_x /= rho;
            u_y /= rho;
            u_z /= rho;

            for (l = 0; l < 19; l++) {
                const Real tmp = (e_x[l] * u_x + e_y[l] * u_y + e_z[l] * u_z);
                feq[l] = w[l] * rho * (1.0 -
                                       (3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
                                       (3.0 *     tmp) +
                                       (9.0/2.0 * tmp * tmp));
            }

            Qo=0;
            Real Sij[3][3],S;
            Real e[19][3];
            int x1, y1, k1;
            for(k1 = 0; k1 < 19; k1++) {
                e[k1][0] = e_x[k1];
                e[k1][1] = e_y[k1];
                e[k1][2] = e_z[k1];
            }

            for(x1 = 0; x1 < 3; x1++) {
                for(y1 = 0; y1 < 3; y1++) {
                    Sij[x1][y1] = 0;
                    for(k1 = 0; k1 < 19; k1++) {
                        Sij[x1][y1] += e[k1][x1] * e[k1][y1] * (fk_step_current_line[sis * 19 + k1] - feq[k1]);
                    }
                    Qo += Sij[x1][y1] * Sij[x1][y1];
                }
            }

            nu = (2.0 / omega - 1.0) / 6.0;
            S = (-nu + sqrt(nu * nu + 18 * CSmago * CSmago * sqrt(Qo))) / 6.0 / CSmago / CSmago;
            omegaNew = 1.0 / (3.0 * (nu + CSmago * CSmago * S) + 0.5);




            for (l = 0; l < 19; l++) {
                fk_step_current_line[sis * 19 + l] =
                        (1.0 - omegaNew) * fk_step_current_line[sis * 19 + l] +
                        omegaNew * feq[l];
            }
        }
    }
}

Real fk_step_current_line_read[STEP_SIZE * 19];//debug
int mpe_dbg_x, mpe_dbg_y;

void stream_and_collide_xy(int nz,
                           Real* current_line,
                           int* walls,
                           int* flags,
                           Real* other_line)
{
    int k, l;

    int x_node_stride = (y_sec + 2) * Z * 19;
    int y_node_stride = Z * 19;

    int total_step = nz / STEP_SIZE;

    int cur_step;
    for(cur_step = 0; cur_step < total_step; cur_step++)
    {
        k = cur_step * STEP_SIZE;
        int z;

        memcpy(fk_step_flags, flags + k, sizeof(int) * STEP_SIZE);
        memcpy(fk_step_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        memcpy(fk_step_current_line, current_line + k * 19, sizeof(Real) * STEP_SIZE * 19);

        memcpy(fk_step_current_line_read, fk_step_current_line, sizeof(Real) * STEP_SIZE * 19);

        k = cur_step * STEP_SIZE;
        int z_st = k - 1, z_ed = k + STEP_SIZE;
        if(k == 0)
        {
            z_st = 0;
        }
        if(z_ed == nz)
        {
            z_ed--;
        }

        int ex, ey;
        for(ex = -1; ex <= 1; ex++)
        {
            for(ey = -1; ey <= 1; ey++)
            {
                memcpy(fk_step_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19,
                        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19],
                        sizeof(Real) * (z_ed - z_st + 1) * 19);
            }
        }

        fk_stream_and_collide_step();

        memcpy(current_line + k * 19, fk_step_current_line, sizeof(Real) * STEP_SIZE * 19);
    }

}

#define fk_nodes_idx(d2, d3, d4, d5) ((((d2) * (y_sec + 2) + (d3)) * Z + (d4)) * 19 + (d5))

void fake_core64_func(int my_core)
{
    int point_cnt = host_flag[STD_POINT_CNT];
    int* xs_ptr = (int*)host_flag[STD_XS_PTR];
    int* ys_ptr = (int*)host_flag[STD_YS_PTR];
    Real* current_head = (Real*)host_flag[STD_CURRENT_HEAD];
    Real* other_head = (Real*)host_flag[STD_OTHER_HEAD];

    int my_start, my_part, max_part;
    my_start = BLOCK_LOW(my_core, 64, point_cnt);
    my_part = BLOCK_SIZE(my_core, 64, point_cnt);
    max_part = BLOCK_SIZE(0, 64, point_cnt);


    int* fk_xs = xs_ptr;
    int* fk_ys = ys_ptr;

    int i, j;
    for(i = my_start; i < my_start + my_part; i++)
    {
        int x = fk_xs[i];
        int y = fk_ys[i];
        mpe_dbg_x = x;
        mpe_dbg_y = y;

        //var for this point
        Real* point_current = current_head + fk_nodes_idx(x, y, 0, 0);
        Real* point_other = other_head + fk_nodes_idx(x, y, 0, 0);
        int* point_walls = host_param.walls + walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = host_param.flags + flags_idx(x, y, 0);

        stream_and_collide_xy(host_param.Z, point_current, point_walls,
                              point_flags, point_other);
    }
}

void fake_core_func()
{
    int core;
    for(core = 0; core < 64; core++)
    {
        fake_core64_func(core);
    }
}

void stream_and_collide_inner(int Xst,
                              int Xed,
                              int Yst,
                              int Yed,
                              int nz,
                              int current,
                              int other)
{
    int i, j, k, l;
    int inv;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew;

    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = STD_LBM_FLAG;
    host_flag[STD_POINT_CNT] = inner_point_cnt;
    host_flag[STD_XS_PTR] = (long)inner_point_xs;
    host_flag[STD_YS_PTR] = (long)inner_point_ys;
    host_flag[STD_CURRENT_HEAD] = (long)(new_nodes + nodes_idx(current, 0, 0, 0, 0));
    host_flag[STD_OTHER_HEAD] = (long)(new_nodes + nodes_idx(other, 0, 0, 0, 0));

    //调用从核计算
    asm volatile ("nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
    athread_handshake();

    //主核emulator
    // fake_core_func();
}

void inner_point_collect(int xs, int xe, int ys, int ye)
{
    int x, y;
    for(x = xs; x <= xe; x++)
    {
        for(y = ys; y <= ye; y++)
        {
            inner_point_xs[inner_point_cnt] = x;
            inner_point_ys[inner_point_cnt] = y;
            inner_point_ns[inner_point_cnt] = x * (y_sec + 2) + y;
            inner_point_cnt++;
        }
    }
}

void lbm_data_init(int _X, int _Y, int _Z, int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec,
                   Real *****_ori_nodes, int ***_ori_flags, int ****_ori_walls, int _STEPS,
                   Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                   Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                   Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                   Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send)
{
    long stcc, edcc;
    stcc = rpcc();


    X = _X;
    Y = _Y;
    Z = _Z;
    Xst = _Xst;
    Xed = _Xed;
    Yst = _Yst;
    Yed = _Yed;
    x_sec = _x_sec;
    y_sec = _y_sec;

    ori_nodes = _ori_nodes;
    ori_walls = _ori_walls;
    ori_flags = _ori_flags;

    STEPS = _STEPS;

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

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    new_flags = &ori_flags[0][0][0];
    new_walls = &ori_walls[0][0][0][0];
    new_nodes = &ori_nodes[0][0][0][0][0];

    //超范围创建，为了方便
    inner_point_cnt = halo_point_cnt = 0;

    halo_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    halo_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    halo_point_ns = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));

    inner_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    inner_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    inner_point_ns = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));

    assert(halo_point_xs && halo_point_ys && halo_point_ns && inner_point_xs && inner_point_ys && inner_point_ns);

    inner_point_collect(1, x_sec, 1, y_sec);//initial version, inner is all

    MPI_Barrier(MPI_COMM_WORLD);

    host_param.master_id = my_rank;
    host_param.host_flag = (long*)&host_flag[0];
    host_param.slave_flag = (long*)&slave_flag[0];
    host_param.iter = STEPS;
    host_param.x_sec = x_sec;
    host_param.y_sec = y_sec;
    host_param.Z = Z;
    host_param.nu = nu;
    host_param.omega = omega;
    host_param.CSmago = CSmago;

    host_param.walls = new_walls;
    host_param.flags = new_flags;

    if(my_rank == 0)
    {
        printf("### MESSAGE FROM MPE ###\n");
        printf("iter = %d, x_sec = %d, y_sec = %d, Z = %d\n", host_param.iter, host_param.x_sec, host_param.y_sec, host_param.Z);
        printf("val: nu = %f, omega = %f, CSmago = %f\n", host_param.nu, host_param.omega, host_param.CSmago);
        printf("ptrs: %x, %x, %x, %x\n",
               host_param.host_flag, host_param.slave_flag, host_param.walls, host_param.flags);
        printf("sizeof(int) = %d, sizeof(Real) = %d\n", sizeof(int), sizeof(Real));
        printf("### MESSAGE FROM CPE ###\n");
    }

    //faker
    int x, y;
    for(x = 0; x < 3; x++)
    {
        for(y = 0; y < 3; y++)
        {
            fk_step_other_line[x][y] = malloc(sizeof(Real) * (STEP_SIZE + 2) * 19);
            assert(fk_step_other_line[x][y]);
        }
    }

    fk_step_current_line = malloc(sizeof(Real) * STEP_SIZE * 19);
    fk_step_walls = malloc(sizeof(int) * (STEP_SIZE * 19));
    fk_step_flags = malloc(sizeof(int) * (STEP_SIZE));

    assert(fk_step_current_line);
    assert(fk_step_walls);
    assert(fk_step_flags);

    edcc = rpcc();
    double ttt = ccts(edcc - stcc);
    if(my_rank == 0)
    {
        printf("total time of lbm_data_init is %f\n", ttt);
    }
}

void main_iter(int* _s)
{
    long stcc, edcc;
    stcc = rpcc();


    athread_init();
    athread_spawn(cpe_athread_daemon, (void*)&host_param);

    int current = 0, other = 1;
    int wait_cnt;
    int n = 0;//10 percent
    MPI_Status sta[16];
    MPI_Request req[16];

    int i;
    for(_s[0] = 0; _s[0] < STEPS; _s[0]++)
    {
        cur_main_iter = _s[0];
        if(my_rank == 0)
            printf("for(s = 0; s < STEPS; s++) s = %d\n", *_s);

        //send init
        hch_timer_start(0);
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
                         ori_nodes,
                         left_send,
                         right_send,
                         up_send,
                         down_send,
                         ld_send,
                         lu_send,
                         rd_send,
                         ru_send);
        hch_timer_stop(0);

        //comm
        hch_timer_start(1);
        bounce_communicate(MPI_COMM_WORLD,
                           NULL,
                           NULL,
                           x_sec,
                           y_sec,
                           Z,
                           &wait_cnt,
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
                           rd_recv );
        hch_timer_stop(1);


        //waitall and update
        hch_timer_start(2);
        for(i = 0; i < wait_cnt; i++) {
            MPI_Wait(&req[i], &sta[i]);
        }
        hch_timer_stop(2);

        hch_timer_start(3);
        bounce_update(X,
                      Y,
                      Z,
                      Xst,
                      Xed,
                      Yst,
                      Yed,
                      my_rank,
                      x_sec,
                      y_sec,
                      other,
                      ori_nodes,
                      left_recv,
                      right_recv,
                      up_recv,
                      down_recv,
                      ld_recv,
                      lu_recv,
                      rd_recv,
                      ru_recv);
        hch_timer_stop(3);

        //stream and collide
        hch_timer_start(4);
        stream_and_collide_inner(Xst, Xed, Yst, Yed, Z, current, other);
        // stream_and_collide_inner(Z, current, other);

        hch_timer_stop(4);

        //        if(*_s == 0)
        //        {
        //            int rank;
        //            for(rank = 0; rank < comm_sz; rank++)
        //            {
        //                MPI_Barrier(MPI_COMM_WORLD);
        //            }
        //        }


        other = current;
        current = (current+1)%2;
    }

    terminate_athread_daemon();

    int x, y;
    for(x = 0; x < 3; x++)
    {
        for(y = 0; y < 3; y++)
        {
            free(fk_step_other_line[x][y]);
        }
    }

    free(fk_step_current_line);
    free(fk_step_walls);
    free(fk_step_flags);

    edcc = rpcc();
    double ttt = ccts(edcc - stcc);
    if(my_rank == 0)
    {
        printf("total time of main_iter is %f\n", ttt);
    }
}


void terminate_athread_daemon()
{
    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = EXIT_FLAG;
    asm volatile ("nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
    wait_slave_flag();

    athread_join();
}

void wait_slave_flag()
{
    while(slave_flag[0] < flag_to_wait);
    flag_to_wait++;
}

void athread_handshake()
{
    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = HANDSHAKE_FLAG;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
}
