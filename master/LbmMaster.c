#include "../header/LbmMaster.h"

volatile long host_flag[FLAG_SIZE];
volatile long slave_flag[FLAG_SIZE];
volatile long local_cc[FLAG_SIZE];
volatile long flag_to_wait;
unsigned long mpe_cc_cur_[100];
unsigned long mpe_cc_total_[100];
int my_rank, comm_sz;


void wait_slave_flag()
{
    while(slave_flag[0] < flag_to_wait);
    flag_to_wait++;
}

void terminate_athread_daemon()
{
    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = EXIT_FLAG;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
    wait_slave_flag();

    athread_join();
}

void athread_handshake()
{
    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = HANDSHAKE_FLAG;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
}


//ori data
int X, Y, Z, Xst, Xed, Yst, Yed, x_sec, y_sec;
int STEPS;
realt *****ori_nodes;
int ***ori_flags, ****ori_walls;
Real ***right_recv;
Real ***left_recv;
Real ***down_recv;
Real ***up_recv;
Real ***right_send;
Real ***left_send;
Real ***down_send;
Real ***up_send;
Real **lu_recv;
Real **ld_recv;
Real **ru_recv;
Real **rd_recv;
Real **lu_send;
Real **ld_send;
Real **ru_send;
Real **rd_send;

//new data

int halo_point_cnt;
int inner_point_cnt;

int *halo_point_xs, *inner_point_xs;
int *halo_point_ys, *inner_point_ys;

//insane
int *insane_halo_point_xs, *insane_inner_point_xs;
int *insane_halo_point_ys, *insane_inner_point_ys;

int insane_halo_point_cnt;
int insane_inner_point_cnt;

int *point_occu_flag;//x_sec * y_sec

//float**
Real* new_nodes;
int* new_flags, *new_walls;

struct lbm_init_param host_param;

#define DBG01 0

// analyse

extern void hch_timer_start(int num);
extern void hch_timer_stop(int num);
extern void hch_timer_manual(long* arr, int cnt, const char* fname, const char* header);

#define nodes_idx(d1, d2, d3, d4, d5) (((((d1) * (x_sec + 2) + (d2)) * (y_sec + 2) + (d3)) * Z + (d4)) * 19 + (d5))
#define flags_idx(d2, d3, d4) (((d2) * (y_sec + 2) + (d3)) * Z + (d4))
#define walls_idx(d2, d3, d4, d5) ((((d2) * y_sec + (d3)) * Z + (d4)) * 19 + (d5))

//faker

// Real* fk_std_other_line[3][3];//in
// Real* fk_std_current_line;
// int* fk_std_walls;
// int* fk_std_flags;

// Real fk_insane_other_line[4][4][(INSANE_SIZE + 2) * 19];//in
// Real fk_insane_current_line[2][2][INSANE_SIZE * 19];
// int fk_insane_walls[2][2][INSANE_SIZE * 19];
// int fk_insane_flags[2][2][INSANE_SIZE];

// Real* fk_insane_other_line[4][4];//in
// Real* fk_insane_current_line[2][2];
// int* fk_insane_walls[2][2];
// int* fk_insane_flags[2][2];



int cur_main_iter;

void fk_stream_and_collide_std_step(Real fk_std_other_line[3][3][(STEP_SIZE + 2) * 19], Real fk_std_current_line[STEP_SIZE * 19], int fk_std_walls[STEP_SIZE * 19], int fk_std_flags[STEP_SIZE])
{
    int inv, l;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew;

    int sis;
#define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
        if(fk_std_flags[sis] == FLUID) {
            for(l = 0; l < 19; l++) {
                inv = dfInv[l];
                fk_std_current_line[sis * 19 + l] = fk_std_other_line[e_x[inv] + 1][e_y[inv] + 1][other_step_index(sis + e_z[inv], l)];
            }
        }
        if(fk_std_flags[sis] == BOUNCE) {

            for(l = 0; l < 19; l++) {
                inv = dfInv[l];
                if(fk_std_walls[sis * 19 + l]) {
                    fk_std_current_line[sis * 19 + l] = fk_std_other_line[1][1][other_step_index(sis, inv)];
                } else
                {
                    fk_std_current_line[sis * 19 + l] = fk_std_other_line[e_x[inv] + 1][e_y[inv] + 1][other_step_index(sis + e_z[inv], l)];
                }
            }
        }

        if(fk_std_flags[sis] == FLUID || fk_std_flags[sis] == BOUNCE)
        {
            rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
            for (l = 0; l < 19; l++) {
                fi = fk_std_current_line[sis * 19 + l];
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
                        Sij[x1][y1] += e[k1][x1] * e[k1][y1] * (fk_std_current_line[sis * 19 + k1] - feq[k1]);
                    }
                    Qo += Sij[x1][y1] * Sij[x1][y1];
                }
            }

            nu = (2.0 / omega - 1.0) / 6.0;
            S = (-nu + sqrt(nu * nu + 18 * CSmago * CSmago * sqrt(Qo))) / 6.0 / CSmago / CSmago;
            omegaNew = 1.0 / (3.0 * (nu + CSmago * CSmago * S) + 0.5);




            for (l = 0; l < 19; l++) {
                fk_std_current_line[sis * 19 + l] =
                        (1.0 - omegaNew) * fk_std_current_line[sis * 19 + l] +
                        omegaNew * feq[l];
            }
        }
    }
}

// Real fk_step_current_line_read[STEP_SIZE * 19];//debug
int mpe_dbg_x, mpe_dbg_y;

void stream_and_collide_xy_std(int nz,
                           Real* current_line,
                           int* walls,
                           int* flags,
                           Real* other_line,
                           Real fk_std_other_line[3][3][(STEP_SIZE + 2) * 19], Real fk_std_current_line[STEP_SIZE * 19], int fk_std_walls[STEP_SIZE * 19], int fk_std_flags[STEP_SIZE])
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

        memcpy(fk_std_flags, flags + k, sizeof(int) * STEP_SIZE);
        memcpy(fk_std_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        memcpy(fk_std_current_line, current_line + k * 19, sizeof(Real) * STEP_SIZE * 19);

        // memcpy(fk_step_current_line_read, fk_std_current_line, sizeof(realt) * STEP_SIZE * 19);

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
                memcpy(fk_std_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19,
                        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19],
                        sizeof(Real) * (z_ed - z_st + 1) * 19);
            }
        }

        fk_stream_and_collide_std_step(fk_std_other_line, fk_std_current_line, fk_std_walls, fk_std_flags);

        memcpy(current_line + k * 19, fk_std_current_line, sizeof(Real) * STEP_SIZE * 19);
    }

}


void fk_stream_and_collide_insane_step(Real fk_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], Real fk_insane_current_line[2][2][INSANE_SIZE * 19], int fk_insane_walls[2][2][INSANE_SIZE * 19], int fk_insane_flags[2][2][INSANE_SIZE])
{
    int inv, l;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew;

    int local_x, local_y;

    int sis;
#define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(local_x = 0; local_x < 2; local_x++)
    {
        for(local_y = 0; local_y < 2; local_y++)
        {
            for(sis = 0; sis < INSANE_SIZE; sis++)
            {
                if(fk_insane_flags[local_x][local_y][sis] == FLUID) {
                    for(l = 0; l < 19; l++) {
                        inv = dfInv[l];
                        fk_insane_current_line[local_x][local_y][sis * 19 + l] = fk_insane_other_line[local_x + e_x[inv] + 1][local_y + e_y[inv] + 1][other_step_index(sis + e_z[inv], l)];
                    }
                }
                if(fk_insane_flags[local_x][local_y][sis] == BOUNCE) {

                    for(l = 0; l < 19; l++) {
                        inv = dfInv[l];
                        if(fk_insane_walls[local_x][local_y][sis * 19 + l]) {
                            fk_insane_current_line[local_x][local_y][sis * 19 + l] = fk_insane_other_line[local_x + 1][local_y + 1][other_step_index(sis, inv)];
                        } else
                        {
                            fk_insane_current_line[local_x][local_y][sis * 19 + l] = fk_insane_other_line[local_x + e_x[inv] + 1][local_y + e_y[inv] + 1][other_step_index(sis + e_z[inv], l)];
                        }
                    }
                }

                if(fk_insane_flags[local_x][local_y][sis] == FLUID || fk_insane_flags[local_x][local_y][sis] == BOUNCE)
                {
                    rho = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
                    for (l = 0; l < 19; l++) {
                        fi = fk_insane_current_line[local_x][local_y][sis * 19 + l];
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
                                Sij[x1][y1] += e[k1][x1] * e[k1][y1] * (fk_insane_current_line[local_x][local_y][sis * 19 + k1] - feq[k1]);
                            }
                            Qo += Sij[x1][y1] * Sij[x1][y1];
                        }
                    }

                    nu = (2.0 / omega - 1.0) / 6.0;
                    S = (-nu + sqrt(nu * nu + 18 * CSmago * CSmago * sqrt(Qo))) / 6.0 / CSmago / CSmago;
                    omegaNew = 1.0 / (3.0 * (nu + CSmago * CSmago * S) + 0.5);




                    for (l = 0; l < 19; l++) {
                        fk_insane_current_line[local_x][local_y][sis * 19 + l] =
                                (1.0 - omegaNew) * fk_insane_current_line[local_x][local_y][sis * 19 + l] +
                                omegaNew * feq[l];
                    }
                }
            }
        }
    }
    
}

void stream_and_collide_xy_insane(int nz,
                           Real* current_line,
                           int* walls,
                           int* flags,
                           Real* other_line,
                           Real fk_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], Real fk_insane_current_line[2][2][INSANE_SIZE * 19], int fk_insane_walls[2][2][INSANE_SIZE * 19], int fk_insane_flags[2][2][INSANE_SIZE])
{
    int k, l;

    int x_node_stride = (y_sec + 2) * Z * 19;
    int y_node_stride = Z * 19;
    int x_wall_stride = (y_sec) * Z * 19;
    int y_wall_stride = Z * 19;
    int x_flag_stride = (y_sec + 2) * Z;
    int y_flag_stride = Z;


    int total_step = nz / INSANE_SIZE;

    int cur_step;
    for(cur_step = 0; cur_step < total_step; cur_step++)
    {
        k = cur_step * INSANE_SIZE;
        int z;

        if(DBG01 && my_rank == 0)
            printf("insane step 1\n");

        int ex, ey;
        for(ex = 0; ex < 2; ex++)
        {
            for(ey = 0; ey < 2; ey++)
            {
                memcpy(fk_insane_flags[ex][ey], flags + (k + ex * x_flag_stride + ey * y_flag_stride), sizeof(int) * INSANE_SIZE);
                memcpy(fk_insane_walls[ex][ey], walls + (k * 19 + ex * x_wall_stride + ey * y_wall_stride), sizeof(int) * INSANE_SIZE * 19);
                memcpy(fk_insane_current_line[ex][ey], current_line + (k * 19 + ex * x_node_stride + ey * y_node_stride), sizeof(Real) * INSANE_SIZE * 19);
            }
        }
        if(DBG01 && my_rank == 0)
            printf("insane step 2\n");

        k = cur_step * INSANE_SIZE;
        int z_st = k - 1, z_ed = k + INSANE_SIZE;
        if(k == 0)
        {
            z_st = 0;
        }
        if(z_ed == nz)
        {
            z_ed--;
        }

        for(ex = -1; ex <= 2; ex++)
        {
            for(ey = -1; ey <= 2; ey++)
            {
                memcpy(fk_insane_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19,
                        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19],
                        sizeof(Real) * (z_ed - z_st + 1) * 19);
            }
        }

        if(DBG01 && my_rank == 0)
            printf("insane step 3\n");

        fk_stream_and_collide_insane_step(fk_insane_other_line, fk_insane_current_line, fk_insane_walls, fk_insane_flags);

        if(DBG01 && my_rank == 0)
            printf("insane step 4\n");

        for(ex = 0; ex < 2; ex++)
        {
            for(ey = 0; ey < 2; ey++)
            {
                memcpy(current_line + (k * 19 + ex * x_node_stride + ey * y_node_stride), fk_insane_current_line[ex][ey], sizeof(Real) * INSANE_SIZE * 19);
            }
        }

        if(DBG01 && my_rank == 0)
            printf("insane step 5\n");
    }

}

#define fk_nodes_idx(d2, d3, d4, d5) ((((d2) * (y_sec + 2) + (d3)) * Z + (d4)) * 19 + (d5))

void std_content(int my_core)
{
    Real fk_std_other_line[3][3][(STEP_SIZE + 2) * 19];//in
    Real fk_std_current_line[STEP_SIZE * 19];
    int fk_std_walls[STEP_SIZE * 19];
    int fk_std_flags[STEP_SIZE];
    //Real fk_std_other_line[3][3][(STEP_SIZE + 2) * 19], Real fk_std_current_line[STEP_SIZE * 19], int fk_std_walls[STEP_SIZE * 19], int fk_std_flags[STEP_SIZE]

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

        stream_and_collide_xy_std(host_param.Z, point_current, point_walls,
                              point_flags, point_other,
                              fk_std_other_line, fk_std_current_line, fk_std_walls, fk_std_flags);
    }
}

void insane_content(int my_core)
{
    Real fk_insane_other_line[4][4][(INSANE_SIZE + 2) * 19];//in
    Real fk_insane_current_line[2][2][INSANE_SIZE * 19];
    int fk_insane_walls[2][2][INSANE_SIZE * 19];
    int fk_insane_flags[2][2][INSANE_SIZE];
    // Real fk_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], Real fk_insane_current_line[2][2][INSANE_SIZE * 19], int fk_insane_walls[2][2][INSANE_SIZE * 19], int fk_insane_flags[2][2][INSANE_SIZE]

    int  insane_point_cnt = host_flag[INSANE_POINT_CNT];
    int* insane_xs_ptr = (int*)host_flag[INSANE_XS_PTR];
    int* insane_ys_ptr = (int*)host_flag[INSANE_YS_PTR];
    Real* current_head = (Real*)host_flag[STD_CURRENT_HEAD];
    Real* other_head = (Real*)host_flag[STD_OTHER_HEAD];

    int insane_my_start, insane_my_part, insane_max_part;
    insane_my_start = BLOCK_LOW(my_core, 64, insane_point_cnt);
    insane_my_part = BLOCK_SIZE(my_core, 64, insane_point_cnt);
    insane_max_part = BLOCK_SIZE(0, 64, insane_point_cnt);

    int* fk_insane_xs = insane_xs_ptr;
    int* fk_insane_ys = insane_ys_ptr;

    int i, j;
    for(i = insane_my_start; i < insane_my_start + insane_my_part; i++)
    {
        int x = fk_insane_xs[i];
        int y = fk_insane_ys[i];
        mpe_dbg_x = x;
        mpe_dbg_y = y;

        Real* point_current = current_head + fk_nodes_idx(x, y, 0, 0);
        Real* point_other = other_head + fk_nodes_idx(x, y, 0, 0);
        int* point_walls = host_param.walls + walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = host_param.flags + flags_idx(x, y, 0);

        if(DBG01 && my_rank == 0)
            printf("insane x = %d, y = %d\n", x, y);
        stream_and_collide_xy_insane(host_param.Z, point_current, point_walls,
                                     point_flags, point_other,
                                     fk_insane_other_line, fk_insane_current_line, fk_insane_walls, fk_insane_flags);
    }
}

void fake_core_func()
{
    int core;
    for(core = 0; core < 64; core++)
    {
        insane_content(core);
        std_content(core);
    }
}

void mix_inner_nowait(
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
    host_flag[INSANE_POINT_CNT] = insane_inner_point_cnt;
    host_flag[INSANE_XS_PTR] = (long)insane_inner_point_xs;
    host_flag[INSANE_YS_PTR] = (long)insane_inner_point_ys;

    //调用从核计算
    asm volatile ("nop":::"memory");
    host_flag[0] = host_flag[0] + 1;

    //主核emulator
    // fake_core_func();
}

void mix_halo(
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
    host_flag[STD_POINT_CNT] = halo_point_cnt;
    host_flag[STD_XS_PTR] = (long)halo_point_xs;
    host_flag[STD_YS_PTR] = (long)halo_point_ys;
    host_flag[STD_CURRENT_HEAD] = (long)(new_nodes + nodes_idx(current, 0, 0, 0, 0));
    host_flag[STD_OTHER_HEAD] = (long)(new_nodes + nodes_idx(other, 0, 0, 0, 0));
    host_flag[INSANE_POINT_CNT] = insane_halo_point_cnt;
    host_flag[INSANE_XS_PTR] = (long)insane_halo_point_xs;
    host_flag[INSANE_YS_PTR] = (long)insane_halo_point_ys;

    //调用从核计算
    asm volatile ("nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();

    //主核emulator
    // fake_core_func();
}

//事实上遍历的范围，xe--，ye--
void insane_halo_point_collect(int xs, int xe, int ys, int ye)
{
    xe--;
    ye--;
    int x, y;

    for(x = xs; x <= xe; x += 2)
    {
        for(y = ys; y <= ye; y += 2)
        {
            if(point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] == 0
             )
            {
                insane_halo_point_xs[insane_halo_point_cnt] = x;
                insane_halo_point_ys[insane_halo_point_cnt] = y;
                insane_halo_point_cnt++;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] = 1;
            }
        }
    }
}

void insane_inner_point_collect(int xs, int xe, int ys, int ye)
{
    xe--;
    ye--;
    int x, y;

    for(x = xs; x <= xe; x += 2)
    {
        for(y = ys; y <= ye; y += 2)
        {
            if(point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] == 0
             )
            {
                insane_inner_point_xs[insane_inner_point_cnt] = x;
                insane_inner_point_ys[insane_inner_point_cnt] = y;
                insane_inner_point_cnt++;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] = 1;
            }
        }
    }
}

void halo_point_collect(int xs, int xe, int ys, int ye)
{
    int x, y;
    for(x = xs; x <= xe; x++)
    {
        for(y = ys; y <= ye; y++)
        {
            if(point_occu_flag[(x - 1) * y_sec + (y - 1)] == 0)
            {
                halo_point_xs[halo_point_cnt] = x;
                halo_point_ys[halo_point_cnt] = y;
                halo_point_cnt++;
                point_occu_flag[(x - 1) * y_sec + (y - 1)] = 1;
            }
        }
    }
}

void inner_point_collect(int xs, int xe, int ys, int ye)
{
    int x, y;
    for(x = xs; x <= xe; x++)
    {
        for(y = ys; y <= ye; y++)
        {
            if(point_occu_flag[(x - 1) * y_sec + (y - 1)] == 0)
            {
                inner_point_xs[inner_point_cnt] = x;
                inner_point_ys[inner_point_cnt] = y;
                inner_point_cnt++;
                point_occu_flag[(x - 1) * y_sec + (y - 1)] = 1;
            }
        }
    }
}

void lbm_data_init(int _X, int _Y, int _Z, int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec,
                   realt *****_ori_nodes, int ***_ori_flags, int ****_ori_walls, int _STEPS,
                   Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                   Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                   Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                   Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send)
{
    //非常重要！！！
    //忘了加这一段，这才是各种玄学出错的根本原因

    int i;
    for(i = 0; i < FLAG_SIZE; i++)
    {
        host_flag[i] = 0;
        slave_flag[i] = 0;
    }
    flag_to_wait = 1;

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
    inner_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    inner_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));

    insane_inner_point_cnt = insane_halo_point_cnt = 0;
    insane_halo_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    insane_halo_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    insane_inner_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    insane_inner_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));

    point_occu_flag = malloc(sizeof(int) * (x_sec) * (y_sec));

    assert(halo_point_xs && halo_point_ys && inner_point_xs && inner_point_ys &&
           insane_halo_point_xs && insane_halo_point_ys && insane_inner_point_xs && insane_inner_point_ys && point_occu_flag);

    memset(point_occu_flag, 0, sizeof(int) * (x_sec) * (y_sec));

    // inner_point_collect(1, x_sec, 1, y_sec);

    insane_inner_point_collect(1 + 2, x_sec - 2, 1 + 2, y_sec - 2);
    inner_point_collect(1 + 2, x_sec - 2, 1 + 2, y_sec - 2);

    insane_halo_point_collect(1, x_sec, 1, 2);
    insane_halo_point_collect(1, 2, 1 + 2, y_sec - 2);
    insane_halo_point_collect(x_sec - 1, x_sec, 1 + 2, y_sec - 2);
    insane_halo_point_collect(1, x_sec, y_sec - 1, y_sec);
    halo_point_collect(1, x_sec, 1, 2);
    halo_point_collect(1, 2, 1 + 2, y_sec - 2);
    halo_point_collect(x_sec - 1, x_sec, 1 + 2, y_sec - 2);
    halo_point_collect(1, x_sec, y_sec - 1, y_sec);

    MPI_Barrier(MPI_COMM_WORLD);

    host_param.my_id = my_rank;
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

    //seq debug
    for(x = 0; x < comm_sz; x++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(x == my_rank)
        {
            printf("rank = %d, halo_insane = %d, halo_std = %d, inner_insane = %d, inner_std = %d\n",
                    my_rank, insane_halo_point_cnt * 4, halo_point_cnt, insane_inner_point_cnt * 4, inner_point_cnt);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

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

        mix_inner_nowait(current, other);

        if(DBG01 && my_rank == 0)
            assert(1 + 1 == 3);

        //send init
        hch_timer_start(MPF1_SEND_INIT);
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
        hch_timer_stop(MPF1_SEND_INIT);

        //comm
        hch_timer_start(MPF1_COMM);
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
        hch_timer_stop(MPF1_COMM);


        //waitall and update
        hch_timer_start(MPF1_WAIT_MPI);
        for(i = 0; i < wait_cnt; i++) {
            MPI_Wait(&req[i], &sta[i]);
        }
        hch_timer_stop(MPF1_WAIT_MPI);

        hch_timer_start(MPF1_UPDATE);
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
        hch_timer_stop(MPF1_UPDATE);

        hch_timer_start(MPF1_WAIT_INNER);
        wait_slave_flag();
        hch_timer_stop(MPF1_WAIT_INNER);

        //stream and collide
        hch_timer_start(MPF1_HALO);
        mix_halo(current, other);
        hch_timer_stop(MPF1_HALO);
        // stream_and_collide_inner(Z, current, other);


        other = current;
        current = (current+1)%2;
    }

    terminate_athread_daemon();

    hch_timer_manual((void*)&slave_flag[1], 11, "cpe_pf.csv", "my_rank, CPF1_WAIT, CPF1_STD_READ1, CPF1_STD_READ2, CPF1_STD_STREAM, CPF1_STD_COLLIDE, CPF1_STD_WRITE, CPF1_INSANE_READ1,CPF1_INSANE_READ2,CPF1_INSANE_STREA,CPF1_INSANE_COLLI,CPF1_INSANE_WRITE");

    // int x, y;
    // for(x = 0; x < 3; x++)
    // {
    //     for(y = 0; y < 3; y++)
    //     {
    //         free(fk_std_other_line[x][y]);
    //     }
    // }

    // free(fk_std_current_line);
    // free(fk_std_walls);
    // free(fk_std_flags);

    edcc = rpcc();
    double ttt = ccts(edcc - stcc);
    if(my_rank == 0)
    {
        printf("total time of main_iter is %f\n", ttt);
    }
}
