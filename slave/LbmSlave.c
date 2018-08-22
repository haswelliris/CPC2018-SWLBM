#include "../header/SwDevice.h"

#define FLUID    0
#define BOUNCE   3

#define Real float

__thread_local Real s_w[19] __attribute__((aligned(32))) = {
        (1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
        (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };
__thread_local int s_e_x[19] __attribute__((aligned(32))) = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0};
__thread_local int s_e_y[19] __attribute__((aligned(32))) = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0};
__thread_local int s_e_z[19] __attribute__((aligned(32)))  = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0};
__thread_local Real s_e_T[19][4] __attribute__((aligned(32))) = {
	{0.,1.,0.,0}, {0.,-1.,0.,0.}, {1.,0.,0.,0.}, {-1.,0.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,-1.,0.},
	{1.,1.,0.,0.}, {-1.,1.,0.,0.}, {1.,-1.,0.,0.}, {-1.,-1.,0.,0.}, {0.,1.,1.,0.}, {0.,1.,-1.,0.},
	{0.,-1.,1.,0.}, {0.,-1.,-1.,0.}, {1.,0.,1.,0.}, {1.,0.,-1.,0.}, {-1.,0.,1.,0.}, {-1.,0.,-1.,0.}, {0.,0.,0.,0.}
};
__thread_local int s_dfInv[19] __attribute__((aligned(32))) = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};


void s_stream_insane_step(Real s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], Real s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE])
{
    int inv, l;

    int local_x, local_y;

    int sis;
#define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(local_x = 0; local_x < 2; local_x++)
    {
        for(local_y = 0; local_y < 2; local_y++)
        {
            for(sis = 0; sis < INSANE_SIZE; sis++)
            {
                if(s_insane_flags[local_x][local_y][sis] == FLUID) {
                    for(l = 0; l < 19; l++) {
                        inv = s_dfInv[l];
                        s_insane_current_line[local_x][local_y][sis * 19 + l] = s_insane_other_line[local_x + s_e_x[inv] + 1][local_y + s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                    }
                }
                if(s_insane_flags[local_x][local_y][sis] == BOUNCE) {

                    for(l = 0; l < 19; l++) {
                        inv = s_dfInv[l];
                        if(s_insane_walls[local_x][local_y][sis * 19 + l]) {
                            s_insane_current_line[local_x][local_y][sis * 19 + l] = s_insane_other_line[local_x + 1][local_y + 1][other_step_index(sis, inv)];
                        } else
                        {
                            s_insane_current_line[local_x][local_y][sis * 19 + l] = s_insane_other_line[local_x + s_e_x[inv] + 1][local_y + s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                        }
                    }
                }
            }
        }
    }
}

void s_collide_insane_step(Real s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], Real s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE])
{

    int local_x, local_y;

    int sis, l;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew, Sij[12];

#define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))

    for(local_x = 0; local_x < 2; local_x++)
    {
        for(local_y = 0; local_y < 2; local_y++)
        {
            for(sis = 0; sis < INSANE_SIZE; sis++)
            {
                if(s_insane_flags[local_x][local_y][sis] == FLUID || s_insane_flags[local_x][local_y][sis] == BOUNCE)
                {
                    floatv4 U = simd_set_floatv4(0., 0., 0., 0.), E;
                    #pragma unroll[19]
                    for (l = 0; l < 19; l++) {
                        // fi = a_current_line[l];
                        fi = s_insane_current_line[local_x][local_y][sis * 19 + l];
                        E = simd_set_floatv4(1., s_e_T[l][0], s_e_T[l][1], s_e_T[l][2]);
                        U = fi * E + U;
                    }
                    simd_store(U, &(Sij[0]));
                    rho = Sij[0], u_x = Sij[1], u_y = Sij[2], u_z = Sij[3];

                    u_x /= rho;
                    u_y /= rho;
                    u_z /= rho;

                    for (l = 0; l < 19; l++) {
                        const Real tmp = (s_e_x[l] * u_x + s_e_y[l] * u_y + s_e_z[l] * u_z);
                        feq[l] = s_w[l] * rho * (1.0 -
                                               (3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
                                               (3.0 *     tmp) +
                                               (9.0/2.0 * tmp * tmp));
                    }

                    Real S;
                    floatv4 vs1, vs2, vs3, vs4, vs5 ,ve, vk1, vk2, vk3, vk4, vk5;
                    int x1, y1, k1;
                    
                    vs1 = simd_set_floatv4(0., 0., 0., 0.);
                    vs2 = simd_set_floatv4(0., 0., 0., 0.);
                    vs3 = simd_set_floatv4(0., 0., 0., 0.);
                    #pragma unroll[19]
                    for(k1 = 0; k1 < 19; k1++) {
                        float tmp = (s_insane_current_line[local_x][local_y][sis * 19 + k1] - feq[k1]);
                        // float tmp = (a_current_line[k1] - feq[k1]);
                        simd_load(ve, &(s_e_T[k1][0]));
                        vk1 = tmp * s_e_T[k1][0];
                        vk2 = tmp * s_e_T[k1][1];
                        vk3 = tmp * s_e_T[k1][2];
                        
                        vs1 = vk1 * ve + vs1;
                        vs2 = vk2 * ve + vs2;
                        vs3 = vk3 * ve + vs3;
                    }
                    simd_store(vs1, &(Sij[0]));
                    simd_store(vs2, &(Sij[4]));
                    simd_store(vs3, &(Sij[8]));

                    //用simd乘加会丢精度，也不会变快
                    //优化流水线会变快，但是会因为计算顺序丢精度，需要调整汇编（但是我不会啊！！
                    Qo=0;
                    #pragma unroll[12]
                    for(x1 = 0; x1 < 12; ++x1)
                    {
                        Qo += Sij[x1] * Sij[x1];
                    }

                    device_param.nu = (2.0 / device_param.omega - 1.0) / 6.0;
                    S = (-device_param.nu + sqrt(device_param.nu * device_param.nu + 18 * device_param.CSmago * device_param.CSmago * sqrt(Qo))) / 6.0 / device_param.CSmago / device_param.CSmago;
                    omegaNew = 1.0 / (3.0 * (device_param.nu + device_param.CSmago * device_param.CSmago * S) + 0.5);

                    
                    
                    for (l = 0; l < 19; l++) {
                        s_insane_current_line[local_x][local_y][sis * 19 + l] =
                                (1.0 - omegaNew) * s_insane_current_line[local_x][local_y][sis * 19 + l] +
                                omegaNew * feq[l];
                    }
                }
            }
        }
    }
}

void s_stream_std(Real s_std_other_line[3][3][(STEP_SIZE + 2) * 19], Real s_std_current_line[STEP_SIZE * 19], int s_std_walls[STEP_SIZE * 19], int s_std_flags[STEP_SIZE])
{
    int inv, l;

    int sis;
    #define other_step_index(d1, d2) (((d1) + 1) * 19 + (d2))
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
        if(s_std_flags[sis] == FLUID) {
            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                s_std_current_line[sis * 19 + l] = s_std_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
            }
        }
        if(s_std_flags[sis] == BOUNCE) {

            for(l = 0; l < 19; l++) {
                inv = s_dfInv[l];
                if(s_std_walls[sis * 19 + l]) {
                    s_std_current_line[sis * 19 + l] = s_std_other_line[1][1][other_step_index(sis, inv)];
                } else 
                {
                    s_std_current_line[sis * 19 + l] = s_std_other_line[s_e_x[inv] + 1][s_e_y[inv] + 1][other_step_index(sis + s_e_z[inv], l)];
                }
            }
        }
    }
}

void s_collide_std(Real s_std_other_line[3][3][(STEP_SIZE + 2) * 19], Real s_std_current_line[STEP_SIZE * 19], int s_std_walls[STEP_SIZE * 19], int s_std_flags[STEP_SIZE])
{
    int l;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew, Sij[12];

    int sis;
    for(sis = 0; sis < STEP_SIZE; sis++)
    {
        if(s_std_flags[sis] == FLUID || s_std_flags[sis] == BOUNCE)
        {
            floatv4 U = simd_set_floatv4(0., 0., 0., 0.), E;
            #pragma unroll[19]
            for (l = 0; l < 19; l++) {
                // fi = a_current_line[l];
                fi = s_std_current_line[sis * 19 + l];
                E = simd_set_floatv4(1., s_e_T[l][0], s_e_T[l][1], s_e_T[l][2]);
                U = fi * E + U;
            }
            simd_store(U, &(Sij[0]));
            rho = Sij[0], u_x = Sij[1], u_y = Sij[2], u_z = Sij[3];

            u_x /= rho;
            u_y /= rho;
            u_z /= rho;

            for (l = 0; l < 19; l++) {
                const Real tmp = (s_e_x[l] * u_x + s_e_y[l] * u_y + s_e_z[l] * u_z);
                feq[l] = s_w[l] * rho * (1.0 -
                                       (3.0/2.0 * (u_x * u_x + u_y * u_y + u_z * u_z)) +
                                       (3.0 *     tmp) +
                                       (9.0/2.0 * tmp * tmp));
            }

            Real S;
			floatv4 vs1, vs2, vs3, vs4, vs5 ,ve, vk1, vk2, vk3, vk4, vk5;
            int x1, y1, k1;
			
			vs1 = simd_set_floatv4(0., 0., 0., 0.);
			vs2 = simd_set_floatv4(0., 0., 0., 0.);
			vs3 = simd_set_floatv4(0., 0., 0., 0.);
            #pragma unroll[19]
			for(k1 = 0; k1 < 19; k1++) {
				float tmp = (s_std_current_line[sis * 19 + k1] - feq[k1]);
				// float tmp = (a_current_line[k1] - feq[k1]);
				simd_load(ve, &(s_e_T[k1][0]));
				vk1 = tmp * s_e_T[k1][0];
				vk2 = tmp * s_e_T[k1][1];
				vk3 = tmp * s_e_T[k1][2];
				
				vs1 = vk1 * ve + vs1;
				vs2 = vk2 * ve + vs2;
				vs3 = vk3 * ve + vs3;
            }
			simd_store(vs1, &(Sij[0]));
			simd_store(vs2, &(Sij[4]));
			simd_store(vs3, &(Sij[8]));

            //用simd乘加会丢精度，也不会变快
            Qo=0;
            #pragma unroll[12]
			for(x1 = 0; x1 < 12; ++x1)
			{
				Qo += Sij[x1] * Sij[x1];
			}


            device_param.nu = (2.0 / device_param.omega - 1.0) / 6.0;
            S = (-device_param.nu + sqrt(device_param.nu * device_param.nu + 18 * device_param.CSmago * device_param.CSmago * sqrt(Qo))) / 6.0 / device_param.CSmago / device_param.CSmago;
            omegaNew = 1.0 / (3.0 * (device_param.nu + device_param.CSmago * device_param.CSmago * S) + 0.5);


            for (l = 0; l < 19; l++) {
               s_std_current_line[sis * 19 + l] =
                       (1.0 - omegaNew) * s_std_current_line[sis * 19 + l] +
                       omegaNew * feq[l];
            }
        }
    }
}

void stream_and_collide_xy_insane(int nz,
                           Real* current_line,
                           int* walls,
                           int* flags,
                           Real* other_line,
                           Real s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], Real s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE])
{
    int k, l;

    int x_node_stride = (device_param.y_sec + 2) * device_param.Z * 19;
    int y_node_stride = device_param.Z * 19;
    int x_wall_stride = (device_param.y_sec) * device_param.Z * 19;
    int y_wall_stride = device_param.Z * 19;
    int x_flag_stride = (device_param.y_sec + 2) * device_param.Z;
    int y_flag_stride = device_param.Z;


    int total_step = nz / INSANE_SIZE;

    int cur_step;
    for(cur_step = 0; cur_step < total_step; cur_step++)
    {
        k = cur_step * INSANE_SIZE;
        int z;

        int ex, ey;
        for(ex = 0; ex < 2; ex++)
        {
            for(ey = 0; ey < 2; ey++)
            {
                async_get(s_insane_flags[ex][ey], flags + (k + ex * x_flag_stride + ey * y_flag_stride), sizeof(int) * INSANE_SIZE);
                async_get(s_insane_walls[ex][ey], walls + (k * 19 + ex * x_wall_stride + ey * y_wall_stride), sizeof(int) * INSANE_SIZE * 19);
                async_get(s_insane_current_line[ex][ey], current_line + (k * 19 + ex * x_node_stride + ey * y_node_stride), sizeof(Real) * INSANE_SIZE * 19);
            }
        }

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
                async_get(s_insane_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19,
                        &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19],
                        sizeof(Real) * (z_ed - z_st + 1) * 19);
            }
        }
        wait_all_async_get();

        s_stream_insane_step(s_insane_other_line, s_insane_current_line, s_insane_walls, s_insane_flags);
        s_collide_insane_step(s_insane_other_line, s_insane_current_line, s_insane_walls, s_insane_flags);

        for(ex = 0; ex < 2; ex++)
        {
            for(ey = 0; ey < 2; ey++)
            {
                async_put(current_line + (k * 19 + ex * x_node_stride + ey * y_node_stride), s_insane_current_line[ex][ey], sizeof(Real) * INSANE_SIZE * 19);
            }
        }
        wait_all_async_put();

    }

}

__thread_local int dbg_print_cnt;

__thread_local int dbg_x, dbg_y;

void s_stream_and_collide_xy(int nz,
                           Real* current_line,
                           int* walls,
                           int* flags,
                           Real* other_line,
                           Real s_std_other_line[3][3][(STEP_SIZE + 2) * 19], Real s_std_current_line[STEP_SIZE * 19], int s_std_walls[STEP_SIZE * 19], int s_std_flags[STEP_SIZE])
{
    int k, l;

    int x_node_stride = (device_param.y_sec + 2) * device_param.Z * 19;
    int y_node_stride = device_param.Z * 19;

    int total_step = device_param.Z / STEP_SIZE;


    int cur_step;
    for(cur_step = 0; cur_step < total_step; cur_step++)
    {
        k = cur_step * STEP_SIZE;
        int z;

        async_get(s_std_flags, flags + k, sizeof(int) * STEP_SIZE);
        async_get(s_std_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        async_get(s_std_current_line, current_line + k * 19, sizeof(Real) * STEP_SIZE * 19);

        // memcpy(s_std_flags, flags + k, sizeof(int) * STEP_SIZE);
        // memcpy(s_std_walls, walls + k * 19, sizeof(int) * STEP_SIZE * 19);
        // memcpy(s_std_current_line, current_line + k * 19, sizeof(int) * STEP_SIZE * 19);

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
                async_get(s_std_other_line[ex + 1][ey + 1] + (z_st - k + 1) * 19, 
                       &other_line[ex * x_node_stride + ey * y_node_stride + z_st * 19], 
                       sizeof(Real) * (z_ed - z_st + 1) * 19);
            }
        }
        wait_all_async_get();

        s_stream_std(s_std_other_line, s_std_current_line, s_std_walls, s_std_flags);
        s_collide_std(s_std_other_line, s_std_current_line, s_std_walls, s_std_flags);

        sync_put(current_line + k * 19, s_std_current_line, sizeof(Real) * STEP_SIZE * 19);
    }

}

#define s_nodes_idx(d2, d3, d4, d5) ((((d2) * (device_param.y_sec + 2) + (d3)) * device_param.Z + (d4)) * 19 + (d5))
#define s_flags_idx(d2, d3, d4) (((d2) * (device_param.y_sec + 2) + (d3)) * device_param.Z + (d4))
#define s_walls_idx(d2, d3, d4, d5) ((((d2) * device_param.y_sec + (d3)) * device_param.Z + (d4)) * 19 + (d5))

void insane_lbm()
{
    Real s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19];//in
    Real s_insane_current_line[2][2][INSANE_SIZE * 19];
    int s_insane_walls[2][2][INSANE_SIZE * 19];
    int s_insane_flags[2][2][INSANE_SIZE];
    // Real s_insane_other_line[4][4][(INSANE_SIZE + 2) * 19], Real s_insane_current_line[2][2][INSANE_SIZE * 19], int s_insane_walls[2][2][INSANE_SIZE * 19], int s_insane_flags[2][2][INSANE_SIZE]

    int  insane_point_cnt = device_in_param[INSANE_POINT_CNT];
    int* insane_xs_ptr = (int*)device_in_param[INSANE_XS_PTR];
    int* insane_ys_ptr = (int*)device_in_param[INSANE_YS_PTR];
    Real* current_head = (Real*)device_in_param[STD_CURRENT_HEAD];
    Real* other_head = (Real*)device_in_param[STD_OTHER_HEAD];

    int insane_my_start, insane_my_part, insane_max_part;
    insane_my_start = BLOCK_HEAD(device_core_id, 64, insane_point_cnt);
    insane_my_part = BLOCK_SIZE(device_core_id, 64, insane_point_cnt);
    insane_max_part = BLOCK_SIZE(0, 64, insane_point_cnt);

    int* s_insane_xs = insane_xs_ptr;
    int* s_insane_ys = insane_ys_ptr;

    int i, j;
    for(i = insane_my_start; i < insane_my_start + insane_my_part; i++)
    {
        int x = s_insane_xs[i];//GLD!!!
        int y = s_insane_ys[i];//GLD!!!

        Real* point_current = current_head + s_nodes_idx(x, y, 0, 0);
        Real* point_other = other_head + s_nodes_idx(x, y, 0, 0);
        int* point_walls = device_param.walls + s_walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = device_param.flags + s_flags_idx(x, y, 0);

        stream_and_collide_xy_insane(device_param.Z, point_current, point_walls,
                                     point_flags, point_other,
                                     s_insane_other_line, s_insane_current_line, s_insane_walls, s_insane_flags);
    }
}

void std_lbm()
{
    Real s_std_other_line[3][3][(STEP_SIZE + 2) * 19];//in
    Real s_std_current_line[STEP_SIZE * 19];
    int s_std_walls[STEP_SIZE * 19];
    int s_std_flags[STEP_SIZE];

    dbg_print_cnt = 0;
    int point_cnt = device_in_param[STD_POINT_CNT];
    int* xs_ptr = (int*)device_in_param[STD_XS_PTR];
    int* ys_ptr = (int*)device_in_param[STD_YS_PTR];
    Real* current_head = (Real*)device_in_param[STD_CURRENT_HEAD];
    Real* other_head = (Real*)device_in_param[STD_OTHER_HEAD];

    int my_start, my_part, max_part;
    my_start = BLOCK_HEAD(device_core_id, 64, point_cnt);
    my_part = BLOCK_SIZE(device_core_id, 64, point_cnt);
    max_part = BLOCK_SIZE(0, 64, point_cnt);

    int* s_xs = xs_ptr;
    int* s_ys = ys_ptr;

    int i, j;
    for(i = my_start; i < my_start + my_part; i++)
    {
        int x = s_xs[i];//GLD!!!
        int y = s_ys[i];//GLD!!!
        dbg_x = x;
        dbg_y = y;

        //var for this point
        Real* point_current = current_head + s_nodes_idx(x, y, 0, 0);
        Real* point_other = other_head + s_nodes_idx(x, y, 0, 0);
        int* point_walls = device_param.walls + s_walls_idx(x - 1, y - 1, 0, 0);
        int* point_flags = device_param.flags + s_flags_idx(x, y, 0);

        s_stream_and_collide_xy(device_param.Z, point_current, point_walls, 
                              point_flags, point_other,
                              s_std_other_line, s_std_current_line, s_std_walls, s_std_flags);
    }
}

void device_run()
{
    std_lbm();
    insane_lbm();
}


