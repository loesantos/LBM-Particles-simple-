//================================================================================================//
//
//                       Functions to Simulate Particle Suspensions
//
//================================================================================================//

//	Lê os dados das partículas (número, densidade e diâmetro)

void read_data_particles ( int*, double*, double* );

//	Inicia as partículas (posição e velocidade)

void inicia_part ( int, double*, double*, double*, double*, double*, double*, int*, int, int, int );

//	Inicia arquivo que guarda variação de quant. de mov. do fluido

void inicia_delta ( double*, double*, double*, int*, int, int, int );

//	Realiza a interação entre fluido e partículas

void inter_fluid_part( double*, double*, double*, double*, double*, double*, double*, double*, 
					double*, int, double, double, double, double, int*, double*, double*, 
					int, int, int );

//	Realiza a interação entre fluido e partículas (com intepolação)

void inter_fluid_part_inter( double*, double*, double*, double*, double*, double*, double*, double*, 
					double*, double*, double*, 
					double*, int, double, double, double, double, int*, double*, double*, 
					int, int, int );

//	 Atualiza a posição das partículas (segunda ordem)

void update_part_second( double*, double*, double*, double*, double*, double*, double*, double*, 
						double*, double*, double*, 
						double*, int, double, int*, int, int, int );	

//	Grava a posição das partículas

void rec_particles ( double*, double*, double*, int, int );		

// 	Verifica a estabilização do escoamento

void verif_stab ( int, int&, int*, double*, double*, double*, double*, double*, double*, int, 
				int, int, int );	

//================================================================================================//



//===================== Lê os dados das partículas ==============================================//

void read_data_particles ( int* num_part, double* rho_part, double* diam_part )
{
	char nome_in[] = "data_particles.txt";

    ifstream f_in( nome_in );
    
    char st_num_part[25];
    
    f_in >> st_num_part;

    f_in >> *num_part;

    cout << "\nNúmero de partículas: " << *num_part << endl;
    
    char st_rho_part[25];
    
    f_in >> st_rho_part;

    f_in >> *rho_part;

    cout << "\nDensidade das partículas: " << *rho_part << endl;

	char st_diam_part[25];
    
    f_in >> st_diam_part;

    f_in >> *diam_part;

    cout << "\nDiâmetro das partículas: " << *diam_part << endl;

	f_in.close();

}

//===================== Inicia as partículas (posição e velocidade) ==============================//

void inicia_part ( int num_part, double* x_part, double* y_part, double* z_part, double* vx_part,
					double* vy_part, double* vz_part, int* ini_meio, int nx, int ny, int nz )

{
	for (int i = 0; i < num_part; i++ )
	{
		int *meio;
		
		do
		{
			x_part[i] = ( double ) ( rand() % nx );
			y_part[i] = ( double ) ( rand() % ny );
			z_part[i] = ( double ) ( rand() % nz );
			
			int x = ( int ) x_part[i];
			int y = ( int ) y_part[i];
			int z = ( int ) z_part[i];
			
			meio = ini_meio + x + y * nx + z * ny * nx;
		}
		//while ( *meio == 0  || x_part[i] < 2 || x_part[i] > nx - 2 
		//		|| y_part[i] < 2 || y_part[i] > ny - 2 );
			
		while ( *meio == 0  || x_part[i] < nx / 4 || x_part[i] > 3 * nx / 4 
				|| y_part[i] < ny / 4 || y_part[i] > 3 * ny / 4 );
			
		vx_part[i] = 0.;
		vy_part[i] = 0.;
		vz_part[i] = 0.;
	}
}

//===================== Inicia arquivo que guarda variação de quant. de mov. do fluido ===========//

void inicia_delta ( double* delta_Mx, double* delta_My, double* delta_Mz, int* ini_meio, int nx, 
					int ny, int nz )

{
	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{
				int *meio = ini_meio + x + y * nx + z * ny * nx;
				
				int coord = *meio - 1;
								
				delta_Mx[coord] = 0.0;
				delta_My[coord] = 0.0;
				delta_Mz[coord] = 0.0;
			}
		}
	}
}

//===================== Interação fluido-partícula ===============================================//

void inter_fluid_part( double* x_part, double* y_part, double* z_part, double* vx_part, 
						double* vy_part, double* vz_part, double* delta_Mx, double* delta_My, 
						double* delta_Mz, int num_part, double mass_part, double diam_part, 
						double area_part, double visc, int* ini_meio, double* ini_f, double* ini_c, 
						int nx, int ny, int nz )
{
	double *vx_viz = new double[8];
	double *vy_viz = new double[8];
	double *vz_viz = new double[8];	
	double *rho_viz = new double[8];

	double *weight = new double[8];

	for ( int i = 0; i < num_part; i++ )
	{
		double sum_vx = 0.0;
		double sum_vy = 0.0;
		double sum_vz = 0.0;
		double sum_rho = 0.0;
		
		double sum_weight = 0.0;
		
		int x_ant = ( int ) x_part[i];			
		int y_ant = ( int ) y_part[i];			
		int z_ant = ( int ) z_part[i];
		
		int viz = 0;
		
		for ( int x = x_ant; x <= x_ant + 1; x++ )
		{
			for ( int y = y_ant; y <= y_ant + 1; y++ )
			{
				for ( int z = z_ant; z <= z_ant + 1; z++ )
				{
					double d_viz = sqrt ( ( x - x_part[i] ) * ( x - x_part[i] ) 
										+ ( y - y_part[i] ) * ( y - y_part[i] )
										+ ( z - z_part[i] ) * ( z - z_part[i] ) );
										
					//weight[viz] = exp ( - d_viz * d_viz );
					
					weight[viz] = sqrt(3) - d_viz;
					
					//weight[viz] = 1. / ( d_viz + 1.e-10 );
															
					int *meio = ini_meio + ( x % nx ) + ( y % ny ) * nx + ( z % nz ) * ny * nx;

					if ( *meio )
					{
						int pto = *meio - 1;
								
						double *f = ini_f + ( pto ) * nvel;	
					
						double vx, vy, vz, rho;
					
						calcula ( f, ini_c, &vx, &vy, &vz, &rho );
						
						vx_viz[viz] = vx;		
						vy_viz[viz] = vy;		
						vz_viz[viz] = vz;
						
						rho_viz[viz] = rho;
					}
					else
					{
						vx_viz[viz] = 0.0;		
						vy_viz[viz] = 0.0;		
						vz_viz[viz] = 0.0;
						
						rho_viz[viz] = 1.0; 
						
						weight[viz] = 0.0;
					}
					
					sum_vx = sum_vx + weight[viz] * vx_viz[viz];
					sum_vy = sum_vy + weight[viz] * vy_viz[viz];
					sum_vz = sum_vz + weight[viz] * vz_viz[viz];
					
					sum_rho = sum_rho + weight[viz] * rho_viz[viz];
															
					sum_weight = sum_weight + weight[viz];
										
					viz++;				
				}
			}
		}
		
		if ( sum_weight == 0.0 ) sum_weight = 1.0;

		double vx_avg = sum_vx / sum_weight;
		double vy_avg = sum_vy / sum_weight;
		double vz_avg = sum_vz / sum_weight;
		
		double rho_avg = sum_rho / sum_weight;
		
		//------------- Força sobre a partícula --------------------------------------------------//
						
		double dvx = vx_avg - vx_part[i];
		double dvy = vy_avg - vy_part[i];
		double dvz = vz_avg - vz_part[i];
		
		double dv2 = dvx*dvx + dvy*dvy +  dvz*dvz;
		
		double mag_dv = sqrt( dv2 );
		
		double Re = diam_part * mag_dv / visc; // Reynolds number
		
		double fd = 1 + 0.15 * pow( Re, 0.687 ) ; // drag factor				
		
		double mag_fdrag = 9.424777961 * ( rho_avg * visc ) * diam_part * fd * mag_dv;
		
		//double cd = 1.; // drag coefficient	
		
		//if ( Re > 2.)	cd = ( 24 / Re ) * ( 1 + 0.15 * pow( Re, 0.687 ) );
		
		//else if ( Re > 1.e-4 )	cd = ( 24 / Re ) * ( 1 + Re * 3 / 16 ) ;
		
		//else cd = 240000.0;				
		
		//double mag_fdrag = 0.5 * rho_avg * cd * area_part * dv2; 
		
		double dir_x_fdrag = 0.0;
		double dir_y_fdrag = 0.0;
		double dir_z_fdrag = 0.0;
		
		if ( mag_dv > 0.0 )
		{
			dir_x_fdrag = dvx / mag_dv;
			dir_y_fdrag = dvy / mag_dv;
			dir_z_fdrag = dvz / mag_dv;
		}
		
		double f_drag_x = mag_fdrag * dir_x_fdrag;
		double f_drag_y = mag_fdrag * dir_y_fdrag;
		double f_drag_z = mag_fdrag * dir_z_fdrag;
		
		//-------------- Altera a velocidade das partículas --------------------------------------//

		double delta_vx_part = f_drag_x / mass_part;
		double delta_vy_part = f_drag_y / mass_part;
		double delta_vz_part = f_drag_z / mass_part;
		
		vx_part[i] = vx_part[i] + delta_vx_part;
		vy_part[i] = vy_part[i] + delta_vy_part;
		vz_part[i] = vz_part[i] + delta_vz_part;
		
		// Força as partículas para longe das paredes

		if ( y_part[i] < 1.0 && vy_part[i] < 0.0 ) vy_part[i] = - vy_part[i]/3.;
		
		if ( y_part[i] > ny-2 && vy_part[i] > 0.0 ) vy_part[i] = - vy_part[i]/3.;
		
		if ( x_part[i] < 1.0 && y_part[i] < ny-1 && vx_part[i] < 0.0 ) vx_part[i] = - vx_part[i]/3.;
		
		if ( x_part[i] > nx-2 && y_part[i] < ny-1 && vx_part[i] > 0.0) vx_part[i] = - vx_part[i]/3.;
		
		//-------------- Guarda a quantidade de movimento a ser alterada no fluido ---------------//
		
		viz = 0;
		
		for ( int x = x_ant; x <= x_ant + 1; x++ )
		{
			for ( int y = y_ant; y <= y_ant + 1; y++ )
			{
				for ( int z = z_ant; z <= z_ant + 1; z++ )
				{		
					int* meio = ini_meio + ( x % nx ) + ( y % ny ) * nx + ( z % nz ) * ny * nx;
					
					int coord = *meio - 1;

					double q_x = mass_part * delta_vx_part;
					double q_y = mass_part * delta_vy_part;
					double q_z = mass_part * delta_vz_part;

					delta_Mx[coord] = delta_Mx[coord] - q_x * weight[viz] / sum_weight;
					delta_My[coord] = delta_My[coord] - q_y * weight[viz] / sum_weight;
					delta_Mz[coord] = delta_Mz[coord] - q_z * weight[viz] / sum_weight;
					
					viz++;					
				}
			}
		}
	}
	
	delete[] vx_viz;
	delete[] vy_viz;
	delete[] vz_viz;
	delete[] rho_viz;
	
	delete[] weight;	
}

//===================== Interação fluido-partícula ===============================================//

void inter_fluid_part_inter( double* x_part, double* y_part, double* z_part, double* vx_part, 
						double* vy_part, double* vz_part, double* vx_part_old, 
						double* vy_part_old, double* vz_part_old,double* delta_Mx, double* delta_My, 
						double* delta_Mz, int num_part, double mass_part, double diam_part, 
						double area_part, double visc, int* ini_meio, double* ini_f, double* ini_c, 
						int nx, int ny, int nz )
{
	double *vx_viz = new double[8];
	double *vy_viz = new double[8];
	double *vz_viz = new double[8];	
	double *rho_viz = new double[8];

	double *weight = new double[8];

	for ( int i = 0; i < num_part; i++ )
	{		
		double sum_weight = 0.0;
		
		int x_0 = ( int ) x_part[i];			
		int y_0 = ( int ) y_part[i];			
		int z_0 = ( int ) z_part[i];
		
		int x_1 = x_0 + 1;
		int y_1 = y_0 + 1;
		int z_1 = z_0 + 1;
				
		for ( int x = x_0; x <= x_1; x++ )
		{
			for ( int y = y_0; y <= y_1; y++ )
			{
				for ( int z = z_0; z <= z_1; z++ )
				{
					int viz = ( x - x_0 ) + ( y - y_0 ) * 2 + ( z - z_0 ) * 4;
					
					double d_viz = sqrt ( ( x - x_part[i] ) * ( x - x_part[i] ) 
										+ ( y - y_part[i] ) * ( y - y_part[i] )
										+ ( z - z_part[i] ) * ( z - z_part[i] ) );
										
					//weight[viz] = exp ( - d_viz * d_viz );
					
					weight[viz] = sqrt(3) - d_viz;
					
					//weight[viz] = 1. / ( d_viz + 1.e-10 );
															
					int *meio = ini_meio + ( x % nx ) + ( y % ny ) * nx + ( z % nz ) * ny * nx;

					if ( *meio )
					{
						int pto = *meio - 1;
								
						double *f = ini_f + ( pto ) * nvel;	
					
						double vx, vy, vz, rho;
					
						calcula ( f, ini_c, &vx, &vy, &vz, &rho );
						
						vx_viz[viz] = vx;		
						vy_viz[viz] = vy;		
						vz_viz[viz] = vz;
						
						rho_viz[viz] = rho;
					}
					else
					{
						vx_viz[viz] = 0.0;		
						vy_viz[viz] = 0.0;		
						vz_viz[viz] = 0.0;
						
						rho_viz[viz] = 1.0; 
						
						weight[viz] = 0.0;
					}										
				}
			}
		}
		
		if ( sum_weight == 0.0 ) sum_weight = 1.0;
		
		//------------- Interpolação -------------------------------------------------------------//
		
		double x_d = ( x_part[i] - x_0 ) / ( x_1 - x_0 );
		double y_d = ( y_part[i] - y_0 ) / ( y_1 - y_0 );
		double z_d = ( z_part[i] - z_0 ) / ( z_1 - z_0 );
		
		//----------------------- vx ------------------------------//

		double cx00 = vx_viz[ 0 ] * ( 1. - x_d ) + vx_viz[ 1 ] * x_d;
		double cx01 = vx_viz[ 4 ] * ( 1. - x_d ) + vx_viz[ 5 ] * x_d;
		double cx10 = vx_viz[ 2 ] * ( 1. - x_d ) + vx_viz[ 3 ] * x_d;
		double cx11 = vx_viz[ 6 ] * ( 1. - x_d ) + vx_viz[ 7 ] * x_d;
		
		double cx0 = cx00 * ( 1. - y_d ) + cx10 * y_d;
		double cx1 = cx01 * ( 1. - y_d ) + cx11 * y_d;
		
		double vx_avg = cx0 * ( 1. - z_d) + cx1 * z_d;
		
		//----------------------- vy ------------------------------//

		double cy00 = vy_viz[ 0 ] * ( 1. - y_d ) + vy_viz[ 1 ] * y_d;
		double cy01 = vy_viz[ 4 ] * ( 1. - y_d ) + vy_viz[ 5 ] * y_d;
		double cy10 = vy_viz[ 2 ] * ( 1. - y_d ) + vy_viz[ 3 ] * y_d;
		double cy11 = vy_viz[ 6 ] * ( 1. - y_d ) + vy_viz[ 7 ] * y_d;
		
		double cy0 = cy00 * ( 1. - y_d ) + cy10 * y_d;
		double cy1 = cy01 * ( 1. - y_d ) + cy11 * y_d;
		
		double vy_avg = cy0 * ( 1. - z_d) + cy1 * z_d;
		
		//----------------------- vz ------------------------------//

		double cz00 = vz_viz[ 0 ] * ( 1. - z_d ) + vz_viz[ 1 ] * z_d;
		double cz01 = vz_viz[ 4 ] * ( 1. - z_d ) + vz_viz[ 5 ] * z_d;
		double cz10 = vz_viz[ 2 ] * ( 1. - z_d ) + vz_viz[ 3 ] * z_d;
		double cz11 = vz_viz[ 6 ] * ( 1. - z_d ) + vz_viz[ 7 ] * z_d;
		
		double cz0 = cz00 * ( 1. - y_d ) + cz10 * y_d;
		double cz1 = cz01 * ( 1. - y_d ) + cz11 * y_d;
		
		double vz_avg = cz0 * ( 1. - z_d) + cz1 * z_d;
		
		//----------------------- rho ------------------------------//

		double c00 = rho_viz[ 0 ] * ( 1. - x_d ) + rho_viz[ 1 ] * x_d;
		double c01 = rho_viz[ 4 ] * ( 1. - x_d ) + rho_viz[ 5 ] * x_d;
		double c10 = rho_viz[ 2 ] * ( 1. - x_d ) + rho_viz[ 3 ] * x_d;
		double c11 = rho_viz[ 6 ] * ( 1. - x_d ) + rho_viz[ 7 ] * x_d;
		
		double c0 = c00 * ( 1. - y_d ) + c10 * y_d;
		double c1 = c01 * ( 1. - y_d ) + c11 * y_d;
		
		double rho_avg = c0 * ( 1. - z_d) + c1 * z_d;
		
		//------------- Força sobre a partícula --------------------------------------------------//
						
		double dvx = vx_avg - vx_part[i];
		double dvy = vy_avg - vy_part[i];
		double dvz = vz_avg - vz_part[i];
		
		double dv2 = dvx*dvx + dvy*dvy +  dvz*dvz;
		
		double mag_dv = sqrt( dv2 );
		
		double Re = diam_part * mag_dv / visc; // Reynolds number
		
		double fd = 1 + 0.15 * pow( Re, 0.687 ) ; // drag factor				
		
		double mag_fdrag = 9.424777961 * ( rho_avg * visc ) * diam_part * fd * mag_dv;
		
		//double cd = 1.; // drag coefficient	
		
		//if ( Re > 2.)	cd = ( 24 / Re ) * ( 1 + 0.15 * pow( Re, 0.687 ) );
		
		//else if ( Re > 1.e-4 )	cd = ( 24 / Re ) * ( 1 + Re * 3 / 16 ) ;
		
		//else cd = 240000.0;				
		
		//double mag_fdrag = 0.5 * rho_avg * cd * area_part * dv2; 
		
		double dir_x_fdrag = 0.0;
		double dir_y_fdrag = 0.0;
		double dir_z_fdrag = 0.0;
		
		if ( mag_dv > 0.0 )
		{
			dir_x_fdrag = dvx / mag_dv;
			dir_y_fdrag = dvy / mag_dv;
			dir_z_fdrag = dvz / mag_dv;
		}
		
		double f_drag_x = mag_fdrag * dir_x_fdrag;
		double f_drag_y = mag_fdrag * dir_y_fdrag;
		double f_drag_z = mag_fdrag * dir_z_fdrag;
		
		//-------------- Altera a velocidade das partículas --------------------------------------//
		
		vx_part_old[i] = vx_part[i];
		vy_part_old[i] = vy_part[i];
		vz_part_old[i] = vz_part[i];

		double delta_vx_part = f_drag_x / mass_part;
		double delta_vy_part = f_drag_y / mass_part;
		double delta_vz_part = f_drag_z / mass_part;
		
		vx_part[i] = vx_part[i] + delta_vx_part;
		vy_part[i] = vy_part[i] + delta_vy_part;
		vz_part[i] = vz_part[i] + delta_vz_part;
		
		// Força as partículas para longe das paredes

		if ( y_part[i] < 1.0 && vy_part[i] < 0.0 ) vy_part[i] = - vy_part[i]/3.;
		
		if ( y_part[i] > ny-2 && vy_part[i] > 0.0 ) vy_part[i] = - vy_part[i]/3.;
		
		if ( x_part[i] < 1.0 && y_part[i] < ny-1 && vx_part[i] < 0.0 ) vx_part[i] = - vx_part[i]/3.;
		
		if ( x_part[i] > nx-2 && y_part[i] < ny-1 && vx_part[i] > 0.0) vx_part[i] = - vx_part[i]/3.;
		
		//-------------- Guarda a quantidade de movimento a ser alterada no fluido ---------------//
		
		for ( int x = x_0; x <= x_1; x++ )
		{
			for ( int y = y_0; y <= y_1; y++ )
			{
				for ( int z = z_0; z <= z_1; z++ )
				{		
					int viz = ( x - x_0 ) + ( y - y_0 ) * 2 + ( z - z_0 ) * 4;
					
					int* meio = ini_meio + ( x % nx ) + ( y % ny ) * nx + ( z % nz ) * ny * nx;
					
					int coord = *meio - 1;

					double q_x = mass_part * delta_vx_part;
					double q_y = mass_part * delta_vy_part;
					double q_z = mass_part * delta_vz_part;

					delta_Mx[coord] = delta_Mx[coord] - q_x * weight[viz] / sum_weight;
					delta_My[coord] = delta_My[coord] - q_y * weight[viz] / sum_weight;
					delta_Mz[coord] = delta_Mz[coord] - q_z * weight[viz] / sum_weight;		
				}
			}
		}
	}
	
	delete[] vx_viz;
	delete[] vy_viz;
	delete[] vz_viz;
	delete[] rho_viz;
	
	delete[] weight;	
}


//===================== Atualiza a posição das partículas ========================================//

void update_part_second( double* x_part, double* y_part, double* z_part, double* vx_part, 
						double* vy_part, double* vz_part, double* vx_part_old, double* vy_part_old, 
						double* vz_part_old, double* delta_Mx, double* delta_My, 
						double* delta_Mz, int num_part, double mass_part, int* ini_meio, 
						int nx, int ny, int nz )
{
	for ( int i = 0; i < num_part; i++ )
	{			
		double x_part_new = x_part[i] + vx_part_old[i] + 0.5 * delta_Mx[i] / mass_part;
		double y_part_new = y_part[i] + vy_part_old[i] + 0.5 * delta_My[i] / mass_part;
		double z_part_new = z_part[i] + vz_part_old[i] + 0.5 * delta_Mz[i] / mass_part;
		
		if ( x_part_new >= nx ) x_part_new = x_part_new - ( double ) nx;
		if ( y_part_new >= ny ) y_part_new = y_part_new - ( double ) ny;
		if ( z_part_new >= nz ) z_part_new = z_part_new - ( double ) nz;
		
		if ( x_part_new < 0. ) x_part_new = ( double ) nx - x_part_new;
		if ( y_part_new < 0. ) y_part_new = ( double ) ny - y_part_new;
		if ( z_part_new < 0. ) z_part_new = ( double ) nz - z_part_new;

		x_part[i] = x_part_new;
		y_part[i] = y_part_new;
		z_part[i] = z_part_new;			
	}
}

//===================== Grava a posição das partículas ===========================================//

void rec_particles ( double* x_part, double* y_part, double* z_part, int num_part, int passo )
{
	char nome_particles[50];

	sprintf ( nome_particles, "particles_%06d.csv", passo );

	ofstream file_particles ( nome_particles );
	
	file_particles << "X,Y,Z" << endl;

	for ( int i = 0; i < num_part; i++ )
	{
		file_particles << x_part[i]  << "," << y_part[i] << "," << z_part[i] << endl;
	}
	
	file_particles.close();
}

//==================== Verifica a estabilização do escoamento ====================================//

void verif_stab ( int passo, int& tmp_stab, int* ini_meio, double* ini_c, double* ini_N, 
				double* vx_avg, double* vy_avg, double* sum_vx, double* sum_vy, int ptos_meio, 
				int nx, int ny, int nz )
{
	#pragma omp parallel for

	for ( int x = 2; x < nx - 2; x++ )
	{
		for ( int y = 2; y < ny - 2; y++ )
		{
			for ( int z = 0; z < nz; z++ )
			{
				int *meio = ini_meio + x + y * nx + z * ny * nx;

				if ( *meio )
				{
					int pto = *meio - 1;
					
					double *N = ini_N + ( pto ) * nvel;
		
					double vx, vy, vz, rho;
					
					calcula ( N, ini_c, &vx, &vy, &vz, &rho );
					
					sum_vx[ pto ] = sum_vx[ pto ] + vx;
					sum_vy[ pto ] = sum_vy[ pto ] + vy;
				}
			}
		}
	}
	
	int tmp_avg = 100;
	
	if ( passo % tmp_avg == 0 && passo > tmp_avg )
	{
		double sum_err = 0.0;
		
		for ( int pto = 0; pto < ptos_meio; pto ++ )
		{
			double vx_avg_new = sum_vx[ pto ] / ( double ) tmp_avg;
			double vy_avg_new = sum_vy[ pto ] / ( double ) tmp_avg;
			
			sum_err = sum_err +( vx_avg_new - vx_avg[ pto ] ) * ( vx_avg_new - vx_avg[ pto ] ) +
						( vy_avg_new - vy_avg[ pto ] ) * ( vy_avg_new - vy_avg[ pto ] );
						
			vx_avg[ pto ] = vx_avg_new;
			vy_avg[ pto ] = vy_avg_new;
			
			sum_vx[ pto ] = 0.0;
			sum_vy[ pto ] = 0.0;
		}
		
		double erro = sum_err / (double) ptos_meio;
		
		if ( passo > 0 && passo % 1000 == 0 ) cout << "\n\nErro = " << erro << endl << endl;
		
		if ( erro < 1.0e-13  )
		{
			tmp_stab = passo;
			
			cout << "\n\nEstabilizou!!!" << endl << endl;
			
			getchar();
		}			
	}
}
//================================================================================================//




