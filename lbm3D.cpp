//====================== LB --> Lattice Boltzmann ================================================//

#include <omp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#define nvel 19
#define dim  3

double one_over_c_s2;

double W[nvel];

#include "functions_3d.cpp"

#include "particles.cpp"

//===================== Inicio do programa principal =============================================//

int main ()
{

    //================= Declaração de variáveis ==================================================//

    char nome_geo[50];  // Nome do arquivo de geometria

    double ftesc;       // Dimensão do pixel

    int npassos;        // Número de passos de simulação

    int numarq;         // Número de arquivos a serem gravados

    double tau;         // Tempo de relaxação

    double visc;        // Viscosidade

    double rho_ini;      // Densidade (partículas/sítio)

    //----------------- Le o arquivo de inicialização --------------------------------------------//

    read_data ( nome_geo, &ftesc, &npassos, &numarq, &tau, &visc, &rho_ini );

    //--------------------------------------------------------------------------------------------//

    int intervalo = npassos/numarq;

    int mais_passos = 0;

    int inicio = 0;

    int fim = npassos;

    double *temp;

    time_t t;

    srand((unsigned) time(&t));
    
    //time_t tempo;
    
    //int t_0 = ( int ) time( &tempo ); // tempo inicial em segundos

    //----------------- Específico para o modelo TRT ---------------------------------------------//

    //double tau_sim = tau;

    //double inv_tau = 1.0 / tau;

    //double tau_ant = ( ( 8.0 - inv_tau ) / ( 8.0 * (  2.0 - inv_tau ) ) );
        
    //================= Define o número de threads ===============================================//
    
    int max_threads = omp_get_max_threads();
	
	cout << "\nNumero max. de threads = " << max_threads << endl;
	
	int set_threads = max_threads / 2;
	
	cout << "\nThreads utilizadas = " << set_threads << endl;
	
	omp_set_num_threads( set_threads );

    //================= Le o arquivo de geometria do meio ========================================//

	int nx, ny, nz;

    ifstream fmatriz( nome_geo );

    string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ )	getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    fmatriz.close();
    
    int pts_in = 0;
    
    int pts_out = 0;
    
    nx = nx + pts_in + pts_out;

    int *ini_meio = new int[ nx * ny * nz ];

    int ptos_meio = read_geo ( nome_geo, ini_meio, pts_in, pts_out );
        
    double phi = ( double ) ptos_meio / ( double )( nx * ny * nz );

    cout << "\nPorosidade = " << phi << "    " << endl;
    
    //getchar();

    //=============== Grava imagem .vtk ==========================================================//

    char st_nome[15] = "meio_mdf.vtk";

    cout << "\nGravando geometria..." << endl;

    rec_geometria ( ini_meio, st_nome, nx, ny, nz );
    
    //============ Número, densidade e diâmetro das partículas ==========================================//

    int num_part; 		// Número de partículas
    
    double rho_part;   	// Densidade da fluido. = 1.0 (dens. do fluido =  0.001 g/cm^3) 
    
    double diam_part;  	// Diâmetro da partícula
    
    read_data_particles ( &num_part, &rho_part, &diam_part );
    
    double raio_part = 0.5 * diam_part;
    
    double area_part = 3.1416 * raio_part * raio_part;	// Área da seção    

    double mass_part = (4./3.) * 3.1415 * pow( raio_part, 3. ) * rho_part; // Massa da partícula
    
    //============ Aloca a memória/ inicializa (partículas) ======================================//
    
    double *x_part = new double[num_part];
    double *y_part = new double[num_part];
    double *z_part = new double[num_part];
    
    double *vx_part = new double[num_part];
    double *vy_part = new double[num_part];
    double *vz_part = new double[num_part];
    
    inicia_part ( num_part, x_part, y_part, z_part, vx_part, vy_part, vz_part, ini_meio, nx,ny,nz );
    
    double *vx_part_old = new double[num_part];
    double *vy_part_old = new double[num_part];
    double *vz_part_old = new double[num_part];
    
    inicia_part ( num_part, x_part, y_part, z_part, vx_part_old, vy_part_old, vz_part_old, ini_meio, 
					nx, ny, nz );
	
	double *delta_Mx = new double[ ptos_meio ];
	double *delta_My = new double[ ptos_meio ];
	double *delta_Mz = new double[ ptos_meio ];
	
	inicia_delta ( delta_Mx, delta_My, delta_Mz, ini_meio, nx, ny, nz );
	
	//============ Guarda as posições iniciais das partículas ====================================//
	
	double *x_0 = new double[num_part];
    double *y_0 = new double[num_part];
    
    for (int i = 0; i < num_part; i++ )
	{
		x_0[i] = x_part[i];
		y_0[i] = y_part[i];
	}

    //============ Aloca a memoria usada pelas outras matrizes ===================================//

    int tam_alloc = ptos_meio * nvel;

    double *ini_N = new double[ tam_alloc ];
        
    double *ini_N_novo = new double[ tam_alloc ];
    
    const int *ini_dir = new int[ tam_alloc ];
    
    //----------------------------------------------//
    
    double *sum_vx = new double[ ptos_meio ];
    double *sum_vy = new double[ ptos_meio ];   
    
    double *vx_avg = new double[ ptos_meio ];
    double *vy_avg = new double[ ptos_meio ];    
    
    for ( int pto = 0; pto < ptos_meio; pto ++ )
    {
		sum_vx[ pto ] = 0.0;
		sum_vy[ pto ] = 0.0;
		vx_avg[ pto ] = 0.0;
		vy_avg[ pto ] = 0.0;
	}

    //============ Define Os Vetores de rede c_i =================================================//

    double *ini_c = new double[nvel*dim];

    def_lattice_d3q19 ( ini_c, W, &one_over_c_s2 );
    
    double *ini_Q = new double[ dim * dim * nvel ];

    calcula_Q ( ini_c, ini_Q );

    //============ Define as direções de propagação ==============================================//

    def_dir_prop ( ini_c, ini_meio, ini_dir, nx, ny, nz );

    //================= Inicialização da rede ====================================================//

    double vx_ini = 0.0;

    double vy_ini = 0.0;

    double vz_ini = 0.0;
    
    int c_int = 110;
    
    char c;

    do
    {
        cout << "\nRecuperar dados? (s/n): ";
        c = getchar();
        cout << "\r                        ";
        c_int = ( int ) c;
    }
    
    while ( c_int != 115 && c != 110 );
	
    //----------------- Inicialização padrão -----------------------------------------------------//

    if ( c_int == 110 )
    {
        for ( int pto = 0; pto < ptos_meio; pto++ )
        {
			//cout << "end = " << ini_N + ( pto ) * nvel << endl;
			
            double* N = ini_N + ( pto ) * nvel;

            dist_eq ( vx_ini, vy_ini, vz_ini, rho_ini, ini_c, N, W, one_over_c_s2 );
                        
        }

        cout << endl;
    }

    //----------------- Inicialização com dados recuperados --------------------------------------//
	
    if ( c_int == 115 )
    {
        ifstream f_read ( "Arq_rec.dat" );

        f_read >> inicio;

        for ( int pto = 0; pto < ptos_meio; pto++ )
        {
            double vx, vy, vz, rho;

            f_read >> vx;

            f_read >> vy;

            f_read >> vz;

            f_read >> rho;

            double *N = ini_N + ( pto ) * nvel;

            dist_eq ( vx, vy, vz, rho, ini_c, N, W, one_over_c_s2 );
        }

        f_read.close();
 
        cout << endl;
    }
	
    //==================== Looping principal =====================================================//

	int tmp_stab = fim;
	
	loop:

    for ( int passo = inicio; passo < fim; passo++ )
    {
        if ( passo == 0 || passo % ( intervalo / 4 ) == 0 ) cout << "\rStep : " << passo;

        //------ Imposição de condições de contorno ----------------------------------------------//
        
        int y_top = ny - 1;
        
        double vx_imp = 0.1;

		bond_eq_rhovel_y ( ini_meio , ini_N, ini_c, y_top, rho_ini,	vx_imp, 0.0, 0.0, nx, ny, nz, 
							W, one_over_c_s2 );		
		
		if ( passo > tmp_stab ) //-------- Interação partícula-fluido ----------------------------//
		{				
	
			inter_fluid_part_inter( x_part, y_part, z_part, vx_part, vy_part, vz_part, 
						vx_part_old, vy_part_old, vz_part_old, delta_Mx, 
						delta_My, delta_Mz, num_part, mass_part, diam_part, area_part, visc, 
						ini_meio, ini_N, ini_c, nx, ny, nz );				

			//-------------- Altera a posição das partículas -------------------------------------//
			
			update_part_second(  x_part, y_part, z_part, vx_part, vy_part, vz_part, 
								vx_part_old, vy_part_old, vz_part_old, delta_Mx,  
								delta_My, delta_Mz, num_part, mass_part, ini_meio, nx, ny, nz );
		}	
		
		#pragma omp parallel for

       	for ( int x = 0; x < nx; x++ )
		{
			for ( int y = 0; y < ny; y++ )
			{
				for ( int z = 0; z < nz; z++ )
				{
					int *meio = ini_meio + x + y * nx + z * ny * nx;

					if ( *meio )
					{
						int pto = *meio - 1;

						double *N = ini_N + ( pto ) * nvel;

						//-------------- Etapa de colisão ----------------------------------------//
						
						double acc_x = 0.0;
						double acc_y = 0.0;
						double acc_z = 0.0;
						
						//--------- Inclui a interação partícula-fluido ---------//
												
						acc_x = acc_x + delta_Mx[pto];
						acc_y = acc_y + delta_My[pto];
						acc_z = acc_z + delta_Mz[pto];

						regularized_collision ( N, ini_c, ini_Q, tau, acc_x, acc_y, acc_z );
						
						delta_Mx[pto] = 0.0; 
						delta_My[pto] = 0.0;
						delta_Mz[pto] = 0.0;

						//------------- Etapa de propagacao --------------------------------------//

						propag_part_site ( ini_dir, ini_N, ini_N_novo, pto );

						//------------------------------------------------------------------------//
					}
				}
			}
		}
		
		//-------------- Calcula e grava o deslocamento médio ------------------------------------// 
		
		double tmp_caract = (double)nx / vx_imp;
		
		int step_dis = tmp_caract / 10;
		
		if ( passo % step_dis == 0 && passo > tmp_stab )
		{
			double mean_sqr_disp = 0.;
			
			for (int i = 0; i < num_part; i++ )
			{
				mean_sqr_disp = mean_sqr_disp + ( x_part[i] - x_0[i] )*( x_part[i] - x_0[i] ) 
				+ ( y_part[i] - y_0[i] )*( y_part[i] - y_0[i] );
			}
			
			mean_sqr_disp = mean_sqr_disp / (double) num_part;
			
			double tmp_rd = ( (double) passo - (double) tmp_stab ) / tmp_caract;
			
			double disp_rd = mean_sqr_disp / (double) ( nx * nx );
			
			ofstream file_MSD ( "AvgDisp.csv", ios::app );
    
			file_MSD << tmp_rd << "," << disp_rd << endl;
    
			file_MSD.close();
			
			if ( tmp_rd > 12. )
			{
				cout << "\n\nDeu o tempo!" << endl;
				
				getchar();
			}
		}

        //-------------- Grava a vorticidade & velocidade ----------------------------------------//
        
        if ( passo % intervalo == 0 && passo > 0 )
        {        
			cout << "\n\nGravando..." << endl;
			
			//rec_vorticity ( ini_meio, ini_N, ini_c, passo, nx, ny, nz );
			
			if ( passo % 1000 == 0 && passo > 0 ) 
				rec_velocity ( ini_meio, ini_N, ini_c, passo, nx, ny, nz );
			
			if ( passo > tmp_stab ) rec_particles ( x_part, y_part, z_part, num_part, passo );
			
			cout << "..." << endl << endl;
		}
		
		//------------------ Verifica a estabilização do escoamento ------------------------------//
		
		if ( passo < tmp_stab )
		{	
			verif_stab ( passo, tmp_stab, ini_meio, ini_c, ini_N, vx_avg, vy_avg, sum_vx, sum_vy, 
							ptos_meio,	nx, ny, nz );				
		}
		/*/
		if ( passo == tmp_stab )
		{
			ofstream file_vx_center ( "vx_center.csv" );
			
			int x_center = nx / 2;
			int z_center = nz / 2;

			for ( int y = 0; y < ny; y++ )
			{
				int *meio = ini_meio + x_center + y * nx + z_center * ny * nx;
				
				if ( *meio )
				{
					int pto = *meio - 1;
					
					double *N = ini_N + ( pto ) * nvel;
		
					double vx, vy, vz, rho;
					
					calcula ( N, ini_c, &vx, &vy, &vz, &rho );
					
					file_vx_center << y  << "," << vx << endl;
				}
				else file_vx_center << y  << "," << 0.0 << endl;
			}
			
			file_vx_center.close();
			
			cout << "\nGravou arquivo de velocidades." << endl;
			
			getchar();			
		}					
		/*/			
        //------------------ Atualiza ( fnovo => f ) ---------------------------------------------//

        temp = ini_N;

        ini_N = ini_N_novo;

        ini_N_novo = temp;

        //------------------- Grava arquivo de recuperação ---------------------------------------// 
    
        if ( passo == tmp_stab )
        {
            int precision = 8;

            rec_recovery( ini_N, ini_c, ptos_meio, passo, precision );
        }
        
    } //===================== Fim do looping principal ===========================================//

    cout << "\n\nMais quantos passos?" << endl << endl;
    cin >> mais_passos;

    if ( mais_passos != 0 )
    {
        inicio = fim;
        fim = inicio + mais_passos;
        goto loop;
    }

}

//================================================================================================//
