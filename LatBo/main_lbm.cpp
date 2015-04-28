// LatBo.cpp : Defines the entry point for the console application.

/*	
	**************************************************************************************************************
	**************************************************************************************************************
	**																											**
	**												LatBo MAIN													**
	**																											**
	**************************************************************************************************************
	**************************************************************************************************************
*/

#include "stdafx.h"
#include "LBM_definitions.h"	// Definitions file
#include "LBM_globalvars.h"			// Global variable references

using namespace std;	// Use the standard namespace

// Entry point
int _tmain( )
{
	/*
	***************************************************************************************************************
	********************************************* GENERAL INITIALISE **********************************************
	***************************************************************************************************************
	*/

	// Fix output format
	cout.precision(4);

	// Timing variables
	clock_t t_start, t_end; // Wall clock variables
	int t = 0;				// Time step counter, initially zero
	double tval = 0;		// Actual value of physical time (for multi-grid do not necessarily have unit time step)
	const int totalloops = (int)(T/deltat);		// Total number of loops to be performed (computed)

	

	/* ***************************************************************************************************************
	*********************************************** LEVEL 0 INITIALISE ***********************************************
	*************************************************************************************************************** */
	
	// Time step
	Grids[0].dt = deltat;

	// Store spacing
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	Grids[0].dx = 2*(Lx/(2*N));
	Grids[0].dy = 2*(Ly/(2*M));
	Grids[0].dz = 2*(Lz/(2*K));

#if (dims == 3)
	// Check that lattice volumes are cubes in 3D
	if ( (Lx/N) != (Ly/M) || (Lx/N) != (Lz/K) ) {
		cout << "Need to have lattice volumes which are cubes -- either change N/M/K or change domain dimensions" << endl;
		exit(EXIT_FAILURE);
	}
#else 
	// 2D so need square lattice cells
	if ( (Lx/N) != (Ly/M) ) {
		cout << "Need to have lattice cells which are squares -- either change N/M or change domain dimensions" << endl;
		exit(EXIT_FAILURE);
	}
#endif

	// Refined indices on L0
	Grids[0].XInd = onespace( Ref_startX, Ref_endX );
	Grids[0].YInd = onespace( Ref_startY, Ref_endY );
	Grids[0].ZInd = onespace( Ref_startZ, Ref_endZ );

	// L0 lattice site coordinates
	Grids[0].XPos = linspace( a_x + Grids[0].dx/2, b_x - Grids[0].dx/2, N );
	Grids[0].YPos = linspace( a_y + Grids[0].dy/2, b_y - Grids[0].dy/2, M );
	Grids[0].ZPos = linspace( a_z + Grids[0].dz/2, b_z - Grids[0].dz/2, K );


	// Initialise L0 macroscopic quantities
	// Velocity field
	Grids[0].u.resize( N*M*K*dims );
	LBM_init_vel(0);

	// Density field
	Grids[0].rho.resize( N*M*K );
	LBM_init_rho(0);

	// Initialise L0 matrices (f, feq) and typing matrix
	Grids[0].f.resize( N*M*K*nVels );
	Grids[0].feq.resize( N*M*K*nVels );
	Grids[0].LatTyp.resize( N*M*K );

	// Typing defined as follows:
	/*
	0 == boundary site
	1 == coarse site
	2 == fine/refined site
	3 == TL to upper (coarser) level
	4 == TL to lower (finer) level
	*/

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k < K; k++) {
				for (int v = 0; v < nVels; v++) {

					// Initialise f to feq
					int idx  = idxmap(i,j,k,v,M,K,nVels);
					Grids[0].f[idx] = LBM_collide(i,j,k,v,0);

				}

				// Label as coarse site
				int idx = idxmap(i,j,k,M,K);
				Grids[0].LatTyp[idx] = 1;

			}
		}
	}
	Grids[0].feq = Grids[0].f; // Make feq = feq too
		
	
	// Relaxation frequency on L0
	// Assign relaxation frequency corrected for grid and time step size
	Grids[0].omega = 1 / ( (nu / (Grids[0].dt*pow(cs,2)) ) + .5 );
	cout << "L0 relaxation time = " << (1/Grids[0].omega) << endl;

	/* ***************************************************************************************************************
	**************************************** REFINED LEVELS INITIALISE ***********************************************
	*************************************************************************************************************** */

	if (Nref != 0) {

		LBM_init_multi();
		
	}

	cout << "Initialisation Complete..." << endl;

	cout << "Initialising LBM time-stepping..." << endl;

	// Write out
	cout << "Writing Output to <Grids.out>..." << endl;
	for (int r = 0; r <= Nref; r++) {
		lbm_write3(r,t);
	}



	/* ***************************************************************************************************************
	********************************************** LBM PROCEDURE *****************************************************
	*************************************************************************************************************** */
	
	// LBM
	do {

		cout << "\n///////////////// Time Step " << t+1 << " /////////////////";

		// Start the clock
		t_start = clock();

		// Call LBM procedure from coarsest level (r = 0)
		LBM_multi(0);

		// Print Time of loop
		t_end = clock();
		clock_t secs = t_end - t_start;
		printf("\nLast time step took %f second(s)\n", ((double)secs)/CLOCKS_PER_SEC);

		// Increment counters
		t++;
		tval += Grids[0].dt;
       
	} while (tval < T);


	/* ***************************************************************************************************************
	*********************************************** POST PROCESS *****************************************************
	*************************************************************************************************************** */

	// None added yet


	return 0;
}

