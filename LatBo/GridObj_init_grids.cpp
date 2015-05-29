/* This file contains the routines necessary to initialise the macroscopic quantities.
*/

#include "stdafx.h"
#include "GridObj.h"
#include "definitions.h"
#include "globalvars.h"
#include <math.h>

using namespace std;

// ***************************************************************************************************

// Initialise velocity method
void GridObj::LBM_init_vel ( ) {

	// Max velocity
	double u_in[3] = {u_0x, u_0y, u_0z};
	
	// Get grid sizes
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();

	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				int idx = idxmap(i,j,k,M_lim,K_lim);
				int idx_x = idxmap(i,j,k,0,M_lim,K_lim,dims);
				int idx_y = idxmap(i,j,k,1,M_lim,K_lim,dims);
				int idx_z = idxmap(i,j,k,2,M_lim,K_lim,dims);
				
				// Uniform flow
				u[idx_x] = u_in[0];
				u[idx_y] = u_in[1];
#if (dims == 3)
				u[idx_z] = u_in[2];
#endif

			}
		}
	}

#ifdef SOLID_ON
	// Perform solid site reset of velocity
	solidSiteReset();
#endif

}


// ***************************************************************************************************

// Initialise density method
void GridObj::LBM_init_rho ( ) {

	// Get grid sizes
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();

	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {		
			for (int k = 0; k < K_lim; k++) {

				// Max velocity
				double u_in[3] = {u_0x, u_0y, u_0z};

				int idx = idxmap(i,j,k,M_lim,K_lim);

				// Uniform Density
				rho[idx] = rho_in;
			}
		}
	}

}

// ***************************************************************************************************

// Initialise wall labels method (L0 only)
void GridObj::LBM_init_wall_lab ( ) {

	int i, idx;

#ifdef SOLID_ON
	for (int i = obj_x_min; i <= obj_x_max; i++) {
		for (int j = obj_y_min; j <= obj_y_max; j++) {
			for (int k = obj_z_min; k <= obj_z_max; k++) {

				idx = idxmap(i,j,k,M,K);
				LatTyp[idx] = 0;

			}
		}
	}
#endif

#ifdef INLET_ON
	i = 0;
	for (int j = 0; j < M; j++) {
		for (int k = 0; k < K; k++) {

			idx = idxmap(i,j,k,M,K);
			LatTyp[idx] = 7;

		}
	}
#endif

#ifdef OUTLET_ON
	i = N-1;
	for (int j = 0; j < M; j++) {
		for (int k = 0; k < K; k++) {

			idx = idxmap(i,j,k,M,K);
			LatTyp[idx] = 8;

		}
	}
#endif

}


// ***************************************************************************************************

// Initialise level 0 grid method
void GridObj::LBM_init_grid( ) {

	// Store spacing
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	dx = 2*(Lx/(2*N));
	dy = 2*(Ly/(2*M));
	dz = 2*(Lz/(2*K));

	// Time step = grid spacing
	dt = dx;

	
	// Checks to make sure grid cell dimensions are suitable
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


	// Checks passed so continue...
	// NODE NUMBERS on L0
	XInd = onespace( 0, N-1 );
	YInd = onespace( 0, M-1 );
	ZInd = onespace( 0, K-1 );

    // Checks to make sure grid size is suitable for refinement
	if (NumLev != 0) {
#if (dims == 3)
		for (int reg = 0; reg < NumReg; reg++) {
			// Check grid is big enough to allow embedded refinement of factor 2
			if (	(
					RefXend[reg]-RefXstart[reg]+1 < 3 || 
					RefYend[reg]-RefYstart[reg]+1 < 3 || 
					RefZend[reg]-RefZstart[reg]+1 < 3
					) || (
					(
					RefXend[reg]-RefXstart[reg]+1 == 3 ||
					RefYend[reg]-RefYstart[reg]+1 == 3 ||
					RefZend[reg]-RefZstart[reg]+1 == 3
					) && NumLev > 1 )
				) {
				cout << "Refined region is too small to support refinement" << endl;
				exit(EXIT_FAILURE);
			}
		}
#else
		for (int reg = 0; reg < NumReg; reg++) {
			// Check grid is big enough to allow embedded refinement of factor 2
			if (	(
					RefXend[reg]-RefXstart[reg]+1 < 3 || 
					RefYend[reg]-RefYstart[reg]+1 < 3
					) || (
					(
					RefXend[reg]-RefXstart[reg]+1 == 3 ||
					RefYend[reg]-RefYstart[reg]+1 == 3
					) && NumLev > 1 )
				) {
				cout << "Refined region is too small to support refinement" << endl;
				exit(EXIT_FAILURE);
			}
		}
#endif
	}

    // Checks passed so continue...
	// L0 lattice site POSITION VECTORS
	XPos = linspace( a_x + dx/2, b_x - dx/2, N );
	YPos = linspace( a_y + dy/2, b_y - dy/2, M );
	ZPos = linspace( a_z + dz/2, b_z - dz/2, K );


	// Define TYPING MATRICES
	LatTyp.resize( N*M*K );

	// Typing defined as follows:
	/*
	0 == boundary site
	1 == coarse site
	2 == fine/refined site
	3 == TL to upper (coarser) level
	4 == TL to lower (finer) level
	*/

	// Label as coarse site
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k < K; k++) {
				
				int idx = idxmap(i,j,k,M,K);
				LatTyp[idx] = 1;

			}
		}
	}


	// Add object-specific labels
	LBM_init_wall_lab();


	// Correct labelling if lower level exists
	if (NumLev > 0) {

		// Correct L0
		for(int reg = 0; reg < NumReg; reg++) {

#if (dims == 3)
			// Add TL to lower labels
			for (size_t i = RefXstart[reg]; i <= RefXend[reg]; i++) {
				for (size_t j = RefYstart[reg]; j <= RefYend[reg]; j++) {
					for (size_t k = RefZstart[reg]; k <= RefZend[reg]; k++) {

						int idx = idxmap(i,j,k,M,K);
						LatTyp[idx] = 4;

					}
				}
			}

			// Add refined labels
			for (size_t i = RefXstart[reg]+1; i <= RefXend[reg]-1; i++) {
				for (size_t j = RefYstart[reg]+1; j <= RefYend[reg]-1; j++) {
					for (size_t k = RefZstart[reg]+1; k <= RefZend[reg]-1; k++) {

						int idx = idxmap(i,j,k,M,K);
						LatTyp[idx] = 2;

					}
				}
			}

#else // 2D CASE
			// Add TL to lower labels
			for (size_t i = RefXstart[reg]; i <= RefXend[reg]; i++) {
				for (size_t j = RefYstart[reg]; j <= RefYend[reg]; j++) {
					size_t k = 0;

						int idx = idxmap(i,j,k,M,K);
						LatTyp[idx] = 4;

				}
			}

			// Add refined labels
			for (size_t i = RefXstart[reg]+1; i <= RefXend[reg]-1; i++) {
				for (size_t j = RefYstart[reg]+1; j <= RefYend[reg]-1; j++) {
					size_t k = 0;

						int idx = idxmap(i,j,k,M,K);
						LatTyp[idx] = 2;

				}
			}

#endif

		}

	}


	// Initialise L0 MACROSCOPIC quantities
	// Velocity field
	u.resize( N*M*K*dims );
	LBM_init_vel();

	// Density field
	rho.resize( N*M*K );
	LBM_init_rho();


	// Initialise L0 POPULATION matrices (f, feq)
	f.resize( N*M*K*nVels );
	feq.resize( N*M*K*nVels );

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k < K; k++) {
				for (int v = 0; v < nVels; v++) {

					// Initialise f to feq
					int idx  = idxmap(i,j,k,v,M,K,nVels);
					f[idx] = LBM_collide( i, j, k, v );

				}
			}
		}
	}
	feq = f; // Make feq = feq too

	
	// Relaxation frequency on L0
	// Assign relaxation frequency corrected for grid and time step size
	omega = 1 / ( (nu / (dt*pow(cs,2)) ) + .5 );

	cout << "L0 relaxation time = " << (1/omega) << endl;

#ifdef SOLID_ON
	// Reynolds Number (object height as length)
	Re = vecnorm(u_0x,u_0y,u_0z) * (YPos[obj_y_max]-YPos[obj_y_min]) / nu;
#else
	Re = 0;
#endif

}
	


// ***************************************************************************************************

// Method to initialise the quantities for a refined subgrid
// Assumes a volumetric setup here. Position of top left site is passed as is spacing and the relaxation frequency.
void GridObj::LBM_init_subgrid (double offsetX, double offsetY, double offsetZ, double dx0, double omega_coarse) {

	// Generate NODE NUMBERS
	XInd = onespace(0, (int)((CoarseLimsX[1] - CoarseLimsX[0] + .5)*2) );
	YInd = onespace(0, (int)((CoarseLimsY[1] - CoarseLimsY[0] + .5)*2) );
#if (dims == 3)
	ZInd = onespace(0, (int)((CoarseLimsZ[1] - CoarseLimsZ[0] + .5)*2) );
#else
	ZInd.insert(ZInd.begin(), 0); // Default for 2D
	// Reset the refined region z-limits if only 2D
	for (int i = 0; i < 2; i++) {
		CoarseLimsZ[i] = 0;
	}
#endif


	// Generate TYPING MATRICES
	// Get grid sizes
	size_t N_lim = XInd.size();
	size_t M_lim = YInd.size();
	size_t K_lim = ZInd.size();

	// Resize
	LatTyp.resize( YInd.size() * XInd.size() * ZInd.size() );

#if (dims == 3)
	// Start with TL from level above
	for (size_t i = 0; i != N_lim; ++i) {
		for (size_t j = 0; j != M_lim; ++j) {
			for (size_t k = 0; k != K_lim; ++k) {

				int idx = idxmap(i,j,k,M_lim,K_lim);
				LatTyp[idx] = 3;

			}
		}
	}

	// Check if lower grids exist
	if (NumLev > level) {

		// Add TL for next level down
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				for (size_t k = 2; k != K_lim-2; ++k) {
				
					int idx = idxmap(i,j,k,M_lim,K_lim);
					LatTyp[idx] = 4;

				}
			}
		}


		// Label rest as fine
		for (size_t i = 3; i != N_lim-3; ++i) {
			for (size_t j = 3; j != M_lim-3; ++j) {
				for (size_t k = 3; k != K_lim-3; ++k) {

					int idx = idxmap(i,j,k,M_lim,K_lim);
					LatTyp[idx] = 2;

				}
			}
		}

	} else {
		// Reached lowest level so label rest as coarse
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				for (size_t k = 2; k != K_lim-2; ++k) {

					int idx = idxmap(i,j,k,M_lim,K_lim);
					LatTyp[idx] = 1;

				}
			}
		}

	}

#else
	// Start with TL from level above
	for (size_t i = 0; i != N_lim; ++i) {
		for (size_t j = 0; j != M_lim; ++j) {
			size_t k = 0;

			int idx = idxmap(i,j,k,M_lim,K_lim);
			LatTyp[idx] = 3;

		}
	}

	// Check if lower grids exist
	if (NumLev > level) {

		// Add TL for next level down
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				size_t k = 0;
				
				int idx = idxmap(i,j,k,M_lim,K_lim);
				LatTyp[idx] = 4;
					
			}
		}


		// Label rest as fine
		for (size_t i = 3; i != N_lim-3; ++i) {
			for (size_t j = 3; j != M_lim-3; ++j) {
				size_t k = 0;

				int idx = idxmap(i,j,k,M_lim,K_lim);
				LatTyp[idx] = 2;

			}
		}

	} else {
		// Reached lowest level so label rest as coarse
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				size_t k = 0;

				int idx = idxmap(i,j,k,M_lim,K_lim);
				LatTyp[idx] = 1;

			}
		}

	}
#endif

	
	// Generate POSITION VECTORS of nodes
	// Define spacing
	dx = dx0/2;
	dy = dx0/2;
	dz = dx0/2;

	// Populate the position vectors
	XPos = linspace(offsetX - dx/2, (offsetX - dx/2) + (XInd.size() - 1) * dx, XInd.size() );
	YPos = linspace(offsetY - dy/2, (offsetY - dy/2) + (YInd.size() - 1) * dy, YInd.size() );
#if dims == 3
	ZPos = linspace(offsetZ - dz/2, (offsetZ - dz/2) + (ZInd.size() - 1) * dz, ZInd.size() );
#else
	ZPos.insert( ZPos.begin(), 1 ); // 2D default
#endif
	

	// Assign MACROSCOPIC quantities
	// Resize
	u.resize(N_lim * M_lim * K_lim * dims);
	rho.resize(N_lim * M_lim * K_lim);

	// Velocity
	LBM_init_vel( );

	// Density
	LBM_init_rho( );


	// Generate POPULATION MATRICES for lower levels
	// Resize
	f.resize(N_lim * M_lim * K_lim * nVels);
	feq.resize(N_lim * M_lim * K_lim * nVels);
			
	for (size_t i = 0; i != N_lim; ++i) {
		for (size_t j = 0; j != M_lim; ++j) {
			for (size_t k = 0; k != K_lim; ++k) {
				for (int v = 0; v < nVels; v++) {
					
					// Initialise f to feq
					int idx  = idxmap(i,j,k,v,M_lim,K_lim,nVels);
					f[idx] = LBM_collide( i, j, k, v );

				}
			}
		}
	}
	feq = f; // Set feq to feq

	// Compute relaxation time
	omega = 1 / ( ( (1/omega_coarse - .5) *2) + .5);
}