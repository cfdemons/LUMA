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
	
	// Wave numbers
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	double k1 = 2*PI*kn / Lx;
	double k2 = 2*PI*km / Ly;
	double k3 = 2*PI*kk / Lz;

	// Get grid sizes
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();

	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {


				// Taylor-Green -- only 2D for now...
				int idx = idxmap(i,j,k,0,M_lim,K_lim,dims);
				u[idx] = -vecnorm(u_in) * 
					cos(k1*XPos[i]) * sin(k2*YPos[j]);

				idx = idxmap(i,j,k,1,M_lim,K_lim,dims);
				u[idx] = vecnorm(u_in) * 
					(k1/k2) * sin(k1*XPos[i]) * cos(k2*YPos[j]);
#if (dims == 3)
				idx = idxmap(i,j,k,2,M_lim,K_lim,dims);
				u[idx] = 0.0;
#endif
			}
		}
	}

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

				// Wave numbers
				double Lx = b_x - a_x;
				double Ly = b_y - a_y;
				double Lz = b_z - a_z;
				double k1 = 2*PI*kn / Lx;
				double k2 = 2*PI*km / Ly;
				double k3 = 2*PI*kk / Lz;

				int idx = idxmap(i,j,k,M_lim,K_lim);

				// Taylor-Green -- only 2D for now...
				rho[idx] = rho_in - 
					(1/pow(cs,2)) * .25 * pow(vecnorm(u_in),2) * 
					(cos(2*k1*XPos[i])+cos(2*k2*YPos[j]));
			}
		}
	}

}


// ***************************************************************************************************

// Initialise level 0 grid method
void GridObj::LBM_init_grid( ) {

	// Time step
	dt = deltat;

	// Store spacing
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	dx = 2*(Lx/(2*N));
	dy = 2*(Ly/(2*M));
	dz = 2*(Lz/(2*K));

	
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
	// Indices on L0
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
	// L0 lattice site coordinates
	XPos = linspace( a_x + dx/2, b_x - dx/2, N );
	YPos = linspace( a_y + dy/2, b_y - dy/2, M );
	ZPos = linspace( a_z + dz/2, b_z - dz/2, K );


	// Initialise L0 macroscopic quantities
	// Velocity field
	u.resize( N*M*K*dims );
	LBM_init_vel();

	// Density field
	rho.resize( N*M*K );
	LBM_init_rho();

	// Initialise L0 matrices (f, feq) and typing matrix
	f.resize( N*M*K*nVels );
	feq.resize( N*M*K*nVels );
	LatTyp.resize( N*M*K );

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

	// Relaxation frequency on L0
	// Assign relaxation frequency corrected for grid and time step size
	omega = 1 / ( (nu / (dt*pow(cs,2)) ) + .5 );

	cout << "L0 relaxation time = " << (1/omega) << endl;

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