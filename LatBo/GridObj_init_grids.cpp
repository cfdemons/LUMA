/* This file contains the routines necessary to initialise the macroscopic quantities and the grids and labels.
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


	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				
				
				for (size_t d = 0; d < dims; d++) {
#ifdef NO_FLOW
					// No flow case
					u(i,j,k,d,M_lim,K_lim,dims) = 0.0;
#else
					// Normal initialise
					u(i,j,k,d,M_lim,K_lim,dims) = u_in[d];
#endif

				}

			}
		}
	}

#if defined SOLID_BLOCK_ON || defined WALLS_ON || defined WALLS_ON_2D
	// Perform solid site reset of velocity
	bc_solid_site_reset();
#endif

}


// ***************************************************************************************************

// Initialise density method
void GridObj::LBM_init_rho ( ) {

	// Get grid sizes
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();

	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {		
			for (int k = 0; k < K_lim; k++) {

				// Uniform Density
				rho(i,j,k,M_lim,K_lim) = rho_in;
			}
		}
	}

}

// ***************************************************************************************************

// Initialise wall and object labels method (L0 only)
void GridObj::LBM_init_wall_lab ( ) {

	// Typing defined as follows:
	/*
	0 == boundary site
	1 == coarse site
	2 == fine/refined site
	3 == TL to upper (coarser) level
	4 == TL to lower (finer) level
	5 == Undefined
	6 == Undefined
	7 == Inlet
	8 == Outlet
	*/

#if defined SOLID_BLOCK_ON || defined WALLS_ON || defined INLET_ON || defined OUTLET_ON || defined WALLS_ON_2D
	int i, j, k;
#endif

#ifdef SOLID_BLOCK_ON

	// Check solid block inside domain
	if (obj_x_max > N || obj_x_min < 0 || obj_y_max > M || obj_y_min < 0 || obj_z_max > K || obj_z_min < 0) {
		// Block outside domain
		std::cout << "Block is placed outside the domain. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	// Loop over object definition
	for (i = obj_x_min; i <= obj_x_max; i++) {
		for (j = obj_y_min; j <= obj_y_max; j++) {
			for (k = obj_z_min; k <= obj_z_max; k++) {

				LatTyp(i,j,k,M,K) = 0;

			}
		}
	}
#endif

#ifdef INLET_ON
	// Left hand face only

	// Check for potential singularity in BC
	if (u_0x == 1) {
		// Singularity so exit
		std::cout << "Inlet BC fails with u_0x = 1, choose something else. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	i = 0;
	for (j = 0; j < M; j++) {
		for (k = 0; k < K; k++) {

			// Inlet site
			LatTyp(i,j,k,M,K) = 7;

		}
	}
#endif

#ifdef OUTLET_ON
	// Right hand face only

	i = N-1;
	for (j = 0; j < M; j++) {
		for (k = 0; k < K; k++) {

			LatTyp(i,j,k,M,K) = 8;

		}
	}
#endif

#if defined WALLS_ON && !defined WALLS_ON_2D && (dims == 3)

	// Front
	k = 0;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {

			LatTyp(i,j,k,M,K) = 0;

		}
	}

	// Back
	k = K-1;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {

			LatTyp(i,j,k,M,K) = 0;

		}
	}

#endif

#if defined WALLS_ON && defined WALLS_ON_2D
	// Top
	j = 0;
	for (i = 0; i < N; i++) {
		for (k = 0; k < K; k++) {

			LatTyp(i,j,k,M,K) = 0;

		}
	}

	// Bottom
	j = M-1;
	for (i = 0; i < N; i++) {
		for (k = 0; k < K; k++) {

			LatTyp(i,j,k,M,K) = 0;

		}
	}
#endif

}


// ***************************************************************************************************

// Initialise level 0 grid method
void GridObj::LBM_init_grid( ) {

	// Store physical spacing
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	dx = 2*(Lx/(2*N));
	dy = 2*(Ly/(2*M));
	dz = 2*(Lz/(2*K));

	// Physical time step = physical grid spacing?
	dt = dx;
	
	// Checks to make sure grid cell dimensions are suitable
#if (dims == 3)
	// Check that lattice volumes are cubes in 3D
	if ( (Lx/N) != (Ly/M) || (Lx/N) != (Lz/K) ) {
		cout << "Need to have lattice volumes which are cubes -- either change N/M/K or change domain dimensions. Exiting." << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	
#else
	// 2D so need square lattice cells
	if ( (Lx/N) != (Ly/M) ) {
		cout << "Need to have lattice cells which are squares -- either change N/M or change domain dimensions. Exiting." << endl;
		system("pause");
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
				cout << "Refined region is too small to support refinement. Exiting." << endl;
				system("pause");
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
				cout << "Refined region is too small to support refinement. Exiting." << endl;
				system("pause");
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
	5 == Undefined
	6 == Undefined
	7 == Inlet
	8 == Outlet
	*/

	// Label as coarse site
	std::fill(LatTyp.begin(), LatTyp.end(), 1);


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

						LatTyp(i,j,k,M,K) = 4;

					}
				}
			}

			// Add refined labels
			for (size_t i = RefXstart[reg]+1; i <= RefXend[reg]-1; i++) {
				for (size_t j = RefYstart[reg]+1; j <= RefYend[reg]-1; j++) {
					for (size_t k = RefZstart[reg]+1; k <= RefZend[reg]-1; k++) {

						LatTyp(i,j,k,M,K) = 2;

					}
				}
			}

#else // 2D CASE
			// Add TL to lower labels
			for (size_t i = RefXstart[reg]; i <= RefXend[reg]; i++) {
				for (size_t j = RefYstart[reg]; j <= RefYend[reg]; j++) {
					size_t k = 0;

						LatTyp(i,j,k,M,K) = 4;

				}
			}

			// Add refined labels
			for (size_t i = RefXstart[reg]+1; i <= RefXend[reg]-1; i++) {
				for (size_t j = RefYstart[reg]+1; j <= RefYend[reg]-1; j++) {
					size_t k = 0;

						LatTyp(i,j,k,M,K) = 2;

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

	// Cartesian force vector
	force_xyz.resize(N*M*K*dims, 0.0);

	// Lattice force vector
	force_i.resize(N*M*K*nVels, 0.0);


	// Initialise L0 POPULATION matrices (f, feq)
	f.resize( N*M*K*nVels );
	feq.resize( N*M*K*nVels );



	// Loop over grid
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k < K; k++) {
				for (int v = 0; v < nVels; v++) {

					// Initialise f to feq
					f(i,j,k,v,M,K,nVels) = LBM_collide( i, j, k, v );

				}
			}
		}
	}
	feq = f; // Make feq = feq too

	// Initialise OTHER parameters

	// Compute kinematic viscosity based on target Reynolds number
#if defined IBM_ON && defined INSERT_CIRCLE_SPHERE
	// If IBM circle use diameter (in lattice units i.e. rescale wrt to physical spacing)
	nu = (ibb_r*2 / dx) * u_0x / Re;
#elif defined IBM_ON && defined INSERT_RECTANGLE_CUBOID
	// If IBM rectangle use y-dimension (in lattice units)
	nu = (ibb_l / dx) * u_0x / Re;
#elif defined SOLID_BLOCK_ON
	// Use object height
	nu = (obj_y_max - obj_y_min) * u_0x / Re;
#else
	// If no object then use domain height (in lattice units)
	nu = M * u_0x / Re;
#endif

	std::cout << "Lattice viscosity = " << nu << std::endl;

	// Relaxation frequency on L0
	// Assign relaxation frequency using lattice viscosity
	omega = 1 / ( (nu / pow(cs,2) ) + .5 );

	// Above is valid for L0 only -- general expression is
	// omega = 1 / ( ( (nu * dt) / (pow(cs,2)*pow(dx,2)) ) + .5 );

	std::cout << "L0 relaxation time = " << (1/omega) << std::endl;


#if (defined USE_MRT && dims == 3)

	// MRT relaxation times (D3Q19)
	mrt_omega.push_back(1.0);	// s0
	mrt_omega.push_back(1.19);	// s1
	mrt_omega.push_back(1.4);	// s2
	mrt_omega.push_back(1.0);	// s3
	mrt_omega.push_back(1.2);	// s4
	mrt_omega.push_back(1.0);	// s5
	mrt_omega.push_back(1.2);	// s6
	mrt_omega.push_back(1.0);	// s7
	mrt_omega.push_back(1.2);	// s8
	mrt_omega.push_back(omega);	// s9
	mrt_omega.push_back(1.4);	// s10
	mrt_omega.push_back(omega);	// s11
	mrt_omega.push_back(1.4);	// s12
	mrt_omega.push_back(omega);	// s13
	mrt_omega.push_back(omega);	// s14
	mrt_omega.push_back(omega);	// s15
	mrt_omega.push_back(1.98);	// s16
	mrt_omega.push_back(1.98);	// s17
	mrt_omega.push_back(1.98);	// s18

#elif defined USE_MRT
	// MRT relaxation times (D2Q9)
	mrt_omega.push_back(1.0);	// s0
	mrt_omega.push_back(1.4);	// s1
	mrt_omega.push_back(1.4);	// s2
	mrt_omega.push_back(1.0);	// s3
	mrt_omega.push_back(1.2);	// s4
	mrt_omega.push_back(1.0);	// s5
	mrt_omega.push_back(1.2);	// s6
	mrt_omega.push_back(omega);	// s7
	mrt_omega.push_back(omega);	// s8
	
#endif


}
	


// ***************************************************************************************************

// Method to initialise the quantities for a refined subgrid
// Assumes a volumetric setup here. Position of top left site is passed as is spacing and the relaxation frequency.
void GridObj::LBM_init_subgrid (double offsetX, double offsetY, double offsetZ, 
								double dx0, double omega_coarse, std::vector<double> mrt_omega_coarse) {

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

	// Typing defined as follows:
	/*
	0 == boundary site
	1 == coarse site
	2 == fine/refined site
	3 == TL to upper (coarser) level
	4 == TL to lower (finer) level
	5 == Undefined
	6 == Undefined
	7 == Inlet
	8 == Outlet
	*/

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

				LatTyp(i,j,k,M_lim,K_lim) = 3;

			}
		}
	}

	// Check if lower grids exist
	if (NumLev > level) {


		// Add TL for next level down
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				for (size_t k = 2; k != K_lim-2; ++k) {
				
					LatTyp(i,j,k,M_lim,K_lim) = 4;

				}
			}
		}


		// Label rest as fine
		for (size_t i = 3; i != N_lim-3; ++i) {
			for (size_t j = 3; j != M_lim-3; ++j) {
				for (size_t k = 3; k != K_lim-3; ++k) {

					LatTyp(i,j,k,M_lim,K_lim) = 2;

				}
			}
		}

	} else {


		// Reached lowest level so label rest as coarse
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				for (size_t k = 2; k != K_lim-2; ++k) {

					LatTyp(i,j,k,M_lim,K_lim) = 1;

				}
			}
		}

	}

#else
	// Start with TL from level above
	for (size_t i = 0; i != N_lim; ++i) {
		for (size_t j = 0; j != M_lim; ++j) {
			size_t k = 0;

			LatTyp(i,j,k,M_lim,K_lim) = 3;

		}
	}

	// Check if lower grids exist
	if (NumLev > level) {

		// Add TL for next level down
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				size_t k = 0;
				
				LatTyp(i,j,k,M_lim,K_lim) = 4;
					
			}
		}


		// Label rest as fine
		for (size_t i = 3; i != N_lim-3; ++i) {
			for (size_t j = 3; j != M_lim-3; ++j) {
				size_t k = 0;

				LatTyp(i,j,k,M_lim,K_lim) = 2;

			}
		}

	} else {
		// Reached lowest level so label rest as coarse
		for (size_t i = 2; i != N_lim-2; ++i) {
			for (size_t j = 2; j != M_lim-2; ++j) {
				size_t k = 0;

				LatTyp(i,j,k,M_lim,K_lim) = 1;

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
	force_xyz.resize(N_lim * M_lim * K_lim * dims);
	force_i.resize(N_lim * M_lim * K_lim * nVels);

	// Velocity
	LBM_init_vel( );

	// Density
	LBM_init_rho( );

	// Cartesian force vector
	for (unsigned int i = 0; i < N_lim*M_lim*K_lim*dims; i++) {
		force_xyz[i] = 0.0;
	}

	// Lattice force vector
	for (unsigned int i = 0; i < N_lim*M_lim*K_lim*nVels; i++) {
		force_i[i] = 0.0;
	}


	// Generate POPULATION MATRICES for lower levels
	// Resize
	f.resize(N_lim * M_lim * K_lim * nVels);
	feq.resize(N_lim * M_lim * K_lim * nVels);
	

	// Loop over grid
	for (size_t i = 0; i != N_lim; ++i) {
		for (size_t j = 0; j != M_lim; ++j) {
			for (size_t k = 0; k != K_lim; ++k) {
				for (int v = 0; v < nVels; v++) {
					
					// Initialise f to feq
					f(i,j,k,v,M_lim,K_lim,nVels) = LBM_collide( i, j, k, v );

				}
			}
		}
	}
	feq = f; // Set feq to feq

	// Compute relaxation time from coarser level assume refinement by factor of 2
	omega = 1 / ( ( (1/omega_coarse - .5) *2) + .5);

#ifdef USE_MRT
	
	// MRT relaxation times on the finer grid related to coarse in same way
	for (size_t i = 0; i < nVels; i++) {
		mrt_omega.push_back( 1 / ( ( (1/mrt_omega_coarse[i] - .5) *2) + .5) );
	}

#endif

}
// ***************************************************************************************************