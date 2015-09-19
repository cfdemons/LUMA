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
void GridObj::LBM_init_bound_lab ( ) {

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
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();

#if defined SOLID_BLOCK_ON || defined WALLS_ON || defined INLET_ON || defined OUTLET_ON || defined WALLS_ON_2D
	int i, j, k;
#endif

#ifdef SOLID_BLOCK_ON

	// Check solid block inside domain
	if (obj_x_max > N || obj_x_min < 0 || obj_y_max > M || obj_y_min < 0 || obj_z_max > K || obj_z_min < 0) {
		// Block outside domain
		*gUtils.logfile << "Block is placed outside the domain. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Loop over object definition
	for (i = obj_x_min; i <= obj_x_max; i++) {
		for (j = obj_y_min; j <= obj_y_max; j++) {
			for (k = obj_z_min; k <= obj_z_max; k++) {

#ifdef BUILD_FOR_MPI
				int local_i, local_j, local_k;

				// Only label if the site is on current rank (allow for periodicity)
#if (dims == 3)
				if( i <= XInd[XInd.size()-2] + 1 && i >= XInd[1] - 1 &&
					j <= YInd[YInd.size()-2] + 1 && j >= YInd[1] - 1 &&
					k <= ZInd[ZInd.size()-2] + 1 && k >= ZInd[1] - 1 ) {
#else
				if( i <= XInd[XInd.size()-2] + 1 && i >= XInd[1] - 1 &&
					j <= YInd[YInd.size()-2] + 1 && j >= YInd[1] - 1 ) {
#endif
						
						// Map global indices to local indices
						local_i = i - XInd[1] + 1;
						local_j = j - YInd[1] + 1;
#if (dims == 3)
						local_k = k - ZInd[1] + 1;
#else
						local_k = k;
#endif

						LatTyp(local_i,local_j,local_k,M_lim,K_lim) = 0;
				}
#else
				LatTyp(i,j,k,M_lim,K_lim) = 0;
#endif

			}
		}
	}
#endif

#ifdef INLET_ON
	// Left hand face only

	// Check for potential singularity in BC
	if (u_0x == 1) {
		// Singularity so exit
		*gUtils.logfile << "Inlet BC fails with u_0x = 1, choose something else. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Search index vector to see if left hand wall on this rank
	for (i = 0; i < N_lim; i++ ) {
		if (XInd[i] == 0) {		// Wall found

			// Label inlet
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Inlet site
					LatTyp(i,j,k,M_lim,K_lim) = 7;

				}
			}
			break;

		}
	}
#endif

#ifdef OUTLET_ON
	// Right hand face only

	// Search index vector to see if right hand wall on this rank
	for (i = 0; i < N_lim; i++ ) {
		if (XInd[i] == N-1) {		// Wall found

			// Label outlet
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					LatTyp(i,j,k,M_lim,K_lim) = 8;

				}
			}
			break;

		}
	}
#endif

#if defined WALLS_ON && !defined WALLS_ON_2D && (dims == 3)

	// Search index vector to see if FRONT wall on this rank
	for (k = 0; k < K_lim; k++ ) {
		if (ZInd[k] == 0) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (j = 0; j < M_lim; j++) {

					LatTyp(i,j,k,M_lim,K_lim) = 0;

				}
			}
			break;

		}
	}

	// Search index vector to see if BACK wall on this rank
	for (k = 0; k < K_lim; k++ ) {
		if (ZInd[k] == K-1) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (j = 0; j < M_lim; j++) {

					LatTyp(i,j,k,M_lim,K_lim) = 0;

				}
			}
			break;

		}
	}

#endif

#if defined WALLS_ON
	// Search index vector to see if TOP wall on this rank
	for (j = 0; j < M_lim; j++ ) {
		if (YInd[j] == 0) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (k = 0; k < K_lim; k++) {

					LatTyp(i,j,k,M_lim,K_lim) = 0;

				}
			}
			break;

		}
	}

	// Search index vector to see if BOTTOM wall on this rank
	for (j = 0; j < M_lim; j++ ) {
		if (YInd[j] == M-1) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (k = 0; k < K_lim; k++) {

					LatTyp(i,j,k,M_lim,K_lim) = 0;

				}
			}
			break;

		}
	}
#endif

}

// ***************************************************************************************************

// Initialise refined labels method (L0 only)
void GridObj::LBM_init_refined_lab ( ) {

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

	// Get local grid sizes (includes overlap)
	size_t N_lim = XPos.size();
	size_t M_lim = YPos.size();
	size_t K_lim = ZPos.size();

	// Correct L0
	for(int reg = 0; reg < NumReg; reg++) {


		// Check that whole refined region exists on a single rank (and not on overlap)
		if	(
				// X starts on rank but finishes outside rank
				( ((int)RefXstart[reg] >= XInd[1] && (int)RefXstart[reg] <= XInd[XInd.size()-2]) && ((int)RefXend[reg] > XInd[XInd.size()-2]) )
			||	
				// X ends on rank but starts outside rank
				( ((int)RefXend[reg] >= XInd[1] && (int)RefXend[reg] <= XInd[XInd.size()-2]) && ((int)RefXstart[reg] < XInd[1]) )

			||

				// Y starts on rank but finishes outside rank
				( ((int)RefYstart[reg] >= YInd[1] && (int)RefYstart[reg] <= YInd[YInd.size()-2]) && ((int)RefYend[reg] > YInd[YInd.size()-2]) )
			||	
				// Y ends on rank but starts outside rank
				( ((int)RefYend[reg] >= YInd[1] && (int)RefYend[reg] <= YInd[YInd.size()-2]) && ((int)RefYstart[reg] < YInd[1]) )
#if (dims == 3)
			||

				// Z starts on rank but finishes outside rank
				( ((int)RefZstart[reg] >= ZInd[1] && (int)RefZstart[reg] <= ZInd[ZInd.size()-2]) && ((int)RefZend[reg] > ZInd[ZInd.size()-2]) )
			||	
				// Z ends on rank but starts outside rank
				( ((int)RefZend[reg] >= ZInd[1] && (int)RefZend[reg] <= ZInd[ZInd.size()-2]) && ((int)RefZstart[reg] < ZInd[1]) )
#endif
			) {

				// Throw an error
				*gUtils.logfile << "Error: Refined region starts and ends on different ranks. Exiting." << std::endl;
				exit(EXIT_FAILURE);
		}



		// Check passed...so label

#if (dims == 3)
		// Add TL to lower labels
		for (size_t i = RefXstart[reg]; i <= RefXend[reg]; i++) {
			for (size_t j = RefYstart[reg]; j <= RefYend[reg]; j++) {
				for (size_t k = RefZstart[reg]; k <= RefZend[reg]; k++) {

#ifdef BUILD_FOR_MPI
					int local_i, local_j, local_k;

					// Only label if the site is on current rank (allow for periodicity)
					if( (int)i <= XInd[XInd.size()-2] + 1 && (int)i >= XInd[1] - 1 &&
						(int)j <= YInd[YInd.size()-2] + 1 && (int)j >= YInd[1] - 1 &&
						(int)k <= ZInd[ZInd.size()-2] + 1 && (int)k >= ZInd[1] - 1 ) {
						
							// Map global indices to local indices
							local_i = i - XInd[1] + 1;
							local_j = j - YInd[1] + 1;
							local_k = k - ZInd[1] + 1;

							LatTyp(local_i,local_j,local_k,M_lim,K_lim) = 4;
					}
#else
					LatTyp(i,j,k,M_lim,K_lim) = 4;
#endif

				}
			}
		}

		// Add refined labels
		for (size_t i = RefXstart[reg]+1; i <= RefXend[reg]-1; i++) {
			for (size_t j = RefYstart[reg]+1; j <= RefYend[reg]-1; j++) {
				for (size_t k = RefZstart[reg]+1; k <= RefZend[reg]-1; k++) {

#ifdef BUILD_FOR_MPI
					int local_i, local_j, local_k;

					// Only label if the site is on current rank (allow for periodicity)
					if( (int)i <= XInd[XInd.size()-2] + 1 && (int)i >= XInd[1] - 1 &&
						(int)j <= YInd[YInd.size()-2] + 1 && (int)j >= YInd[1] - 1 &&
						(int)k <= ZInd[ZInd.size()-2] + 1 && (int)k >= ZInd[1] - 1 ) {
						
							// Map global indices to local indices
							local_i = i - XInd[1] + 1;
							local_j = j - YInd[1] + 1;
							local_k = k - ZInd[1] + 1;

							LatTyp(local_i,local_j,local_k,M_lim,K_lim) = 2;
					}
#else
					LatTyp(i,j,k,M_lim,K_lim) = 2;
#endif

				}
			}
		}




#else // 2D CASE
		// Add TL to lower labels
		for (size_t i = RefXstart[reg]; i <= RefXend[reg]; i++) {
			for (size_t j = RefYstart[reg]; j <= RefYend[reg]; j++) {
				size_t k = 0;

#ifdef BUILD_FOR_MPI
					int local_i, local_j;

					// Only label if the site is on current rank (allow for periodicity)
					if( (int)i <= XInd[XInd.size()-2] + 1 && (int)i >= XInd[1] - 1 &&
						(int)j <= YInd[YInd.size()-2] + 1 && (int)j >= YInd[1] - 1 ) {
						
							// Map global indices to local indices
							local_i = i - XInd[1] + 1;
							local_j = j - YInd[1] + 1;

							LatTyp(local_i,local_j,k,M_lim,K_lim) = 4;
					}
#else
					LatTyp(i,j,k,M_lim,K_lim) = 4;
#endif

			}
		}

		// Add refined labels
		for (size_t i = RefXstart[reg]+1; i <= RefXend[reg]-1; i++) {
			for (size_t j = RefYstart[reg]+1; j <= RefYend[reg]-1; j++) {
				size_t k = 0;

#ifdef BUILD_FOR_MPI
					int local_i, local_j;

					// Only label if the site is on current rank (allow for periodicity)
					if( (int)i <= XInd[XInd.size()-2] + 1 && (int)i >= XInd[1] - 1 &&
						(int)j <= YInd[YInd.size()-2] + 1 && (int)j >= YInd[1] - 1 ) {
						
							// Map global indices to local indices
							local_i = i - XInd[1] + 1;
							local_j = j - YInd[1] + 1;

							LatTyp(local_i,local_j,k,M_lim,K_lim) = 2;
					}
#else
					LatTyp(i,j,k,M_lim,K_lim) = 2;
#endif

			}
		}

#endif

	}

}

// ***************************************************************************************************

// Non-MPI initialise level 0 grid wrapper
void GridObj::LBM_init_grid( ) {

	// Set local size to total grid size
	std::vector<unsigned int> local_size;
	local_size.push_back(N);
	local_size.push_back(M);
	local_size.push_back(K);

	// Set default value for these MPI-specific settings
	std::vector< std::vector<unsigned int> > GlobalLimsInd;
	GlobalLimsInd.resize(1, std::vector<unsigned int>(1) );
	GlobalLimsInd[0][0] = 0;
	std::vector< std::vector<double> > GlobalLimsPos;
	GlobalLimsPos.resize(1, std::vector<double>(1) );
	GlobalLimsPos[0][0] = 0.0;

	// Call general initialiser
	LBM_init_grid( local_size, GlobalLimsInd, GlobalLimsPos );

}


// ***************************************************************************************************

// Initialise level 0 grid method
void GridObj::LBM_init_grid( std::vector<unsigned int> local_size, 
							std::vector< std::vector<unsigned int> > GlobalLimsInd, 
							std::vector< std::vector<double> > GlobalLimsPos ) {
	
	// Store physical spacing (All L0 grids are the same spacing)
	// Global dimensions
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	dx = 2*(Lx/(2*N));
	dy = 2*(Ly/(2*M));
	dz = 2*(Lz/(2*K));

	// Physical time step = physical grid spacing?
	dt = dx;
	


	////////////////////////////
	// Check input parameters //
	////////////////////////////

#if (dims == 3)
	// Check that lattice volumes are cubes in 3D
	if ( (Lx/N) != (Ly/M) || (Lx/N) != (Lz/K) ) {
		*gUtils.logfile << "Need to have lattice volumes which are cubes -- either change N/M/K or change domain dimensions. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}
	
#else
	// 2D so need square lattice cells
	if ( (Lx/N) != (Ly/M) ) {
		*gUtils.logfile << "Need to have lattice cells which are squares -- either change N/M or change domain dimensions. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

#endif

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
				*gUtils.logfile << "Refined region is too small to support refinement. Exiting." << std::endl;
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
				*gUtils.logfile << "Refined region is too small to support refinement. Exiting." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
#endif
	}


    ///////////////////
	// Checks passed //
	///////////////////

	// Get local grid sizes (includes overlap)
	size_t N_lim = local_size[0];
	size_t M_lim = local_size[1];
#if (dims == 3)
	size_t K_lim = local_size[2];
#else
	size_t K_lim = 1;
#endif
	

	// NODE NUMBERS on L0
	/* When using MPI:
	 * Node numbers should be specified in the global system.
	 * The overlapping sites have the same number as an index
	 * on its overlapping grid. The setup assumes periodicity
	 * on all edges, even where the edge of the domain is a 
	 * boundary.
	 */
#ifdef BUILD_FOR_MPI	

	// Build index vectors
	XInd = gUtils.onespace( (int)GlobalLimsInd[0][my_rank], (int)GlobalLimsInd[1][my_rank] - 1 );
	YInd = gUtils.onespace( (int)GlobalLimsInd[2][my_rank], (int)GlobalLimsInd[3][my_rank] - 1 );
	ZInd = gUtils.onespace( (int)GlobalLimsInd[4][my_rank], (int)GlobalLimsInd[5][my_rank] - 1 );

	// Add overlap indices to both ends of the vector taking into account periodicity
	XInd.insert( XInd.begin(), (XInd[0]-1 + N) % N ); XInd.insert( XInd.end(), (XInd[XInd.size()-1]+1 + N) % N );
	YInd.insert( YInd.begin(), (YInd[0]-1 + M) % M ); YInd.insert( YInd.end(), (YInd[YInd.size()-1]+1 + M) % M );
#if (dims == 3)
	ZInd.insert( ZInd.begin(), (ZInd[0]-1 + K) % K ); ZInd.insert( ZInd.end(), (ZInd[ZInd.size()-1]+1 + K) % K );
#endif

#else
	// When not builiding for MPI indices are straightforward
	XInd = gUtils.onespace( 0, N-1 );
	YInd = gUtils.onespace( 0, M-1 );
	ZInd = gUtils.onespace( 0, K-1 );
#endif



	// L0 lattice site POSITION VECTORS
	/* When using MPI:
	 * As with the indices, the overlap is assume periodic
	 * in all directions.
	 */

#ifdef BUILD_FOR_MPI

	// Create position vectors excluding overlap
	XPos = gUtils.linspace( GlobalLimsPos[0][my_rank] + dx/2, GlobalLimsPos[1][my_rank] - dx/2, N_lim-2 );
	YPos = gUtils.linspace( GlobalLimsPos[2][my_rank] + dy/2, GlobalLimsPos[3][my_rank] - dy/2, M_lim-2 );
	ZPos = gUtils.linspace( GlobalLimsPos[4][my_rank] + dz/2, GlobalLimsPos[5][my_rank] - dz/2, K_lim-2 );

	// Add overlap sites taking into account periodicity
	XPos.insert( XPos.begin(), fmod(XPos[0]-dx + Lx, Lx) ); XPos.insert( XPos.end(), fmod(XPos[XPos.size()-1]+dx + Lx, Lx) );
	YPos.insert( YPos.begin(), fmod(YPos[0]-dy + Ly, Ly) ); YPos.insert( YPos.end(), fmod(YPos[YPos.size()-1]+dy + Ly, Ly) );
#if (dims == 3)
	ZPos.insert( ZPos.begin(), fmod(ZPos[0]-dz + Lz, Ly) ); ZPos.insert( ZPos.end(), fmod(ZPos[ZPos.size()-1]+dz + Lz, Lz) );
#endif

#else
	// When not builiding for MPI positions are straightforward
	XPos = gUtils.linspace( a_x + dx/2, b_x - dx/2, N );
	YPos = gUtils.linspace( a_y + dy/2, b_y - dy/2, M );
	ZPos = gUtils.linspace( a_z + dz/2, b_z - dz/2, K );
#endif

	

	// Define TYPING MATRICES
	LatTyp.resize( N_lim*M_lim*K_lim );

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

	// Add boundary-specific labels
	LBM_init_bound_lab();

	// Refined region labelling
	if (NumLev != 0) {
		LBM_init_refined_lab();
	}

	

	// Initialise L0 MACROSCOPIC quantities
	// Velocity field
	u.resize( N_lim*M_lim*K_lim*dims );
	LBM_init_vel();

	// Density field
	rho.resize( N_lim*M_lim*K_lim );
	LBM_init_rho();

	// Cartesian force vector
	force_xyz.resize(N_lim*M_lim*K_lim*dims, 0.0);

	// Lattice force vector
	force_i.resize(N_lim*M_lim*K_lim*nVels, 0.0);


	// Initialise L0 POPULATION matrices (f, feq)
	f.resize( N_lim*M_lim*K_lim*nVels );
	feq.resize( N_lim*M_lim*K_lim*nVels );



	// Loop over grid
	for (size_t i = 0; i < N_lim; i++) {
		for (size_t j = 0; j < M_lim; j++) {
			for (size_t k = 0; k < K_lim; k++) {
				for (size_t v = 0; v < nVels; v++) {

					// Initialise f to feq
					f(i,j,k,v,M_lim,K_lim,nVels) = LBM_collide( i, j, k, v );

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

	// Relaxation frequency on L0
	// Assign relaxation frequency using lattice viscosity
	omega = 1 / ( (nu / pow(cs,2) ) + .5 );

	/* Above is valid for L0 only when dx = 1 -- general expression is:
	 * omega = 1 / ( ( (nu * dt) / (pow(cs,2)*pow(dx,2)) ) + .5 );
	 */

#ifdef USE_MRT
	double tmp[] = mrt_relax;
	for (int i = 0; i < nVels; i++) {
		mrt_omega.push_back(tmp[i]);
	}
#endif

}
	


// ***************************************************************************************************

// Method to initialise the quantities for a refined subgrid
// Assumes a volumetric setup here. Position of top left site is passed as is spacing and the relaxation frequency.
void GridObj::LBM_init_subgrid (double offsetX, double offsetY, double offsetZ, 
								double dx0, double omega_coarse, std::vector<double> mrt_omega_coarse) {

	// Generate NODE NUMBERS
	XInd = gUtils.onespace(0, (int)((CoarseLimsX[1] - CoarseLimsX[0] + .5)*2) );
	YInd = gUtils.onespace(0, (int)((CoarseLimsY[1] - CoarseLimsY[0] + .5)*2) );
#if (dims == 3)
	ZInd = gUtils.onespace(0, (int)((CoarseLimsZ[1] - CoarseLimsZ[0] + .5)*2) );
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
	XPos = gUtils.linspace(offsetX - dx/2, (offsetX - dx/2) + (XInd.size() - 1) * dx, XInd.size() );
	YPos = gUtils.linspace(offsetY - dy/2, (offsetY - dy/2) + (YInd.size() - 1) * dy, YInd.size() );
#if dims == 3
	ZPos = gUtils.linspace(offsetZ - dz/2, (offsetZ - dz/2) + (ZInd.size() - 1) * dz, ZInd.size() );
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
	
	// MRT relaxation times on the finer grid related to coarse in same way as SRT
	for (size_t i = 0; i < nVels; i++) {
		mrt_omega.push_back( 1 / ( ( (1/mrt_omega_coarse[i] - .5) *2) + .5) );
	}

#endif

}
// ***************************************************************************************************