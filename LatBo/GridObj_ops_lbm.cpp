/* This file holds all the code for the core LBM operations including collision,
streaming and macroscopic calulcation.
*/

#include "stdafx.h"
#include "GridObj.h"
#include "definitions.h"
#include "globalvars.h"

using namespace std;


// ***************************************************************************************************

// LBM multi-grid kernel applicable for both single and multi-grid
void GridObj::LBM_multi ( ) {

	// Loop twice as refinement ratio per level is 2
	int count = 1;
	do {

		// Collision on Lr
		if (level == 0) {

			// Collide on whole grid
			LBM_collide(false);

			

		} else {

			// Collide on core only (excludes upper transition layer)
			LBM_collide(true);

		}

		// Check if lower level exists
		if (NumLev > level) {

			size_t regions = subGrid.size();
			for (size_t reg = 0; reg < regions; reg++) {

				// Explode
				LBM_explode(reg);

				// Call same routine for lower level
				subGrid[reg].LBM_multi();

			}

			// Stream
			LBM_stream();
			
			for (size_t reg = 0; reg < regions; reg++) {

				// Coalesce
				LBM_coalesce(reg);

			}
			

		} else {

			// Stream
			LBM_stream();
			

		}

		// Apply boundary conditions
		LBM_boundary(0);

		// Update macroscopic quantities
		LBM_macro();
		

		// Check if on L0 and if so drop out as only need to loop once on coarsest level
		if (level == 0) {
			break;
		}

		// Increment counter
		count++;

	} while (count < 3);

}

// ***************************************************************************************************

// Collision operator
// Excludes the upper TL sites if core_flag set to true
void GridObj::LBM_collide( bool core_flag ) {

	/*
	Loop through the lattice points to compute the new distribution functions.
	Equilibrium based on:
	       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4
	*/

	// Declarations and Grid size
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();
	int i_low, j_low, k_low, i_high, j_high, k_high;

	// Respond to core flag
	if (core_flag) {
		// Ignore TL sites to upper level
		i_low = 2; i_high = N_lim-2;
		j_low = 2; j_high = M_lim-2;
#if (dims == 3)
		k_low = 2; k_high = K_lim-2;
#else
		k_low = 0; k_high = K_lim; // if 2D set to default
#endif

	} else {
		i_low = 0; i_high = N_lim;
		j_low = 0; j_high = M_lim;
		k_low = 0; k_high = K_lim;
	}

	// Create temporary lattice to prevent overwriting useful populations and initialise with same values as
	// pre-collision f grid.
	vector<double> f_new;
	f_new.resize( f.size() );
	// Initialise with current f values
	f_new = f;


	// Loop over lattice sites
	for (int i = i_low; i < i_high; i++) {
		for (int j = j_low; j < j_high; j++) {
			for (int k = k_low; k < k_high; k++) {

				// Get index
				int idx_ijk = idxmap(i,j,k,M_lim,K_lim);


				// Ignore refined sites
				if (LatTyp[idx_ijk] == 2) {
					// Do nothing as taken care of on lower level grid

				} else {


					for (int v = 0; v < nVels; v++) {
						
						// Get index
						int idx_f = idxmap(i,j,k,v,M_lim,K_lim,nVels);
						
						// Get feq value by calling overload of collision function
						feq[idx_f] = LBM_collide( i, j, k, v );						
						
						// Recompute distribution function f
						f_new[idx_f] = (-omega * (f[idx_f] - feq[idx_f]) ) + f[idx_f];
						
					}

				}

			}
		}
	}

	// Update f from fnew
	f = f_new;

}


// Overload of collision function to allow calculation of feq only for initialisation
double GridObj::LBM_collide( int i, int j, int k, int v ) {

	/* LBGK equilibrium function is represented as:
		feq_i = rho * w_i * ( 1 + u_a c_ia / cs^2 + Q_iab u_a u_b / 2*cs^4 )
	where
		Q_iab = c_ia c_ib - cs^2 * delta_ab
	and
		delta_ab is the Kronecker delta.
	*/

	// Declare single feq value and intermediate values A and B
	double feq, A, B;

	// Other declarations
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();
	
	// Compute the parts of the expansion for feq

		// 3D case -- IF CONTIGUOUS ONLY NEED 1 INDEX AND USE POINTER ARITHMETIC?
		int idx0 = idxmap(i,j,k,0,M_lim,K_lim,dims);
		int idx1 = idxmap(i,j,k,1,M_lim,K_lim,dims);
		int idx2 = idxmap(i,j,k,2,M_lim,K_lim,dims);

#if (dims == 3)
		// Compute c_ia * u_a which is actually the dot product of c and u
		A = (c[0][v] * u[idx0]) + (c[1][v] * u[idx1]) + (c[2][v] * u[idx2]);

		/*
		Compute second term in the expansion
		Q_iab u_a u_b = 
		(c_x^2 - cs^2)u_x^2 + (c_y^2 - cs^2)u_y^2 + (c_z^2 - cs^2)u_z^2
		+ 2c_x c_y u_x u_y + 2c_x c_z u_x u_z + 2c_y c_z u_y u_z
		*/

		B =	(pow(c[0][v],2) - pow(cs,2)) * pow(u[idx0],2) + 
			(pow(c[1][v],2) - pow(cs,2)) * pow(u[idx1],2) +
			(pow(c[2][v],2) - pow(cs,2)) * pow(u[idx2],2) +
			2 * c[0][v]*c[1][v] * u[idx0] * u[idx1] + 
			2 * c[0][v]*c[2][v] * u[idx0] * u[idx2] + 
			2 * c[1][v]*c[2][v] * u[idx1] * u[idx2];
#else
		// 2D versions of the above
		A = (c[0][v] * u[idx0]) + (c[1][v] * u[idx1]);

		B =	(pow(c[0][v],2) - pow(cs,2)) * pow(u[idx0],2) + 
			(pow(c[1][v],2) - pow(cs,2)) * pow(u[idx1],2) +
			2 * c[0][v]*c[1][v] * u[idx0] * u[idx1];
#endif
	
	
	// Compute f^eq
	int idx_ijk = idxmap(i,j,k,M_lim,K_lim);
	feq = rho[idx_ijk] * w[v] * (1 + (A / pow(cs,2)) + (B / (2*pow(cs,4))));

	return feq;

}


// ***************************************************************************************************

// Streaming operator
// Applies periodic BCs on level 0
void GridObj::LBM_stream( ) {

	// Declarations
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();
	int dest_x, dest_y, dest_z;

	// Create temporary lattice of zeros to prevent overwriting useful populations
	vector<double> f_new( f.size(), 0.0 );

	// Stream one lattice site at a time
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				for (int v = 0; v < nVels; v++) {

					// Get index
					int idx_ijk = idxmap(i,j,k,M_lim,K_lim);

					// If fine site then do not stream in any direction
					if (LatTyp[idx_ijk] == 2) {
						break;
					}

					// Only apply periodic BCs on coarsest level
					if (level == 0) {
						// Compute destination coordinates
						dest_x = (i+c[0][v] + N_lim) % N_lim;
						dest_y = (j+c[1][v] + M_lim) % M_lim;
						dest_z = (k+c[2][v] + K_lim) % K_lim;
					} else {
						// No periodic BCs
						dest_x = i+c[0][v];
						dest_y = j+c[1][v];
						dest_z = k+c[2][v];
					}

					// Get destination index
					int idx_dest = idxmap(dest_x,dest_y,dest_z,M_lim,K_lim);

					// If destination off-grid, do not stream
					if (	(dest_x >= N_lim || dest_x < 0) ||
							(dest_y >= M_lim || dest_y < 0) ||
							(dest_z >= K_lim || dest_z < 0)
						) {
						// Do nothing
				
					} else {

						// Check destination site type and decide whether to stream or not
						if ( (LatTyp[idx_dest] == 2) || // Fine -- ignore
							( (LatTyp[idx_dest] == 4) && (LatTyp[idx_ijk] == 4) ) // TL lower level to TL lower level -- done on lower grid stream so ignore.
							) {
														
					
						} else {

							// Get index
							int idx_destv = idxmap(dest_x,dest_y,dest_z,v,M_lim,K_lim,nVels);
							int idx_ijkv = idxmap(i,j,k,v,M_lim,K_lim,nVels);

							// Stream population
							f_new[idx_destv] = f[idx_ijkv];

						}

					}

				}

			}
		}
	}

	// Replace old grid with new grid
	f = f_new;

}


// ***************************************************************************************************

// Macroscopic quantity calculation
void GridObj::LBM_macro( ) {

	// Declarations
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();
	double rho_temp = 0;
	double fux_temp = 0;
	double fuy_temp = 0;
	double fuz_temp = 0;

	// Loop over lattice
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// Reset temporary variables
				rho_temp = 0; fux_temp = 0; fuy_temp = 0; fuz_temp = 0;

				// Get index
				int idx_rho = idxmap(i,j,k,M_lim,K_lim);

				for (int v = 0; v < nVels; v++) {

					// Get index
					int idx_f = idxmap(i,j,k,v,M_lim,K_lim,nVels);

					// Sum up to find momentum
					fux_temp += c[0][v] * f[idx_f];
					fuy_temp += c[1][v] * f[idx_f];
					fuz_temp += c[2][v] * f[idx_f];

					// Sum up to find density
					rho_temp += f[idx_f];

				}

				// Assign density
				rho[idx_rho] = rho_temp;

				// Get index
				int idx_ux = idxmap(i,j,k,0,M_lim,K_lim,dims);
				int idx_uy = idxmap(i,j,k,1,M_lim,K_lim,dims);
				int idx_uz = idxmap(i,j,k,2,M_lim,K_lim,dims);

				// Assign velocity
				u[idx_ux] = fux_temp / rho_temp;
				u[idx_uy] = fuy_temp / rho_temp;
#if (dims == 3)
				u[idx_uz] = fuz_temp / rho_temp;
#endif

			}
		}
	}

}


// ***************************************************************************************************

// Explosion operation
void GridObj::LBM_explode( int RegionNumber ) {

	// Declarations
	int y_start, x_start, z_start;
	int N_fine = subGrid[RegionNumber].YPos.size();
	int M_fine = subGrid[RegionNumber].YPos.size();
	int K_fine = subGrid[RegionNumber].ZPos.size();
	int N_coarse = XPos.size();
	int M_coarse = YPos.size();
	int K_coarse = ZPos.size();

	// Loop over coarse grid (just region of interest)
	for (size_t i = subGrid[RegionNumber].CoarseLimsX[0]; i <= subGrid[RegionNumber].CoarseLimsX[1]; i++) {
		for (size_t j = subGrid[RegionNumber].CoarseLimsY[0]; j <= subGrid[RegionNumber].CoarseLimsY[1]; j++) {
			for (size_t k = subGrid[RegionNumber].CoarseLimsZ[0]; k <= subGrid[RegionNumber].CoarseLimsZ[1]; k++) {

				int idx_coarse = idxmap(i,j,k,M_coarse,K_coarse);

				// If TL to lower level and point belongs to region then partitioning required
				if (LatTyp[idx_coarse] == 4) {

					// Lookup indices for lower level
					x_start = subGrid[RegionNumber].CoarseLimsX[0];
					y_start = subGrid[RegionNumber].CoarseLimsY[0];
					z_start = subGrid[RegionNumber].CoarseLimsZ[0];

					// Find indices of fine site
					vector<int> idx_fine = indmapref(i, x_start, j, y_start, k, z_start);
					int fi = idx_fine[0];
					int fj = idx_fine[1];
					int fk = idx_fine[2];

					// Update fine grid values according to Rohde et al.
					for (int v = 0; v < nVels; v++) {

						// Get coarse site value
						int idx_coarsef = idxmap(i,j,k,v,M_coarse,K_coarse,nVels);
						double coarse_f = f[idx_coarsef];

#if (dims == 3)
						// 3D Case -- cube of 8 cells
							
						// Flatten indices for each cell in the cube
						int idx1 = idxmap(fi,	fj,		fk,		v,M_fine,K_fine,nVels);
						int idx2 = idxmap(fi+1,	fj,		fk,		v,M_fine,K_fine,nVels);
						int idx3 = idxmap(fi,	fj+1,	fk,		v,M_fine,K_fine,nVels);
						int idx4 = idxmap(fi+1,	fj+1,	fk,		v,M_fine,K_fine,nVels);
						int idx5 = idxmap(fi,	fj,		fk+1,	v,M_fine,K_fine,nVels);
						int idx6 = idxmap(fi+1,	fj,		fk+1,	v,M_fine,K_fine,nVels);
						int idx7 = idxmap(fi,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
						int idx8 = idxmap(fi+1,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
							
						// Copy coarse to fine
						subGrid[RegionNumber].f[idx1] = coarse_f;
						subGrid[RegionNumber].f[idx2] = coarse_f;
						subGrid[RegionNumber].f[idx3] = coarse_f;
						subGrid[RegionNumber].f[idx4] = coarse_f;
						subGrid[RegionNumber].f[idx5] = coarse_f;
						subGrid[RegionNumber].f[idx6] = coarse_f;
						subGrid[RegionNumber].f[idx7] = coarse_f;
						subGrid[RegionNumber].f[idx8] = coarse_f;

#else

						// 2D Case -- square of 4 cells

						// Flatten indices
						int idx1 = idxmap(fi,	fj,		v,M_fine,nVels);
						int idx2 = idxmap(fi+1,	fj,		v,M_fine,nVels);
						int idx3 = idxmap(fi,	fj+1,	v,M_fine,nVels);
						int idx4 = idxmap(fi+1,	fj+1,	v,M_fine,nVels);
							
						// Copy coarse to fine
						subGrid[RegionNumber].f[idx1] = coarse_f;
						subGrid[RegionNumber].f[idx2] = coarse_f;
						subGrid[RegionNumber].f[idx3] = coarse_f;
						subGrid[RegionNumber].f[idx4] = coarse_f;

#endif
					}

				}

			}
		}
	}


}


// ***************************************************************************************************

// Coalesce operation -- called from coarse level
void GridObj::LBM_coalesce( int RegionNumber ) {

	// Declarations
	int y_start, x_start, z_start;
	int N_fine = subGrid[RegionNumber].YPos.size();
	int M_fine = subGrid[RegionNumber].YPos.size();
	int K_fine = subGrid[RegionNumber].ZPos.size();
	int N_coarse = XPos.size();
	int M_coarse = YPos.size();
	int K_coarse = ZPos.size();

	// Loop over coarse grid (only region of interest)
	for (size_t i = subGrid[RegionNumber].CoarseLimsX[0]; i <= subGrid[RegionNumber].CoarseLimsX[1]; i++) {
		for (size_t j = subGrid[RegionNumber].CoarseLimsY[0]; j <= subGrid[RegionNumber].CoarseLimsY[1]; j++) {
			for (size_t k = subGrid[RegionNumber].CoarseLimsZ[0]; k <= subGrid[RegionNumber].CoarseLimsZ[1]; k++) {

				int idx_coarse = idxmap(i,j,k,M_coarse,K_coarse);

				// If TL to lower level then fetch values from lower level
				if (LatTyp[idx_coarse] == 4) {

					// Lookup indices for lower level
					x_start = subGrid[RegionNumber].CoarseLimsX[0];
					y_start = subGrid[RegionNumber].CoarseLimsY[0];
					z_start = subGrid[RegionNumber].CoarseLimsZ[0];

					// Find indices of fine site
					vector<int> idx_fine = indmapref(i, x_start, j, y_start, k, z_start);
					int fi = idx_fine[0];
					int fj = idx_fine[1];
					int fk = idx_fine[2];

					// Loop over directions
					for (int v = 0; v < nVels; v++) {

						// Get coarse site index
						int idx_coarsef = idxmap(i,j,k,v,M_coarse,K_coarse,nVels);

						// Check to see if f value is missing on coarse level
						if (f[idx_coarsef] == 0) {
																										
#if (dims == 3)
							// 3D Case -- cube of 8 cells
							
							// Flatten indices for each cell in the cube
							int idx1 = idxmap(fi,	fj,		fk,		v,M_fine,K_fine,nVels);
							int idx2 = idxmap(fi+1,	fj,		fk,		v,M_fine,K_fine,nVels);
							int idx3 = idxmap(fi,	fj+1,	fk,		v,M_fine,K_fine,nVels);
							int idx4 = idxmap(fi+1,	fj+1,	fk,		v,M_fine,K_fine,nVels);
							int idx5 = idxmap(fi,	fj,		fk+1,	v,M_fine,K_fine,nVels);
							int idx6 = idxmap(fi+1,	fj,		fk+1,	v,M_fine,K_fine,nVels);
							int idx7 = idxmap(fi,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
							int idx8 = idxmap(fi+1,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
							
							// Average the values
							f[idx_coarsef] = (
								subGrid[RegionNumber].f[idx1] + 
								subGrid[RegionNumber].f[idx2] + 
								subGrid[RegionNumber].f[idx3] + 
								subGrid[RegionNumber].f[idx4] +
								subGrid[RegionNumber].f[idx5] + 
								subGrid[RegionNumber].f[idx6] + 
								subGrid[RegionNumber].f[idx7] + 
								subGrid[RegionNumber].f[idx8]
								) / pow(2, dims);

#else

							// 2D Case -- square of 4 cells

							// Flatten indices
							int idx1 = idxmap(fi,	fj,		v,M_fine,nVels);
							int idx2 = idxmap(fi+1,	fj,		v,M_fine,nVels);
							int idx3 = idxmap(fi,	fj+1,	v,M_fine,nVels);
							int idx4 = idxmap(fi+1,	fj+1,	v,M_fine,nVels);
							
							// Average the values
							f[idx_coarsef] = (
								subGrid[RegionNumber].f[idx1] + 
								subGrid[RegionNumber].f[idx2] + 
								subGrid[RegionNumber].f[idx3] + 
								subGrid[RegionNumber].f[idx4] 
								) / pow(2, dims);

#endif
						}
					
					}

				}

			}
		}
	}

}

// ***************************************************************************************************