/* This file holds all the code for the core LBM operations including collision,
streaming and macroscopic calulcation.
*/

#include "stdafx.h"
#include "LBM_definitions.h"
#include "LBM_globalvars.h"

using namespace std;

// ***************************************************************************************************

// LBM algorithm applicable for both single and multi-grid
// Supply an integer r indicating from which level the algorithm is to be executed.
void LBM_multi (int r) {

	// Global variables accessible through header file references

	// Loop twice as refinement ratio per level is 2
	int count = 1;
	do {

		// Collision on Lr
		if (r == 0) {

			// Collide on whole grid
			LBM_collide(r, false);

		} else {

			// Collide on core only (excludes upper transition layer)
			LBM_collide(r, true);

		}

		// Check to see if higher level exists and on first loop
		if (r-1 >= 0 && count == 1) {

			// Explode and update the upper TL on level r+1
			LBM_explode(r);

		}

		// Check if lower level exists
		if (Nref > r) {

			// Call same routine for lower level
			LBM_multi(r+1);

			// Stream
			LBM_stream(r);
			

			// Coalesce
			LBM_coalesce(r);

		} else {

			// Stream
			LBM_stream(r);
			

		}

		// Apply boundary conditions
		LBM_boundary(r,0);

		// Update macroscopic quantities
		LBM_macro(r);
		

		// Check if on L0 and if so drop out as only need to loop once on coarsest level
		if (r == 0) {
			break;
		}

		// Increment counter
		count++;

	} while (count < 3);


}


// ***************************************************************************************************

// Collision operator
// Operates on level r and excludes the upper TL sites if core_flag set to true
void LBM_collide( int r, bool core_flag ) {

	/*
	Loop through the lattice points to compute the new distribution functions.
	Equilibrium based on:
	       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4
	*/
	
	// Declarations
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();
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
	f_new.resize( Grids[r].f.size() );
	// Initialise with current f values
	f_new = Grids[r].f;


	// Loop over lattice sites
	for (int i = i_low; i < i_high; i++) {
		for (int j = j_low; j < j_high; j++) {
			for (int k = k_low; k < k_high; k++) {

				// Get index
				int idxijk = idxmap(i,j,k,M_lim,K_lim);


				// Ignore refined sites
				if (Grids[r].LatTyp[idxijk] == 2) {
					// Do nothing as taken care of on lower level grid

				} else {


					for (int v = 0; v < nVels; v++) {
						
						// Get index
						int idxf = idxmap(i,j,k,v,M_lim,K_lim,nVels);
						
						// Get feq value by calling overload of collision function
						Grids[r].feq[idxf] = LBM_collide(i,j,k,v,r);						
						
						// Recompute distribution function f
						f_new[idxf] = (-Grids[r].omega * (Grids[r].f[idxf] - Grids[r].feq[idxf]) ) + Grids[r].f[idxf];
						
					}

				}

			}
		}
	}

	// Update f from fnew
	Grids[r].f = f_new;


}

// Overload of collision function to allow calculation of feq only for initialisation
double LBM_collide( int i, int j, int k, int v, int r ) {

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
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();
	
	// Compute the parts of the expansion for feq

		// 3D case -- IF CONTIGUOUS ONLY NEED 1 INDEX AND USE POINTER ARITHMETIC
		int idx0 = idxmap(i,j,k,0,M_lim,K_lim,dims);
		int idx1 = idxmap(i,j,k,1,M_lim,K_lim,dims);
		int idx2 = idxmap(i,j,k,2,M_lim,K_lim,dims);

#if (dims == 3)
		// Compute c_ia * u_a which is actually the dot product of c and u
		A = (c[0][v] * Grids[r].u[idx0]) + (c[1][v] * Grids[r].u[idx1]) + (c[2][v] * Grids[r].u[idx2]);

		/*
		Compute second term in the expansion
		Q_iab u_a u_b = 
		(c_x^2 - cs^2)u_x^2 + (c_y^2 - cs^2)u_y^2 + (c_z^2 - cs^2)u_z^2
		+ 2c_x c_y u_x u_y + 2c_x c_z u_x u_z + 2c_y c_z u_y u_z
		*/

		B =	(pow(c[0][v],2) - pow(cs,2)) * pow(Grids[r].u[idx0],2) + 
			(pow(c[1][v],2) - pow(cs,2)) * pow(Grids[r].u[idx1],2) +
			(pow(c[2][v],2) - pow(cs,2)) * pow(Grids[r].u[idx2],2) +
			2 * c[0][v]*c[1][v] * Grids[r].u[idx0] * Grids[r].u[idx1] + 
			2 * c[0][v]*c[2][v] * Grids[r].u[idx0] * Grids[r].u[idx2] + 
			2 * c[1][v]*c[2][v] * Grids[r].u[idx1] * Grids[r].u[idx2];
#else
		// 2D versions of the above
		A = (c[0][v] * Grids[r].u[idx0]) + (c[1][v] * Grids[r].u[idx1]);

		B =	(pow(c[0][v],2) - pow(cs,2)) * pow(Grids[r].u[idx0],2) + 
			(pow(c[1][v],2) - pow(cs,2)) * pow(Grids[r].u[idx1],2) +
			2 * c[0][v]*c[1][v] * Grids[r].u[idx0] * Grids[r].u[idx1];
#endif
	
	
	// Compute f^eq
	int idxijk = idxmap(i,j,k,M_lim,K_lim);
	feq = Grids[r].rho[idxijk] * w[v] * (1 + (A / pow(cs,2)) + (B / (2*pow(cs,4))));

	return feq;


}


// ***************************************************************************************************

// Streaming operator
// Applies periodic BCs on level 0
void LBM_stream(int r) {

	// Declarations
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();
	int dest_x, dest_y, dest_z;

	// Create temporary lattice of zeros to prevent overwriting useful populations
	vector<double> f_new( Grids[r].f.size(), 0.0 );

	// Stream one lattice site at a time
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				for (int v = 0; v < nVels; v++) {

					// Get index
					int idx_ijk = idxmap(i,j,k,M_lim,K_lim);

					// If fine site then do not stream in any direction
					if (Grids[r].LatTyp[idx_ijk] == 2) {
						break;
					}

					// Only apply periodic BCs on coarsest level
					if (r == 0) {
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
						if ( (Grids[r].LatTyp[idx_dest] == 2) ||
							( (Grids[r].LatTyp[idx_dest] == 4) 
							&& (Grids[r].LatTyp[idx_ijk] == 4) ) ) {
							// Fine -- ignore
							// TL lower level to TL lower level -- done on lower grid stream so ignore too.
					
						} else {

							// Get index
							int idx_destv = idxmap(dest_x,dest_y,dest_z,v,M_lim,K_lim,nVels);
							int idx_ijkv = idxmap(i,j,k,v,M_lim,K_lim,nVels);

							// Stream population
							f_new[idx_destv] = Grids[r].f[idx_ijkv];

						}

					}

				}

			}
		}
	}

	// Replace old grid with new grid
	Grids[r].f = f_new;

}


// ***************************************************************************************************

// Macroscopic quantity calculation
void LBM_macro(int r) {

	// Declarations
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();
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
				int idxrho = idxmap(i,j,k,M_lim,K_lim);

				for (int v = 0; v < nVels; v++) {

					// Get index
					int idxf = idxmap(i,j,k,v,M_lim,K_lim,nVels);

					// Sum up to find momentum
					fux_temp += c[0][v] * Grids[r].f[idxf];
					fuy_temp += c[1][v] * Grids[r].f[idxf];
					fuz_temp += c[2][v] * Grids[r].f[idxf];

					// Sum up to find density
					rho_temp += Grids[r].f[idxf];

				}

				// Assign density
				Grids[r].rho[idxrho] = rho_temp;

				// Get index
				int idx_ux = idxmap(i,j,k,0,M_lim,K_lim,dims);
				int idx_uy = idxmap(i,j,k,1,M_lim,K_lim,dims);
				int idx_uz = idxmap(i,j,k,2,M_lim,K_lim,dims);

				// Assign velocity
				Grids[r].u[idx_ux] = fux_temp / rho_temp;
				Grids[r].u[idx_uy] = fuy_temp / rho_temp;
#if (dims == 3)
				Grids[r].u[idx_uz] = fuz_temp / rho_temp;
#endif

			}
		}
	}


}


// ***************************************************************************************************

// Explosion operation -- called from fine level
void LBM_explode( int r ) {

	// Declarations
	int y_start, x_start, z_start;
	int M_fine = Grids[r].YPos.size();
	int K_fine = Grids[r].ZPos.size();
	int M_coarse = Grids[r-1].YPos.size();
	int K_coarse = Grids[r-1].ZPos.size();

	// Loop over coarse grid
	for (size_t i = 0; i != Grids[r-1].XPos.size(); ++i) {
		for (size_t j = 0; j != Grids[r-1].YPos.size(); ++j) {
			for (size_t k = 0; k != Grids[r-1].ZPos.size(); ++k) {

				int idx_coarse = idxmap(i,j,k,M_coarse,K_coarse);

				// If TL to lower level then partitioning required
				if (Grids[r-1].LatTyp[idx_coarse] == 4) {

					// Lookup indices for lower level. L0 start is different from other 
					// levels as it only has a single TL.
					if (r-1 == 0) {
						x_start = Grids[0].XInd[0];
						y_start = Grids[0].YInd[0];
						z_start = Grids[0].ZInd[0];
					} else {
						x_start = Grids[r-1].XInd[2];
						y_start = Grids[r-1].YInd[2];
#if (dims == 3)
						z_start = Grids[r-1].ZInd[2];
#else
						z_start = Grids[r-1].ZInd[0]; // if 2D set to default
#endif
					}
					
					// Update fine grid values according to Rohde et al.
					for (int v = 0; v < nVels; v++) {

						// Get coarse site value
						int idx_coarsef = idxmap(i,j,k,v,M_coarse,K_coarse,nVels);
						double coarse_f = Grids[r-1].f[idx_coarsef];
						
						// Find indices of fine site
						vector<int> idx_fine = indmapref(i, x_start, j, y_start, k, z_start);
						int fi = idx_fine[0];
						int fj = idx_fine[1];
						int fk = idx_fine[2];

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
						Grids[r].f[idx1] = coarse_f;
						Grids[r].f[idx2] = coarse_f;
						Grids[r].f[idx3] = coarse_f;
						Grids[r].f[idx4] = coarse_f;
						Grids[r].f[idx5] = coarse_f;
						Grids[r].f[idx6] = coarse_f;
						Grids[r].f[idx7] = coarse_f;
						Grids[r].f[idx8] = coarse_f;

#else

						// 2D Case -- square of 4 cells

						// Flatten indices
						int idx1 = idxmap(fi,	fj,		v,M_fine,nVels);
						int idx2 = idxmap(fi+1,	fj,		v,M_fine,nVels);
						int idx3 = idxmap(fi,	fj+1,	v,M_fine,nVels);
						int idx4 = idxmap(fi+1,	fj+1,	v,M_fine,nVels);
							
						// Copy coarse to fine
						Grids[r].f[idx1] = coarse_f;
						Grids[r].f[idx2] = coarse_f;
						Grids[r].f[idx3] = coarse_f;
						Grids[r].f[idx4] = coarse_f;

#endif
					}

				}

			}
		}
	}


}


// ***************************************************************************************************

// Coalesce operation -- called from coarse level
void LBM_coalesce( int r ) {

	// Declarations
	int y_start, x_start, z_start;
	int M_coarse = Grids[r].YPos.size();
	int K_coarse = Grids[r].ZPos.size();
	int M_fine = Grids[r+1].YPos.size();
	int K_fine = Grids[r+1].ZPos.size();

	// Loop over coarse grid
	for (size_t i = 0; i != Grids[r].XPos.size(); ++i) {
		for (size_t j = 0; j != Grids[r].YPos.size(); ++j) {
			for (size_t k = 0; k != Grids[r].ZPos.size(); ++k) {

				int idx_coarse = idxmap(i,j,k,M_coarse,K_coarse);

				// If TL to lower level then fetch values from lower level
				if (Grids[r].LatTyp[idx_coarse] == 4) {

					// Lookup indices for lower level. L0 start is different from other 
					// levels as it only has a single TL.
					if (r == 0) {
						x_start = Grids[0].XInd[0];
						y_start = Grids[0].YInd[0];
						z_start = Grids[0].ZInd[0];
					} else {
						x_start = Grids[r].XInd[2];
						y_start = Grids[r].YInd[2];
#if (dims == 3)			
						z_start = Grids[r].ZInd[2];
#else					
						z_start = Grids[r].ZInd[0]; // if 2D set to default
#endif
					}

					// Loop over directions
					for (int v = 0; v < nVels; v++) {

						// Get coarse site index
						int idx_coarsef = idxmap(i,j,k,v,M_coarse,K_coarse,nVels);

						// Check to see if f value is missing on coarse level
						if (Grids[r].f[idx_coarsef] != 0) {
							// If not do nothing

						} else {

							// Find indices of corresponding fine site
							vector<int> idx_fine = indmapref(i, x_start, j, y_start, k, z_start);
							int fi = idx_fine[0];
							int fj = idx_fine[1];
							int fk = idx_fine[2];
												
#if (dims == 3)

							// 3D Case -- cube of 8 cells
							
							// Now flatten indices for each cell in the cube
							int idx1 = idxmap(fi,	fj,		fk,		v,M_fine,K_fine,nVels);
							int idx2 = idxmap(fi+1,	fj,		fk,		v,M_fine,K_fine,nVels);
							int idx3 = idxmap(fi,	fj+1,	fk,		v,M_fine,K_fine,nVels);
							int idx4 = idxmap(fi+1,	fj+1,	fk,		v,M_fine,K_fine,nVels);
							int idx5 = idxmap(fi,	fj,		fk+1,	v,M_fine,K_fine,nVels);
							int idx6 = idxmap(fi+1,	fj,		fk+1,	v,M_fine,K_fine,nVels);
							int idx7 = idxmap(fi,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
							int idx8 = idxmap(fi+1,	fj+1,	fk+1,	v,M_fine,K_fine,nVels);
							
							// Average the values
							Grids[r].f[idx_coarsef] = (
								Grids[r+1].f[idx1] + Grids[r+1].f[idx2] + Grids[r+1].f[idx3] + Grids[r+1].f[idx4] 
								+ Grids[r+1].f[idx5] + Grids[r+1].f[idx6] + Grids[r+1].f[idx7] + Grids[r+1].f[idx8]
								) / pow(2, dims);

#else

							// 2D Case -- square of 4 cells

							// Flatten indices
							int idx1 = idxmap(fi,	fj,		v,M_fine,nVels);
							int idx2 = idxmap(fi+1,	fj,		v,M_fine,nVels);
							int idx3 = idxmap(fi,	fj+1,	v,M_fine,nVels);
							int idx4 = idxmap(fi+1,	fj+1,	v,M_fine,nVels);
							
							// Average the values
							Grids[r].f[idx_coarsef] = (
								Grids[r+1].f[idx1] + Grids[r+1].f[idx2] + Grids[r+1].f[idx3] + Grids[r+1].f[idx4] 
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