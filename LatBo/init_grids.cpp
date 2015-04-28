/* This file contains the routines necessary to initialise the macroscopic quantities.
*/

#include "stdafx.h"
#include "LBM_definitions.h"
#include "LBM_globalvars.h"
#include <math.h>

using namespace std;

// ***************************************************************************************************

// Initialise velocity on level r
void LBM_init_vel ( int r ) {

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
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();

	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {


				// Only a 2D expression used here for now...
				int idx = idxmap(i,j,k,0,M_lim,K_lim,dims);
				Grids[r].u[idx] = -vecnorm(u_in) * 
					cos(k1*Grids[r].XPos[i]) * sin(k2*Grids[r].YPos[j]);

				idx = idxmap(i,j,k,1,M_lim,K_lim,dims);
				Grids[r].u[idx] = vecnorm(u_in) * 
					(k1/k2) * sin(k1*Grids[r].XPos[i]) * cos(k2*Grids[r].YPos[j]);
			}
		}
	}

}


// ***************************************************************************************************

// Initialise density on level r
void LBM_init_rho ( int r ) {

	// Get grid sizes
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();

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

				// Only a 2D expression used here for now...
				Grids[r].rho[idx] = rho_in - 
					(1/pow(cs,2)) * .25 * pow(vecnorm(u_in),2) * 
					(cos(2*k1*Grids[r].XPos[i])+cos(2*k2*Grids[r].YPos[j]));
			}
		}
	}

}


// ***************************************************************************************************

// Generate all the quantities for the extra (refined) grids
// Assumes a volumetric setup
void LBM_init_multi ( ) {

	// Generate NODE NUMBERS
	// Check grid core to make sure it can support a finer grid, if not throw exception
	try {

#if (dims == 3)

		if ( (Grids[0].XInd.size() < 3 || Grids[0].YInd.size() < 3 || Grids[0].ZInd.size() < 3) || 
				((Grids[0].XInd.size() == 3 || Grids[0].YInd.size() == 3 || Grids[0].ZInd.size() == 3) 
				&& Nref > 1) ) {
			throw;

#else

		if ( (Grids[0].XInd.size() < 3 || Grids[0].YInd.size() < 3) || 
				((Grids[0].XInd.size() == 3 || Grids[0].YInd.size() == 3) 
				&& Nref > 1) ) {
			throw;

#endif

		} else {

			// First level unique as doesn't need TL
			Grids[1].XInd = onespace(0, Grids[0].XInd.size()*2 -1);
			Grids[1].YInd = onespace(0, Grids[0].YInd.size()*2 -1);

#if (dims == 3)
			Grids[1].ZInd = onespace(0, Grids[0].ZInd.size()*2 -1);
#else
			Grids[1].ZInd.insert(Grids[1].ZInd.begin(), 1); // Default for 2D
#endif
				
			// Lower levels
			if (Nref > 1) {
				for (int r = 2; r <= Nref; r++) {
					// Increase refinement by factor of 2
					Grids[r].XInd = onespace(0, (Grids[r-1].XInd.size()-4)*2 -1);
					Grids[r].YInd = onespace(0, (Grids[r-1].YInd.size()-4)*2 -1);
#if (dims == 3)
					Grids[r].ZInd = onespace(0, (Grids[r-1].ZInd.size()-4)*2 -1);
#else
					Grids[r].ZInd.insert(Grids[r].ZInd.begin(), 1); // Default for 2D
#endif
				}
			}
		}
	}
	catch (int err) {
		cout << "Refined Region is too small to support refinement: error " + int2str(err) << endl;
	}

	// Correct L0 TYPING MATRIX to add suitable labels
#if (dims == 3)
	for (int i = Grids[0].XInd[0]; i <= Grids[0].XInd[Grids[0].XInd.size()-1]; i++) {
		for (int j = Grids[0].YInd[0]; j <= Grids[0].YInd[Grids[0].YInd.size()-1]; j++) {		
			for (int k = Grids[0].ZInd[0]; k <= Grids[0].ZInd[Grids[0].ZInd.size()-1]; k++) {

				int idx = idxmap(i,j,k,M,K);
				Grids[0].LatTyp[idx] = 4;

			}
		}
	}
	for (int i = Grids[0].XInd[1]; i <= Grids[0].XInd[Grids[0].XInd.size()-2]; i++) {
		for (int j = Grids[0].YInd[1]; j <= Grids[0].YInd[Grids[0].YInd.size()-2]; j++) {
			for (int k = Grids[0].ZInd[1]; k <= Grids[0].ZInd[Grids[0].ZInd.size()-2]; k++) {
		
				int idx = idxmap(i,j,k,M,K);
				Grids[0].LatTyp[idx] = 2;
			}
		}
	}

#else
	for (int i = Grids[0].XInd[0]; i <= Grids[0].XInd[Grids[0].XInd.size()-1]; i++) {
		for (int j = Grids[0].YInd[0]; j <= Grids[0].YInd[Grids[0].YInd.size()-1]; j++) {		
			int k = 0;

			int idx = idxmap(i,j,k,M,K);
			Grids[0].LatTyp[idx] = 4;

		}
	}
	for (int i = Grids[0].XInd[1]; i <= Grids[0].XInd[Grids[0].XInd.size()-2]; i++) {
		for (int j = Grids[0].YInd[1]; j <= Grids[0].YInd[Grids[0].YInd.size()-2]; j++) {
			int k = 0;
		
			int idx = idxmap(i,j,k,M,K);
			Grids[0].LatTyp[idx] = 2;
			
		}
	}
#endif

	// Generate lower level TYPING MATRICES
#if (dims == 3)
	for (int r = 1; r <= Nref; r++) {

		// Get grid sizes
		int M_lim = Grids[r].YInd.size();
		int K_lim = Grids[r].ZInd.size();

		// Resize
		Grids[r].LatTyp.resize(Grids[r].YInd.size()*Grids[r].XInd.size()*Grids[r].ZInd.size());

		// Start with TL from level above
		for (size_t i = 0; i != Grids[r].XInd.size(); ++i) {
			for (size_t j = 0; j != Grids[r].YInd.size(); ++j) {
				for (size_t k = 0; k != Grids[r].ZInd.size(); ++k) {

					int idx = idxmap(i,j,k,M_lim,K_lim);
					Grids[r].LatTyp[idx] = 3;

				}
			}
		}

		// Check if lower grids exist
		if (Nref > r) {

			// Add TL for next level down
			for (size_t i = 2; i != Grids[r].XInd.size()-2; ++i) {
				for (size_t j = 2; j != Grids[r].YInd.size()-2; ++j) {
					for (size_t k = 2; k != Grids[r].ZInd.size()-2; ++k) {
				
						int idx = idxmap(i,j,k,M_lim,K_lim);
						Grids[r].LatTyp[idx] = 4;

					}
				}
			}


			// Label rest as fine
			for (size_t i = 3; i != Grids[r].XInd.size()-3; ++i) {
				for (size_t j = 3; j != Grids[r].YInd.size()-3; ++j) {
					for (size_t k = 3; k != Grids[r].ZInd.size()-3; ++k) {

						int idx = idxmap(i,j,k,M_lim,K_lim);
						Grids[r].LatTyp[idx] = 2;

					}
				}
			}

		} else {
			// Reached lowest level so label rest as coarse
			for (size_t i = 2; i != Grids[r].XInd.size()-2; ++i) {
				for (size_t j = 2; j != Grids[r].YInd.size()-2; ++j) {
					for (size_t k = 2; k != Grids[r].ZInd.size()-2; ++k) {

						int idx = idxmap(i,j,k,M_lim,K_lim);
						Grids[r].LatTyp[idx] = 1;

					}
				}
			}

		}
	}

#else
	for (int r = 1; r <= Nref; r++) {

		// Get grid sizes
		int M_lim = Grids[r].YInd.size();
		int K_lim = Grids[r].ZInd.size();

		// Resize
		Grids[r].LatTyp.resize(Grids[r].YInd.size()*Grids[r].XInd.size()*Grids[r].ZInd.size());

		// Start with TL from level above
		for (size_t i = 0; i != Grids[r].XInd.size(); ++i) {
			for (size_t j = 0; j != Grids[r].YInd.size(); ++j) {
				size_t k = 0;

				int idx = idxmap(i,j,k,M_lim,K_lim);
				Grids[r].LatTyp[idx] = 3;

			}
		}

		// Check if lower grids exist
		if (Nref > r) {

			// Add TL for next level down
			for (size_t i = 2; i != Grids[r].XInd.size()-2; ++i) {
				for (size_t j = 2; j != Grids[r].YInd.size()-2; ++j) {
					size_t k = 0;
				
					int idx = idxmap(i,j,k,M_lim,K_lim);
					Grids[r].LatTyp[idx] = 4;
					
				}
			}


			// Label rest as fine
			for (size_t i = 3; i != Grids[r].XInd.size()-3; ++i) {
				for (size_t j = 3; j != Grids[r].YInd.size()-3; ++j) {
					size_t k = 0;

					int idx = idxmap(i,j,k,M_lim,K_lim);
					Grids[r].LatTyp[idx] = 2;

				}
			}

		} else {
			// Reached lowest level so label rest as coarse
			for (size_t i = 2; i != Grids[r].XInd.size()-2; ++i) {
				for (size_t j = 2; j != Grids[r].YInd.size()-2; ++j) {
					size_t k = 0;

					int idx = idxmap(i,j,k,M_lim,K_lim);
					Grids[r].LatTyp[idx] = 1;

				}
			}

		}
	}
#endif
		
	// Generate POSITION VECTORS of nodes
	int count, first_idx, last_idx;

	// Define spacing
	Grids[1].dx = Grids[0].dx/2;
	Grids[1].dy = Grids[0].dy/2;
	Grids[1].dz = Grids[0].dz/2;

	// L0 unusual due to no extra TL
	count = 0;
	first_idx = Grids[0].XInd[0];
	last_idx = Grids[0].XInd[ Grids[0].XInd.size()-1 ];
	for (int i = first_idx; i <= last_idx; i++) {
			
		// Call position mapping function
		double val1 = posmapref(Grids[0].XPos[i],1,'x','+');
		double val2 = posmapref(Grids[0].XPos[i],1,'x','-');
				
		// Assign to Position vectors
		Grids[1].XPos.insert( Grids[1].XPos.begin() + count, val2 );
		Grids[1].XPos.insert( Grids[1].XPos.begin() + count+1, val1 );

		// Increment counter
		count+=2;

	}

	count = 0;
	first_idx = Grids[0].YInd[0];
	last_idx = Grids[0].YInd[ Grids[0].YInd.size()-1 ];
	for (int j = first_idx; j <= last_idx; j++) {

		// Call position mapping function
		double val1 = posmapref(Grids[0].YPos[j],1,'y','+');
		double val2 = posmapref(Grids[0].YPos[j],1,'y','-');
				
		// Assign to Position vectors
		Grids[1].YPos.insert( Grids[1].YPos.begin() + count, val2 );
		Grids[1].YPos.insert( Grids[1].YPos.begin() + count+1, val1 );

		// Increment counter
		count +=2;

	}

#if (dims == 3)
	count = 0;
	first_idx = Grids[0].ZInd[0];
	last_idx = Grids[0].ZInd[ Grids[0].ZInd.size()-1 ];
	for (int k = first_idx; k <= last_idx; k++) {
		
		// Call position mapping function
		double val1 = posmapref(Grids[0].ZPos[k],1,'z','+');
		double val2 = posmapref(Grids[0].ZPos[k],1,'z','-');
				
		// Assign to Position vectors
		Grids[1].ZPos.insert( Grids[1].ZPos.begin() + count, val2 );
		Grids[1].ZPos.insert( Grids[1].ZPos.begin() + count+1, val1 );

		// Increment counter
		count+=2;

	}

#else
	Grids[1].ZPos.insert( Grids[1].ZPos.begin(), 1 ); // 2D default
#endif


	// Now do lower levels
	for (int r = 2; r <= Nref; r++ ) {

		// Spacing
		Grids[r].dx = Grids[r-1].dx/2;
		Grids[r].dy = Grids[r-1].dy/2;
		Grids[r].dz = Grids[r-1].dz/2;

		count = 0;
		first_idx = Grids[r-1].XInd[2];
		last_idx = Grids[r-1].XInd[ Grids[r-1].XInd.size()-3 ];
		for (int i = first_idx; i <= last_idx; i++) {
			
			// Call position mapping function
			double val1 = posmapref(Grids[r-1].XPos[i],r,'x','+');
			double val2 = posmapref(Grids[r-1].XPos[i],r,'x','-');
				
			// Assign to Position vectors
			Grids[r].XPos.insert( Grids[r].XPos.begin() + count, val2 );
			Grids[r].XPos.insert( Grids[r].XPos.begin() + count+1, val1 );

			// Increment counter
			count+=2;

		}
		
		count = 0;
		first_idx = Grids[r-1].YInd[2];
		last_idx = Grids[r-1].YInd[ Grids[r-1].YInd.size()-3 ];
		for (int j = first_idx; j <= last_idx; j++) {
			
			// Call position mapping function
			double val1 = posmapref(Grids[r-1].YPos[j],r,'y','+');
			double val2 = posmapref(Grids[r-1].YPos[j],r,'y','-');
				
			// Assign to Position vectors
			Grids[r].YPos.insert( Grids[r].YPos.begin() + count, val2 );
			Grids[r].YPos.insert( Grids[r].YPos.begin() + count+1, val1 );

			// Increment counter
			count+=2;

		}

#if (dims == 3)
		count = 0;
		first_idx = Grids[r-1].ZInd[2];
		last_idx = Grids[r-1].ZInd[ Grids[r-1].ZInd.size()-3 ];
		for (int k = first_idx; k <= last_idx; k++) {
			
			// Call position mapping function
			double val1 = posmapref(Grids[r-1].ZPos[k],r,'z','+');
			double val2 = posmapref(Grids[r-1].ZPos[k],r,'z','-');
				
			// Assign to Position vectors
			Grids[r].ZPos.insert( Grids[r].ZPos.begin() + count, val2 );
			Grids[r].ZPos.insert( Grids[r].ZPos.begin() + count+1, val1 );

			// Increment counter
			count+=2;

		}
#else
	Grids[1].ZPos.insert( Grids[1].ZPos.begin(), 1 ); // 2D default
#endif

	}


	// Assign MACROSCOPIC quantities
	for (int r = 1; r <= Nref; r++) {

		// Get grid sizes
		int N_lim = Grids[r].XPos.size();
		int M_lim = Grids[r].YPos.size();
		int K_lim = Grids[r].ZPos.size();

		// Resize
		Grids[r].u.resize(N_lim * M_lim * K_lim * dims);
		Grids[r].rho.resize(N_lim * M_lim * K_lim);

		// Velocity
		LBM_init_vel(r);

		// Density
		LBM_init_rho(r);

	}

	// Generate POPULATION MATRICES for lower levels
	for (int r = 1; r <= Nref; r++ ) {

		// Get grid sizes
		int N_lim = Grids[r].XPos.size();
		int M_lim = Grids[r].YPos.size();
		int K_lim = Grids[r].ZPos.size();

		// Resize
		Grids[r].f.resize(N_lim * M_lim * K_lim * nVels);
		Grids[r].feq.resize(N_lim * M_lim * K_lim * nVels);
			
		for (size_t i = 0; i != Grids[r].XInd.size(); ++i) {
			for (size_t j = 0; j != Grids[r].YInd.size(); ++j) {
				for (size_t k = 0; k != Grids[r].ZInd.size(); ++k) {
					for (int v = 0; v < nVels; v++) {


						// Initialise f to feq
						int idx  = idxmap(i,j,k,v,M_lim,K_lim,nVels);
						Grids[r].f[idx] = LBM_collide(i,j,k,v,r);

					}
				}
			}
		}
		Grids[r].feq = Grids[r].f; // Set feq to feq

		// Compute relaxation time
		Grids[r].omega = 1 / ( ( (1/Grids[r-1].omega - .5) *2) + .5);
	}

	

}