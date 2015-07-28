/* This file holds all the code for applying standard LBM boundary conditions such as bounce-back, inlet and outlet. */

#include "stdafx.h"
#include "GridObj.h"
#include <vector>
#include "definitions.h"
#include "globalvars.h"

size_t getOpposite(size_t direction);

// ***************************************************************************************************
/*
	Boundary condition application routine
	Supply an integer r indicating from which level the algorithm is to be executed plus a flag to specify
	which type of condition should be applied.
	0 == apply solid site (bounceback)
	1 == apply solid site (Zou-He)
	2 == apply inlet and outlet (Zou-He & extrapolation)

	Boundary label types are:
	0 == solid site (no-slip)
	7 == inlet site
	8 == outlet site
*/
void GridObj::LBM_boundary (int bc_type_flag) {

	// Get grid sizes
	size_t N_lim = XPos.size();
	size_t M_lim = YPos.size();
	size_t K_lim = ZPos.size();


	// Loop over grid, identify BC required & apply BC
	for (size_t i = 0; i < N_lim; i++) {
		for (size_t j = 0; j < M_lim; j++) {
			for (size_t k = 0; k < K_lim; k++) {


				/*	******************************************************************************************
					************************************** BOUNCEBACK ****************************************
					*********************************** for solid sites **************************************
					******************************************************************************************
				*/
    
				if (LatTyp(i,j,k,M_lim,K_lim) == 0 && bc_type_flag == 0) {

					// For each direction
					for (size_t v = 0; v < nVels; v++) {

						// Identify where it streams to
						size_t src_x = i+c[0][v];
						size_t src_y = j+c[1][v];
						size_t src_z = k+c[2][v];

						// If this site is off-grid or another boundary site then retain current value 
						if (	(src_x >= N_lim || src_x < 0) ||
								(src_y >= M_lim || src_y < 0) ||
								(src_z >= K_lim || src_z < 0) ||
								(
								LatTyp(src_x,src_y,src_z,M_lim,K_lim) == 0 || 
								LatTyp(src_x,src_y,src_z,M_lim,K_lim) == 7 || 
								LatTyp(src_x,src_y,src_z,M_lim,K_lim) == 8) 
								) {

						// If site is fluid site then need to apply reverse
						} else if (LatTyp(src_x,src_y,src_z,M_lim,K_lim) == 1 || LatTyp(src_x,src_y,src_z,M_lim,K_lim) == 2) {

							// Get reverse direction
							size_t v_rev = getOpposite(v);
							
							// Assign reverse velocity to current site
							f(i,j,k,v,M_lim,K_lim,nVels) = f(src_x,src_y,src_z,v_rev,M_lim,K_lim,nVels);
							
						}
					}


				/*	******************************************************************************************
					*************************************** ZOU - HE *****************************************
					*********************************** for inlet sites **************************************
					******************************************************************************************
				*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == 7 && bc_type_flag == 2) {

					// !! FOR NOW ASSUME THIS IS LEFT HAND WALL !!

					// Apply inlet Zou-He
					bc_applyZouHe(LatTyp(i,j,k,M_lim,K_lim), i, j, k, M_lim, K_lim);
    

				/*	******************************************************************************************
					************************************* EXTRAPOLATION **************************************
					************************************ for outlet sites ************************************
					******************************************************************************************
				*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == 8 && bc_type_flag == 2) {

					// !! FOR NOW ASSUME THIS IS RIGHT HAND WALL !!
#if dims == 3

					// In 3D, extrapolate populations [1 7 9 15 16]
					for (size_t v = 1; v < 17; v++) {
                
						// Make all this generic in future release
						if (v == 1 || v == 7 || v == 9 || v == 15 || v == 16) {

							float y2 = (float)f(i-1,j,k,v,M_lim,K_lim,nVels);
							float y1 = (float)f(i-2,j,k,v,M_lim,K_lim,nVels);
							float x1 = 0.0;
							float x2 = (float)dx;
							float x3 = 2 * x2;
                
							float lin_m = (y2 - y1) / (x2 - x1);
							float lin_c = y1;
                
							f(i,j,k,v,M_lim,K_lim,nVels) = lin_m * x3 + lin_c;
						}
					}

#else
					// In 2D, extrapolate populations [7 1 5]
					for (size_t v = 1; v < 8; v++) {

						if (v == 7 || v == 1 || v == 5) {
                
							float y2 = (float)f(i-1,j,k,v,M_lim,K_lim,nVels);
							float y1 = (float)f(i-2,j,k,v,M_lim,K_lim,nVels);
							float x1 = 0.0;
							float x2 = (float)dx;
							float x3 = 2 * x2;
                
							float lin_m = (y2 - y1) / (x2 - x1);
							float lin_c = y1;
                
							f(i,j,k,v,M_lim,K_lim,nVels) = lin_m * x3 + lin_c;

						}
					}
#endif
                


				/*	******************************************************************************************
					*************************************** ZOU - HE *****************************************
					***********************************  for solid sites *************************************
					******************************************************************************************
				*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == 0 && bc_type_flag == 1) {

					// IMPLEMENT THIS LATER AS NOT TRIVIAL -- SEE HARTING & HECHT PAPER + ERRATA

				}

			}
		}
	}
    
}

// ***************************************************************************************************

// Routine to apply Zou-He boundary conditions
void GridObj::bc_applyZouHe(int label, int i, int j, int k, int M_lim, int K_lim) {

	/* Zou-He velocity boundary condition computed from the following equations
	rho = sum ( fi )
	rho*ux = sum( fi * cxi )
	rho*uy = sum( fi * cyi )
	rho*uz = sum( fi * czi )
	(fi - feq)_in = (fi - feq)_out ------ normal to wall
	
	+ transverse momentum corrections for 3D.

	3 populations (2D) or 5 populations (3D) will be unknown for the boundary site
	*/

	// Get references for f values to make the following a bit neater and easier to read
	// but does make it a lot slower to execute
	ivector<double> ftmp;
	for (size_t n = 0; n < nVels; n++) {
		ftmp.push_back(f(i,j,k,n,M_lim,K_lim,nVels));
	}

#if dims == 3

	/*	Implement using equations
		rho_in = sum( fi )
		rho_in * ux = (f0 + f6 + f8 + f14 + f17) - (f1 + f7 + f9 + f15 + f16)
		rho_in * uy = (f2 + f6 + f9 + f10 + f12) - (f3 + f7 + f8 + f11 + f13)
		rho_in * uz = (f4 + f10 + f13 + f14 + f16) - (f5 + f11 + f12 + f15 + f17)
		f0 - feq0 = f1 - feq1 (equilibrium normal to boundary)

		Plus transverse momentum corrections (Hecht & Harting)
	*/
	            
	// Find density on wall corresponding to given velocity
    double rho_w = (1.0 / (1.0 - u_0x)) * ( (
        ftmp[18] + ftmp[2] + ftmp[3] + ftmp[4] + ftmp[5] + ftmp[10] + ftmp[11] + ftmp[12] + ftmp[13] 
		) + 2.0 * (
        ftmp[1] + ftmp[7] + ftmp[9] + ftmp[15] + ftmp[16]
		) );

	// Find f0
	ftmp[0] = ftmp[1] + (1.0/3.0) * rho_w * u_0x;

	// Compute transverse momentum corrections
	double Nxy = 0.5 * ( ftmp[2] + ftmp[10] + ftmp[12] - ( ftmp[3] + ftmp[11] + ftmp[13] ) ) - (1.0/3.0) * rho_w * u_0y;
	double Nxz = 0.5 * ( ftmp[4] + ftmp[10] + ftmp[13] - ( ftmp[5] + ftmp[11] + ftmp[12] ) ) - (1.0/3.0) * rho_w * u_0z;

	// Compute f6, f9, f14 and f18
	ftmp[6] = ftmp[7] + (2.0 * w[7] / pow(cs,2)) * rho_w * (u_0x + u_0y) - Nxy;
	ftmp[8] = ftmp[9] + (2.0 * w[9] / pow(cs,2)) * rho_w * (u_0x - u_0y) + Nxy;
	ftmp[14] = ftmp[15] + (2.0 * w[15] / pow(cs,2)) * rho_w * (u_0x + u_0z) - Nxz;
	ftmp[17] = ftmp[16] + (2.0 * w[16] / pow(cs,2)) * rho_w * (u_0x - u_0z) + Nxz;


#else

	// 2D Zou-He for a left hand inlet

	// Implement using 4 equations
    // rho_in = sum( fi )
    // rho_in * ux = (f6 + f0 + f4) - (f7 + f1 + f5)
    // rho_in * uy = (f4 + f2 + f7) - (f5 + f3 + f6)
    // f0 - feq0 = f1 - feq1 (equilibrium normal to boundary)

    // Find density on wall corresponding to given velocity
    double rho_w = (1.0 / (1.0 - u_0x)) * 
		( ftmp[8] + ftmp[2] + ftmp[3] + 
			2.0 * (
			ftmp[7] + ftmp[1] + ftmp[5]
		) );
            
    // Find f0 using equations above
    ftmp[0] = ftmp[1] + (2.0/3.0) * rho_w * u_0x;
            
    // Find f4 using equations above
    ftmp[4] = 0.5 * ( (rho_w * u_0x) - 
        (ftmp[0] + ftmp[2]) + ftmp[1] + 2.0*ftmp[5] + ftmp[3] );
            
    // Find f6 using equations above
    ftmp[6] = 0.5 * ( (rho_w * u_0x) - 
        (ftmp[0] + ftmp[3]) + ftmp[2] + 2.0*ftmp[7] + ftmp[1] );

#endif

	// Apply new f values to grid
	for (size_t n = 0; n < nVels; n++) {
		f(i,j,k,n,M_lim,K_lim,nVels) = ftmp[n];
	}

}


// ***************************************************************************************************

// Routine to reset the velocity at solid sites to zero
void GridObj::bc_solid_site_reset( ) {

	// Get grid sizes
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();


	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				
				// Reset solid site velocities to zero
				if (LatTyp(i,j,k,M_lim,K_lim) == 0) {
					
					u(i,j,k,0,M_lim,K_lim,dims) = 0.0;
					u(i,j,k,1,M_lim,K_lim,dims) = 0.0;
#if (dims == 3)
					u(i,j,k,2,M_lim,K_lim,dims) = 0.0;
#endif
				}

			}
		}
	}

}


// ***************************************************************************************************
// ***************************************************************************************************
// Routine to compute the opposite direction of the one supplied based on D2Q9 or D3Q19 numbering
size_t getOpposite(size_t direction) {

	size_t direction_opposite;

	// If rest particle then opposite is simply itself
	if (direction == nVels-1) {
		
		direction_opposite = direction;

	} else {
		
		/*	If direction is even, then opposite is direction+1.
			If direction is odd, then opposite is direction-1.
			e.g. direction 0 (+x direction) has opposite 1 (-x direction) --> +1
			however, direction 1 has opposite 0 --> -1
			Hence we can add (-1 ^ direction) so it alternates between +/-1
		*/

		direction_opposite = direction + (int)pow(-1,direction);

	}

	return direction_opposite;

}

// ***************************************************************************************************
// ***************************************************************************************************