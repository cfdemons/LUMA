/* This file holds all the code for applying standard LBM boundary conditions such as bounce-back, inlet and outlet. */

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include <vector>
#include "../inc/definitions.h"
#include "../inc/globalvars.h"


// ***************************************************************************************************
/*
	Boundary condition application routine.
	Supply an integer to specify which type of condition should be applied:
	0 == apply all boundary conditions simultaneously
	1 == apply solid wall conditions only
	2 == apply inlet conditions only
	3 == apply outlet conditions only
	4 == apply inlet and outlet simultaneously

	Recognised boundary label types are:
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
					************************************* Solid Sites ****************************************
					******************************************************************************************	*/
    
				if (LatTyp(i,j,k,M_lim,K_lim) == 0 && (bc_type_flag == 0 || bc_type_flag == 1) ) {

					// For each outgoing direction
					for (size_t v_outgoing = 0; v_outgoing < nVels; v_outgoing++) {

						// Identify site where the population will be streamed to
						size_t dest_x = i+c[0][v_outgoing];
						size_t dest_y = j+c[1][v_outgoing];
						size_t dest_z = k+c[2][v_outgoing];

						// If this site is off-grid then cannot apply the BC in preparation for streaming
						if (	(dest_x >= N_lim || dest_x < 0) ||
								(dest_y >= M_lim || dest_y < 0) ||
								(dest_z >= K_lim || dest_z < 0)
						   ) {
							   continue;	// Move on to next direction

						// If site is fluid site then apply BC by overwriting reverse population with expected incoming value
						} else if (LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 1 || LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 2) {

							// Get incoming direction
							size_t v_incoming = GridUtils::getOpposite(v_outgoing);
							
							// Assign incoming population to outgoing at current site
							f(i,j,k,v_outgoing,M_lim,K_lim,nVels) = f(dest_x,dest_y,dest_z,v_incoming,M_lim,K_lim,nVels);
							
						}
					}


				/*	******************************************************************************************
					************************************* Inlet Sites ****************************************
					******************************************************************************************	*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == 7 && (bc_type_flag == 0 || bc_type_flag == 2 || bc_type_flag == 4) ) {

					// !! FOR NOW ASSUME THIS IS LEFT HAND WALL !!

					// Choose option
#if (defined INLET_ON && !defined INLET_DO_NOTHING && !defined INLET_REGULARISED)
					// Apply inlet Zou-He
					bc_applyZouHe(LatTyp(i,j,k,M_lim,K_lim), i, j, k, M_lim, K_lim);

#elif (defined INLET_ON && !defined INLET_DO_NOTHING && defined INLET_REGULARISED)
					// Apply regularised BC
					bc_applyRegularised(LatTyp(i,j,k,M_lim,K_lim), i, j, k, M_lim, K_lim);
#endif

    
				/*	******************************************************************************************
					************************************ Outlet Sites ****************************************
					******************************************************************************************	*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == 8 && (bc_type_flag == 0 || bc_type_flag == 3 || bc_type_flag == 4) ) {

					// !! FOR NOW ASSUME THIS IS RIGHT HAND WALL !!

					// Apply extrapolation
					bc_applyExtrapolation(LatTyp(i,j,k,M_lim,K_lim), i, j, k, M_lim, K_lim);

				}

			// End of lattice site loop
			}
		}
	}
    
}

// ***************************************************************************************************
// Routine to apply Extrapolation outlet boundary condition
void GridObj:: bc_applyExtrapolation(int label, int i, int j, int k, int M_lim, int K_lim) {

	/* When using MPI, an outlet may appear at an edge of the grid which 
	 * is not the right-hand edge. As a remedy to this we first check
	 * whether the extrapolation has access to the values. If it does not
	 * the extrapolation is not performed. Although this will result in 
	 * the boundary condition not being applied to the outlet on one of 
	 * the MPI ranks, this outlet will be overwritten during MPI comms.
	 * Furthermore, the erroneous values will likely stream to an inlet 
	 * or other boundary site where the incoming values are assumed 
	 * unknown and are overwritten by the boundary condition at this site.
	 * Therefore, the impact of not updating these sites is anticipated to
	 * be zero.
	 */

	if (i-1 < 0 || i-2 < 0) {

		// Cannot update the outlet using extrapolation

	} else {

#if dims == 3

		// In 3D, extrapolate populations [1 7 9 15 16]
		for (size_t v = 1; v < 17; v++) {
                
			// Make all this generic in future release
			if (v == 1 || v == 7 || v == 9 || v == 15 || v == 16) {

				/*
				float y2 = (float)f(i-1,j,k,v,M_lim,K_lim,nVels);
				float y1 = (float)f(i-2,j,k,v,M_lim,K_lim,nVels);
				float x1 = 0.0;
				float x2 = (float)dx;
				float x3 = 2 * x2;
                
				float lin_m = (y2 - y1) / (x2 - x1);
				float lin_c = y1;
                
				f(i,j,k,v,M_lim,K_lim,nVels) = lin_m * x3 + lin_c;
				*/

				// Just copy value for now as the linear extrapolation doesn't work
				f(i,j,k,v,M_lim,K_lim,nVels) = f(i-1,j,k,v,M_lim,K_lim,nVels);

			}
		}

#else
		// In 2D, extrapolate populations [7 1 5]
		for (size_t v = 1; v < 8; v++) {

			if (v == 7 || v == 1 || v == 5) {
                
				/*
				float y2 = (float)f(i-1,j,k,v,M_lim,K_lim,nVels);
				float y1 = (float)f(i-2,j,k,v,M_lim,K_lim,nVels);
				float x1 = 0.0;
				float x2 = (float)dx;
				float x3 = 2 * x2;
                
				float lin_m = (y2 - y1) / (x2 - x1);
				float lin_c = y1;
                
				f(i,j,k,v,M_lim,K_lim,nVels) = lin_m * x3 + lin_c;
				*/

				// Just copy value for now as the linear extrapolation doesn't work
				f(i,j,k,v,M_lim,K_lim,nVels) = f(i-1,j,k,v,M_lim,K_lim,nVels);

			}
		}
#endif

	}

}


// ***************************************************************************************************
// Routine to apply Zou-He boundary conditions
void GridObj::bc_applyZouHe(int label, int i, int j, int k, int M_lim, int K_lim) {

	/* Zou-He velocity boundary condition computed from the following equations
	 * rho = sum ( fi )
	 * rho*ux = sum( fi * cxi )
	 * rho*uy = sum( fi * cyi )
	 * rho*uz = sum( fi * czi )
	 * (fi - feq)_in = (fi - feq)_out ------ normal to wall
	 * 
	 * + transverse momentum corrections for 3D.
	 * 
	 * 3 populations (2D) or 5 populations (3D) will be unknown for the boundary site
	 */

	// Get references for f values to make the following a bit neater and easier to read
	// but does make it slower
	IVector<double> ftmp;
	for (size_t n = 0; n < nVels; n++) {
		ftmp.push_back(f(i,j,k,n,M_lim,K_lim,nVels));
	}

#if dims == 3

	/* Implement using equations
	 * rho_in = sum( fi )
	 * rho_in * ux = (f0 + f6 + f8 + f14 + f17) - (f1 + f7 + f9 + f15 + f16)
	 * rho_in * uy = (f2 + f6 + f9 + f10 + f12) - (f3 + f7 + f8 + f11 + f13)
	 * rho_in * uz = (f4 + f10 + f13 + f14 + f16) - (f5 + f11 + f12 + f15 + f17)
	 * f0 - feq0 = f1 - feq1 (equilibrium normal to boundary)
	 * 
	 * Plus transverse momentum corrections (Hecht & Harting)
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

	/* 2D Zou-He for a left hand inlet
	 *
	 * Implement using 4 equations
     * rho_in = sum( fi )
     * rho_in * ux = (f6 + f0 + f4) - (f7 + f1 + f5)
     * rho_in * uy = (f4 + f2 + f7) - (f5 + f3 + f6)
     * f0 - feq0 = f1 - feq1 (equilibrium normal to boundary)
	 */


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

// Routine to apply Regularised boundary conditions
void GridObj::bc_applyRegularised(int label, int i, int j, int k, int M_lim, int K_lim) {

	/* According to Latt & Chopard 2008 and the cited thesis by Latt 2007 we define the regularised
	 * boundary condition as folows:
	 *
	 * 1) Apply off-equilibrium bounceback to the unknown populations.
	 * 2) Compute off-equilibrium momentum flux tensor components PI^neq_ab = sum( c_ia c_ib f^neq_i ).
	 * 3) Substitute off-equilibrium definitions.
	 * 4) Compute regularised off-equilibrium part from f^neq = (w_i / (2*cs^4)) Q_iab * PI^neq_ab
	 * where Q_iab = dot(c_i, c_i) - cs^2 delta_ab.
	 * 5) Finally replace all populations on the inlet boundary node as f_i = f^eq_i + f^neq_i
	 * 
	 */

	// Get references for f values to make the following a bit neater and easier to read
	// but does make it slower
	IVector<double> ftmp;
	for (size_t n = 0; n < nVels; n++) {
		ftmp.push_back(f(i,j,k,n,M_lim,K_lim,nVels));
	}

#if (dims == 3)
	// 3D Regularised BC for left hand inlet //

	// Compute density on the wall based on desired velocity
	double rho_wall = ( 1.0 / (1.0 - u_0x) ) * 
		(ftmp[2] + ftmp[3] + ftmp[4] + ftmp[5] + ftmp[10] + ftmp[11] + ftmp[12] + ftmp[13] + ftmp[18] +
		2*(ftmp[1] + ftmp[7] + ftmp[9] + ftmp[15] + ftmp[16]) );

	// Compute off-equilibrium momentum flux tensor components
	double Sxx, Syy, Szz, Sxy, Sxz, Syz;
	Sxx = 2*(ftmp[1] + ftmp[7] + ftmp[9] + ftmp[15] + ftmp[16]) - rho_wall*( (1.0/3.0) - u_0x + pow(u_0x,2) );
	Syy = ftmp[2] + ftmp[3] + ftmp[10] + ftmp[11] + ftmp[12] + ftmp[13] + 2*(ftmp[7] + ftmp[9]) - rho_wall*( (1.0/3.0) - (1.0/3.0)*u_0x );
	Szz = ftmp[4] + ftmp[5] + ftmp[10] + ftmp[11] + ftmp[12] + ftmp[13] + 2*(ftmp[15] + ftmp[16]) - rho_wall*( (1.0/3.0) - (1.0/3.0)*u_0x );
	Sxy = 2*(ftmp[7] - ftmp[9]);
	Sxz = 2*(ftmp[15] - ftmp[16]);
	Syz = ftmp[10] + ftmp[11] - ftmp[12] - ftmp[13];

	// Compute regularised off-equilibrium components, overwriting ftmp as we have finished with it
	for (unsigned int v = 0; v < nVels; v++) {

		ftmp[v] = (w[v] / (2*pow(cs,4))) * (
			( (pow(c[0][v],2) - pow(cs,2)) * Sxx ) + 
			( (pow(c[1][v],2) - pow(cs,2)) * Syy ) +
			( (pow(c[2][v],2) - pow(cs,2)) * Szz ) +
			( 2*c[0][v]*c[1][v]*Sxy ) +
			( 2*c[0][v]*c[2][v]*Sxz ) +
			( 2*c[1][v]*c[2][v]*Syz )
			);

	}


#else
	// 2D Regularised BC for left hand inlet //

	// Compute density on the wall based on desired velocity
	double rho_wall = ( 1.0 / (1.0 - u_0x) ) * 
		(ftmp[2] + ftmp[3] + ftmp[8] + 2*(ftmp[1] + ftmp[5] + ftmp[7]) );

	// Compute off-equilibrium momentum flux tensor components
	double Sxx, Sxy, Syy;
	Sxx = 2*(ftmp[1] + ftmp[5] + ftmp[7]) - rho_wall*( (1.0/3.0) - u_0x + pow(u_0x,2) );
	Sxy = 2*(ftmp[5] - ftmp[7]);
	Syy = ftmp[2] + ftmp[3] + 2*(ftmp[5] + ftmp[7]) - rho_wall*( (1.0/3.0) - (1.0/3.0)*u_0x );

	// Compute regularised off-equilibrium components, overwriting ftmp as we have finished with it
	for (unsigned int v = 0; v < nVels; v++) {

		ftmp[v] = (w[v] / (2*pow(cs,4))) * (
			( (pow(c[0][v],2) - pow(cs,2)) * Sxx ) + 
			( 2*c[0][v]*c[1][v]*Sxy ) +
			( (pow(c[1][v],2) - pow(cs,2)) * Syy )
			);

	}

#endif

	// Set macroscopic quantities for calculation of feq and also ensuing collision step
	rho(i,j,k,M_lim,K_lim) = rho_wall;
	u(i,j,k,0,M_lim,K_lim,dims) = u_0x;
	u(i,j,k,1,M_lim,K_lim,dims) = u_0y;
#if (dims == 3)
	u(i,j,k,2,M_lim,K_lim,dims) = u_0z;
#endif
	
	// Overwrite all the populations on the node
	for (unsigned int v = 0; v < nVels; v++) {

		// Find feq corresponding to prescribed inlet macroscopic conditions and add regularised f_neq
		f(i,j,k,v,M_lim,K_lim,nVels) = LBM_collide(i,j,k,v) + ftmp[v];

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