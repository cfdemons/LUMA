/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) 2015, 2016
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * distribution without written consent.
 *
 */

/* This file holds all the code for applying standard LBM boundary conditions such as bounce-back, inlet and outlet. */

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/definitions.h"
#include "../inc/globalvars.h"
#include "../inc/BFLBody.h"
#include "../inc/ObjectManager.h"
#include <numeric>


// ***************************************************************************************************
/*
	Boundary condition application routine.
	Supply an integer to specify which type of condition should be applied:
	0 == apply all boundary conditions simultaneously
	1 == apply solid wall and symmetry conditions only
	2 == apply inlet conditions only
	3 == apply outlet conditions only
	4 == apply inlet and outlet simultaneously
	5 == apply BFL conditions

	Recognised boundary label types are:
	0 == solid site (no-slip)
	5 == BFL site
	6 == symmetry site (free-slip)
	7 == inlet site
	8 == outlet site
	10 == solid site (no-slip -- site refined)
	16 == symmetry site (free-slip -- site refined)
	17 == inlet site (site refined)
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
    
				// Apply to solid sites
				if ( (LatTyp(i,j,k,M_lim,K_lim) == 0 || LatTyp(i,j,k,M_lim,K_lim) == 10)
					&& (bc_type_flag == 0 || bc_type_flag == 1) ) {

					// Apply half-way bounce-back
					bc_applyBounceBack(LatTyp(i,j,k,M_lim,K_lim),i,j,k,N_lim,M_lim,K_lim);

				} else if ( (LatTyp(i,j,k,M_lim,K_lim) == 6 || LatTyp(i,j,k,M_lim,K_lim) == 16)
					&& (bc_type_flag == 0 || bc_type_flag == 1) ) {

					// Apply half-way specular reflection
					bc_applySpecReflect(LatTyp(i,j,k,M_lim,K_lim),i,j,k,N_lim,M_lim,K_lim);


				/*	******************************************************************************************
					************************************* Inlet Sites ****************************************
					******************************************************************************************	*/
    
				} else if ( (LatTyp(i,j,k,M_lim,K_lim) == 7 || LatTyp(i,j,k,M_lim,K_lim) == 17) 
					&& (bc_type_flag == 0 || bc_type_flag == 2 || bc_type_flag == 4) ) {

					// Choose option
#if (defined INLET_ON && defined INLET_REGULARISED)
					// Apply regularised BC
					bc_applyRegularised(LatTyp(i,j,k,M_lim,K_lim), i, j, k, N_lim, M_lim, K_lim);
#endif
					// Do-Nothing BC is applied by default in a different way

    
				/*	******************************************************************************************
					************************************ Outlet Sites ****************************************
					******************************************************************************************	*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == 8 && (bc_type_flag == 0 || bc_type_flag == 3 || bc_type_flag == 4) ) {

					// !! RIGHT HAND WALL ONLY !!
#ifdef OUTLET_ON
					// Apply extrapolation
					bc_applyExtrapolation(LatTyp(i,j,k,M_lim,K_lim), i, j, k, M_lim, K_lim);
#endif



				/*	******************************************************************************************
					************************************** BFL Sites *****************************************
					******************************************************************************************	*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == 5 && (bc_type_flag == 0 || bc_type_flag == 5) ) {

					bc_applyBfl(i,j,k);
						

				}	// End of BC type control flow
			
			}	// End of lattice site loop
		}
	}
    
}


// ***************************************************************************************************
// Routine to apply half-way bounceback boundary condition
void GridObj::bc_applyBounceBack(int label, int i, int j, int k, int N_lim, int M_lim, int K_lim) {

	int dest_x, dest_y, dest_z;

	// For each outgoing direction
	for (size_t v_outgoing = 0; v_outgoing < nVels; v_outgoing++) {

		// Identify site where the population will be streamed to //

		// Consider periodic stream (serial/non-MPI)
#if (defined PERIODIC_BOUNDARIES && !defined BUILD_FOR_MPI)

		dest_x = (i+c[0][v_outgoing] + N_lim) % N_lim;
		dest_y = (j+c[1][v_outgoing] + M_lim) % M_lim;
		dest_z = (k+c[2][v_outgoing] + K_lim) % K_lim;
#else
		dest_x = i+c[0][v_outgoing];
		dest_y = j+c[1][v_outgoing];
		dest_z = k+c[2][v_outgoing];
#endif

		// If this site is off-grid then cannot apply the BC in preparation for streaming
		if (	(dest_x >= N_lim || dest_x < 0) ||
				(dest_y >= M_lim || dest_y < 0) ||
				(dest_z >= K_lim || dest_z < 0)
			) {
				continue;	// Move on to next direction

#ifdef BUILD_FOR_MPI
		// When using MPI, equivalent to off-grid is when destination is in 
		// periodic recv layer with periodic boundaries disabled.
		} else if ( GridUtils::isOnRecvLayer(XPos[dest_x],YPos[dest_y],ZPos[dest_z]) 
			&& GridUtils::isOverlapPeriodic(dest_x,dest_y,dest_z,*this) ) {

#if (!defined PERIODIC_BOUNDARIES)
			continue;	// Periodic boundaries disabled so do not try to apply BC
#endif

			/* If periodic boundaries are enabled then this is a valid stream 
			 * so allow flow through to next section */

#endif	// BUILD_FOR_MPI

		}	// End of exclusions

		/* Not been filtered by above exclusions so try to apply boundary condition. */

		// Only apply if destination is a fluid site or precomputed inlet
		if (	LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 1 || 
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 3 ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 4 ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 7
			) {

			// Get incoming direction
			size_t v_incoming = GridUtils::getOpposite(v_outgoing);
							
			// Overwriting outgoing population with expected incoming value
			f(i,j,k,v_outgoing,M_lim,K_lim,nVels) = f(dest_x,dest_y,dest_z,v_incoming,M_lim,K_lim,nVels);
							
		}
	}

}

// ***************************************************************************************************
// Routine to apply Extrapolation outlet boundary condition
void GridObj::bc_applyExtrapolation(int label, int i, int j, int k, int M_lim, int K_lim) {

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
// Routine to apply Regularised boundary conditions
void GridObj::bc_applyRegularised(int label, int i, int j, int k, int N_lim, int M_lim, int K_lim) {

	int dest_x, dest_y, dest_z, n;

	/** To allow a generalised application of the regularised BC we need to 
	 * implement the following:
	 * 1. Check normal directions to find orientation of normal wall;
	 * 2. Store the unknown and known directions in two vectors
	 * 3. Apply the boundary condition.
	 * Note: for the corners, the unknown distributions are too numerous 
	 * to evaluate the density so must use extrapolation for the 4 corners (2D) 
	 * or 8 corners (3D) from the bulk flow to find the density. Hence, we 
	 * handle this case after the normal cases.
	 */

	// Get references for f values to make the following a bit neater and easier to read
	// but does make it slower
	IVector<double> ftmp;
	for (n = 0; n < nVels; n++) {
		ftmp.push_back(f(i,j,k,n,M_lim,K_lim,nVels));
	}	

	// Loop to find orientation of wall (do not include rest direction)
	for (size_t v_outgoing = 0; v_outgoing < nVels - 1; v_outgoing++) {

		// Identify site in this direction //
		dest_x = i+c[0][v_outgoing];
		dest_y = j+c[1][v_outgoing];
		dest_z = k+c[2][v_outgoing];

		// If this site is off-grid then cannot determine orientation of wall
		if (GridUtils::isOffGrid(dest_x,dest_y,dest_z,N_lim,M_lim,K_lim,*this)) {

			continue;

		}

		// Only apply if destination is a fluid site or a solid site (i.e. an inlet site which feeds the domain)
		if (	LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 1 || 
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 3 ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 4 ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 0 ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 10
			) {

			// Declarations
			double rho_wall;
			double f_plus = 0, f_zero = 0;
			std::vector<double> x_plus, x_minus, y_plus, y_minus, z_plus, z_minus;
			double u_normal;
			int normal_dir;

			// First 4 (2D) or 6 (3D) are normals so apply straight wall version of density
			if (v_outgoing < 2 * dims) {

				/* According to Latt & Chopard 2008 and the cited thesis by Latt 2007 we define the regularised
				 * boundary condition as folows:
				 *
				 * 1) Apply off-equilibrium bounceback to the unknown populations.
				 * 2) Compute off-equilibrium stress components PI^neq_ab = sum( c_ia c_ib f^neq_i ).
				 * 3) Substitute off-equilibrium definitions.
				 * 4) Compute regularised off-equilibrium part from f^neq = (w_i / (2*cs^4)) Q_iab * PI^neq_ab
				 * where Q_iab = dot(c_i, c_i) - cs^2 delta_ab.
				 * 5) Finally replace all populations on the inlet boundary node as f_i = f^eq_i + f^neq_i
				 * 
				 */				

				// Find nature of normal and assign normal velocity
				for (n = 0; n < dims; n++) {

					// Is it normal x, y or z vector?
					if (abs(c[n][v_outgoing]) == 1) {
						normal_dir = n;

						// Get appropriate normal velocity (+ve u_normal = incoming flow)
						switch (n) {
						case 0:
							u_normal = c[normal_dir][v_outgoing] * u_0x;
							break;
						case 1:
							u_normal = c[normal_dir][v_outgoing] * u_0y;
							break;
						case 2:
							u_normal = c[normal_dir][v_outgoing] * u_0z;
							break;
						}

						break;
					}
				}


				// Loop through direction vectors and compute the relevant momenta
				// from the known (outgoing) populations
				for (n = 0; n < nVels; n++) {

					// Check against normal
					if (c[normal_dir][n] == -c[normal_dir][v_outgoing]) {
			
						// Add to outgoing momentum
						f_plus += ftmp[n];

					} else if (c[normal_dir][n] == 0) {

						// Tangential and rest momentum
						f_zero += ftmp[n];

					}

				}

				// Compute wall density from above momenta
				rho_wall = (1.0 / (1.0 - u_normal)) * (2.0 * f_plus + f_zero);


			// Otherwise, it must be a corner/edge
			} else {

				// Use second-order extrapolation along a line to approximate density at corner 
				// as insufficient populations to compute it using the usual method
				int dest_x2 = dest_x + c[0][v_outgoing];
				int dest_y2 = dest_y + c[1][v_outgoing];
				int dest_z2 = dest_z + c[2][v_outgoing];
				rho_wall = 2 * rho(dest_x,dest_y,dest_z,M_lim,K_lim) - rho(dest_x2,dest_y2,dest_z2,M_lim,K_lim);

			}



			// Now density has been found continue... //

			// Set macroscopic quantities to desired values
			rho(i,j,k,M_lim,K_lim) = rho_wall;
			u(i,j,k,0,M_lim,K_lim,dims) = u_0x;
			u(i,j,k,1,M_lim,K_lim,dims) = u_0y;
#if (dims == 3)
			u(i,j,k,2,M_lim,K_lim,dims) = u_0z;
#endif

			// Update feq to match the desired density and velocity and
			// apply off-equilibrium bounce-back to unknown populations
			for (n = 0; n < nVels; n++) {
				feq(i,j,k,n,M_lim,K_lim,nVels) = LBM_collide(i,j,k,n,M_lim,K_lim);


				// Different "if conditions" for corners/edges and normal walls
				if (v_outgoing < 2 * dims) {
					// Unknowns on normal wall share normal component
					if (c[normal_dir][n] == c[normal_dir][v_outgoing]) {
						ftmp[n] = ftmp[GridUtils::getOpposite(n)] - 
							(feq(i,j,k,GridUtils::getOpposite(n),M_lim,K_lim,nVels) - feq(i,j,k,n,M_lim,K_lim,nVels));
					}
				} else {
					// Unknowns on corners/edges are ones which do not share any component of outgoing vector
					if (c[0][n] != c[0][v_outgoing] && c[1][n] != c[1][v_outgoing]
#if (dims == 3)
					&& c[2][n] != c[2][v_outgoing]
#endif
					) {
						ftmp[n] = ftmp[GridUtils::getOpposite(n)] - 
							(feq(i,j,k,GridUtils::getOpposite(n),M_lim,K_lim,nVels) - feq(i,j,k,n,M_lim,K_lim,nVels));
					}
				}
			}

			// Get off-equilibrium stress components
			double Sxx = 0, Syy = 0, Sxy = 0;	// 2D & 3D
			double Szz = 0, Sxz = 0, Syz = 0;	// Just 3D
			for (n = 0; n < nVels; n++) {
				Sxx += c[0][n] * c[0][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,nVels));
				Syy += c[1][n] * c[1][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,nVels));
				Szz += c[2][n] * c[2][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,nVels));
				Sxy += c[0][n] * c[1][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,nVels));
				Sxz += c[0][n] * c[2][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,nVels));
				Syz += c[1][n] * c[2][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,nVels));
			}

			// Compute regularised non-equilibrium components and add to feq to get new populations
			for (int v = 0; v < nVels; v++) {

				f(i,j,k,v,M_lim,K_lim,nVels) = 

					LBM_collide(i,j,k,v,M_lim,K_lim) +
						
					(w[v] / (2 * pow(cs,4))) * (
					( (pow(c[0][v],2) - pow(cs,2)) * Sxx ) + 
					( (pow(c[1][v],2) - pow(cs,2)) * Syy ) +
					( (pow(c[2][v],2) - pow(cs,2)) * Szz ) +
					( 2 * c[0][v] * c[1][v] * Sxy ) +
					( 2 * c[0][v] * c[2][v] * Sxz ) +
					( 2 * c[1][v] * c[2][v] * Syz )
					);
			}

			// Do not check anymore directions as BC applied but increment 
			// so we can check whether BC has been applied
			v_outgoing++;
			break;

		}

		/** If it reaches here and has gone through all the directions then the wall orientation cannot be found
		 * based on the fact that it cannot find an adjacent site within the centre of the domain on which to 
		 * base the BC. In this case it must be a buffered layer of inlet sites (the "second row" of sites you get
		 * when embedding a sub-grid in the inlet) so just set it to default values as it doesn't affect the domain
		 * anyway. */
		if (v_outgoing == nVels) {

			// Set macroscopic quantities to default values
			rho(i,j,k,M_lim,K_lim) = rho_in;
			u(i,j,k,0,M_lim,K_lim,dims) = u_0x;
			u(i,j,k,1,M_lim,K_lim,dims) = u_0y;
#if (dims == 3)
			u(i,j,k,2,M_lim,K_lim,dims) = u_0z;
#endif

			// Set f to default values
			for (int v = 0; v < nVels; v++) {
				feq(i,j,k,v,M_lim,K_lim,nVels) = LBM_collide(i,j,k,v,M_lim,K_lim);
				f(i,j,k,v,M_lim,K_lim,nVels) = feq(i,j,k,v,M_lim,K_lim,nVels);
			}
		}


	}

}

// ***************************************************************************************************
// Routine to apply Bfl boundary condition
void GridObj::bc_applyBfl(int i, int j, int k) {

	// Declarations
	int N_lim = static_cast<int>( XInd.size() );
	int M_lim = static_cast<int>( YInd.size() );
	int K_lim = static_cast<int>( ZInd.size() );
	MarkerData* m_data;
	ObjectManager* objMan = ObjectManager::getInstance();

	// For each even outgoing direction (saves BC being applied twice otherwise)
	for (size_t v_outgoing = 0; v_outgoing < nVels; v_outgoing+=2) {

		// Get oppposite direction
		size_t v_incoming = GridUtils::getOpposite(v_outgoing);

		// Identify site where the population will be streamed to (no periodicity)
		int dest_i = i+c[0][v_outgoing];
		int dest_j = j+c[1][v_outgoing];
		int dest_k = k+c[2][v_outgoing];

		// If this site is off-grid then cannot apply the BC in preparation for streaming
		if (	(dest_i >= N_lim || dest_i < 0) ||
				(dest_j >= M_lim || dest_j < 0) ||
				(dest_k >= K_lim || dest_k < 0)
			) {
				continue;	// Move on to next direction

#ifdef BUILD_FOR_MPI
		// When using MPI, equivalent to off-grid is destination is if in periodic recv layer with 
		// periodic boundaries disabled.
		} else if (GridUtils::isOnRecvLayer(XPos[dest_i],YPos[dest_j],ZPos[dest_k]) 
			&& GridUtils::isOverlapPeriodic(dest_i,dest_j,dest_k,*this)) {

#if (!defined PERIODIC_BOUNDARIES)
			continue;	// Periodic boundaries disabled so do not try to apply BC
#endif

#endif	// BUILD_FOR_MPI

		}	// End of exclusions

		/* Not been filtered by above exclusions so try to apply boundary condition. */

		// Get marker ID
#ifdef BUILD_FOR_MPI
		// Convert locals to global for marker access
		std::vector<int> globals; GridUtils::local_to_global(i,j,k,this,globals);
		m_data = BFLBody::getMarkerData(globals[0],globals[1],globals[2],&(objMan->pBody[0]));
#else
		m_data = BFLBody::getMarkerData(i,j,k,&(objMan->pBody[0]));
#endif

		/** Apply BC in pairs  -- BC 1 **/

		// Get value of q (outgoing direction, source store)
		double q = objMan->pBody[0].Q[v_outgoing][m_data->ID];

		// Choose which implementation is appropriate //     
		
		// Half-way Bounce Back
		if (q == 0) {

			f(i,j,k,v_outgoing,M_lim,K_lim,nVels) = 
				objMan->f_prestream(i,j,k,v_incoming,M_lim,K_lim,nVels);

		// q less than 0.5 BFL bounce back (post stream)
		} else if (q < 0.5 && q > 0) {

			f(i,j,k,v_outgoing,M_lim,K_lim,nVels) = 
				(1 - 2 * q) * 
				objMan->f_prestream(i - c[0][v_outgoing],j - c[1][v_outgoing],k - c[2][v_outgoing],v_incoming,M_lim,K_lim,nVels) +
				2 * q * 
				objMan->f_prestream(i,j,k,v_incoming,M_lim,K_lim,nVels);

		// q greater than or equal to 0.5 BFL bounce back (post stream)
		} else if (q >= 0.5 && q < 1) {

			f(i,j,k,v_outgoing,M_lim,K_lim,nVels) = 
				1 / (2 * q) *
				objMan->f_prestream(i,j,k,v_incoming,M_lim,K_lim,nVels) +
				(2 * q - 1) / (2 * q) * 
				objMan->f_prestream(i,j,k,v_outgoing,M_lim,K_lim,nVels);
                          
		}		


		/** Apply BC in pairs  -- BC 2 **/

		// Get value of q (incoming direction, destination store)
		q = objMan->pBody[0].Q[v_incoming + nVels][m_data->ID];

		// Half-way Bounce Back
		if (q == 0) {

			f(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,nVels) = 
				objMan->f_prestream(dest_i,dest_j,dest_k,v_outgoing,M_lim,K_lim,nVels);

		// q less than 0.5 BFL bounce back (post stream)
		} else if (q < 0.5 && q > 0) {

			f(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,nVels) = 
				(1 - 2 * q) * 
				objMan->f_prestream(dest_i - c[0][v_incoming],dest_j - c[1][v_incoming],dest_k - c[2][v_incoming],v_outgoing,M_lim,K_lim,nVels) +
				2 * q * 
				objMan->f_prestream(dest_i,dest_j,dest_k,v_outgoing,M_lim,K_lim,nVels);

		// q greater than or equal to 0.5 BFL bounce back (post stream)
		} else if (q >= 0.5 && q < 1) {

			f(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,nVels) = 
				1 / (2 * q) *
				objMan->f_prestream(dest_i,dest_j,dest_k,v_outgoing,M_lim,K_lim,nVels) +
				(2 * q - 1) / (2 * q) * 
				objMan->f_prestream(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,nVels);
                          
		}

		delete m_data;
	}

}

// ***************************************************************************************************
// Routine to apply Non-reflecting boundary condition
void GridObj::bc_applyNrbc(int i, int j, int k) {

	// Jon's NRBC to go here

}


// ***************************************************************************************************
// Routine to apply half-way specular reflection condition for free-slip
void GridObj::bc_applySpecReflect(int label, int i, int j, int k, int N_lim, int M_lim, int K_lim) {

	int dest_x, dest_y, dest_z;

	// Check the 4 (2D) or 6 (3D) normals to get orientation of the wall
	for (size_t normal = 0; normal < (dims * 2); normal++) {

		// Identify site where the population will be streamed to //

		// Consider periodic stream (serial/non-MPI)
#if (defined PERIODIC_BOUNDARIES && !defined BUILD_FOR_MPI)

		dest_x = (i+c[0][normal] + N_lim) % N_lim;
		dest_y = (j+c[1][normal] + M_lim) % M_lim;
		dest_z = (k+c[2][normal] + K_lim) % K_lim;
#else
		dest_x = i+c[0][normal];
		dest_y = j+c[1][normal];
		dest_z = k+c[2][normal];
#endif

		// If this site is off-grid then cannot apply the BC in preparation for streaming
		if (	(dest_x >= N_lim || dest_x < 0) ||
				(dest_y >= M_lim || dest_y < 0) ||
				(dest_z >= K_lim || dest_z < 0)
			) {
				continue;	// Move on to next direction

#ifdef BUILD_FOR_MPI
		// When using MPI, equivalent to off-grid is when destination is in 
		// periodic recv layer with periodic boundaries disabled.
		} else if ( GridUtils::isOnRecvLayer(XPos[dest_x],YPos[dest_y],ZPos[dest_z]) 
			&& GridUtils::isOverlapPeriodic(dest_x,dest_y,dest_z,*this) ) {

#if (!defined PERIODIC_BOUNDARIES)
			continue;	// Periodic boundaries disabled so do not try to apply BC
#endif

			/* If periodic boundaries are enabled then this is a valid stream 
			 * so allow flow through to next section */

#endif	// BUILD_FOR_MPI

		}	// End of exclusions

		/* Not been filtered by above exclusions so try to apply boundary condition. */

		// Only apply if destination is a fluid site
		if (	LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 1 || 
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 3 ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == 4
			) {

			// This normal is valid so apply specular reflection to sites getting 
			// correct direction depending on the normal
			for (size_t v_outgoing = 0; v_outgoing < nVels; v_outgoing++) {

				f(i,j,k,v_outgoing,M_lim,K_lim,nVels) = 
					f(dest_x,dest_y,dest_z,GridUtils::dir_reflect[normal][v_outgoing],M_lim,K_lim,nVels);
			}
			break; // Do not check for any other normals
		}
	}

}


// ***************************************************************************************************
// HELPER ROUTINES
// ***************************************************************************************************
// Routine to reset the velocity at solid sites to zero
void GridObj::bc_solidSiteReset( ) {

	// Get grid sizes
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();


	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				
				// Reset solid site velocities to zero
				if (LatTyp(i,j,k,M_lim,K_lim) == 0 || LatTyp(i,j,k,M_lim,K_lim) == 10) {
					
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