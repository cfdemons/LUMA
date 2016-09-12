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
#include "../inc/BFLBody.h"
#include "../inc/ObjectManager.h"
#include <numeric>


// ***************************************************************************************************
// Method to apply suitable boundary conditions
void GridObj::LBM_boundary (int bc_type_flag) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Loop over grid, identify BC required & apply BC
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {


				/*	******************************************************************************************
					************************************* Solid Sites ****************************************
					******************************************************************************************	*/
    
				// Apply to solid sites
				if ( (LatTyp(i,j,k,M_lim,K_lim) == eSolid || LatTyp(i,j,k,M_lim,K_lim) == eRefinedSolid)
					&& (bc_type_flag == eBCAll || bc_type_flag == eBCSolidSymmetry) ) {

					// Apply half-way bounce-back
					bc_applyBounceBack(LatTyp(i,j,k,M_lim,K_lim),i,j,k);

#ifdef L_LD_OUT
					// Compute lift and drag contribution of this site
					objman->computeLiftDrag(i, j, k, this);
#endif

				} else if ( (LatTyp(i,j,k,M_lim,K_lim) == eSymmetry || LatTyp(i,j,k,M_lim,K_lim) == eRefinedSymmetry)
					&& (bc_type_flag == eBCAll || bc_type_flag == eBCSolidSymmetry)) {

					// Apply half-way specular reflection
					bc_applySpecReflect(LatTyp(i,j,k,M_lim,K_lim),i,j,k);


				/*	******************************************************************************************
					************************************* Inlet Sites ****************************************
					******************************************************************************************	*/
    
				} else if ( (LatTyp(i,j,k,M_lim,K_lim) == eInlet || LatTyp(i,j,k,M_lim,K_lim) == eRefinedInlet) 
					&& (bc_type_flag == eBCAll || bc_type_flag == eBCInlet || bc_type_flag == eBCInletOutlet) ) {

					// Choose option
#if (defined L_INLET_ON && defined L_INLET_REGULARISED)
					// Apply regularised BC
					bc_applyRegularised(LatTyp(i,j,k,M_lim,K_lim), i, j, k);
#endif
					// Do-Nothing BC is applied by default in a different way

    
				/*	******************************************************************************************
					************************************ Outlet Sites ****************************************
					******************************************************************************************	*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == eOutlet && 
					(bc_type_flag == eBCAll || bc_type_flag == eBCOutlet || bc_type_flag == eBCInletOutlet) ) {

					// !! RIGHT HAND WALL ONLY !!
#ifdef L_OUTLET_ON
					// Apply extrapolation
					bc_applyExtrapolation(LatTyp(i,j,k,M_lim,K_lim), i, j, k);
#endif



				/*	******************************************************************************************
					************************************** BFL Sites *****************************************
					******************************************************************************************	*/
    
				} else if (LatTyp(i,j,k,M_lim,K_lim) == eBFL && (bc_type_flag == eBCAll || bc_type_flag == eBCBFL) ) {

					bc_applyBfl(i,j,k);
						

				}	// End of BC type control flow
			
			}	// End of lattice site loop
		}
	}
    
}


// ***************************************************************************************************
// Routine to apply half-way bounceback boundary condition
void GridObj::bc_applyBounceBack(int label, int i, int j, int k) {

	int dest_x, dest_y, dest_z;

	// For each outgoing direction
	for (int v_outgoing = 0; v_outgoing < L_nVels; v_outgoing++) {

		// Identify site where the population will be streamed to //

		// Consider periodic stream (serial/non-MPI)
#if (defined L_PERIODIC_BOUNDARIES && !defined L_BUILD_FOR_MPI)

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

#ifdef L_BUILD_FOR_MPI
		// When using MPI, equivalent to off-grid is when destination is in 
		// periodic recv layer with periodic boundaries disabled.
		} else if ( GridUtils::isOnRecvLayer(XPos[dest_x],YPos[dest_y],ZPos[dest_z]) 
			&& GridUtils::isOverlapPeriodic(dest_x,dest_y,dest_z,*this) ) {

#if (!defined L_PERIODIC_BOUNDARIES)
			continue;	// Periodic boundaries disabled so do not try to apply BC
#endif

			/* If periodic boundaries are enabled then this is a valid stream 
			 * so allow flow through to next section */

#endif	// L_BUILD_FOR_MPI

		}	// End of exclusions

		/* Not been filtered by above exclusions so try to apply boundary condition. */

		// Only apply if destination is a fluid site or precomputed inlet
		if (	LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eFluid || 
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eTransitionToCoarser ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eTransitionToFiner ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eInlet
			) {

			// Get incoming direction
			int v_incoming = GridUtils::getOpposite(v_outgoing);
							
			// Overwriting outgoing population with expected incoming value
			f(i,j,k,v_outgoing,M_lim,K_lim,L_nVels) = f(dest_x,dest_y,dest_z,v_incoming,M_lim,K_lim,L_nVels);
							
		}
	}

}

// ***************************************************************************************************
// Routine to apply Extrapolation outlet boundary condition
void GridObj::bc_applyExtrapolation(int label, int i, int j, int k) {

	/* When using MPI, an outlet may appear at an edge of the grid which 
	 * is not the right-hand edge. As a remedy to this we first check
	 * whether the extrapolation has access to the values. If it does not
	 * the extrapolation is not performed. Although this will result in 
	 * the boundary condition not being applied to the outlet on one of 
	 * the MPI ranks, this outlet will be overwritten during MPI comms.
	 * Furthermore, the erroneous values will likely stream to an inlet 
	 * or other boundary site where the incoming values are assumed 
	 * unknown and are overwritten by the boundary condition at this site.
	 * Therefore, there is no impact of not updating these sites.
	 */

	if (i - 1 >= 0 && i - 2 >= 0)
	{
		for (size_t v = 0; v < L_nVels; v++) {

			// If right-hand wall incoming population then copy value of f from the upstream site
			if (c[0][v] == -1) f(i,j,k,v,M_lim,K_lim,L_nVels) = f(i-1,j,k,v,M_lim,K_lim,L_nVels);

		}
	}

}


// ***************************************************************************************************
// Routine to apply Regularised boundary conditions
void GridObj::bc_applyRegularised(int label, int i, int j, int k) {

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
	for (n = 0; n < L_nVels; n++) {
		ftmp.push_back(f(i,j,k,n,M_lim,K_lim,L_nVels));
	}	

	// Loop to find orientation of wall (do not include rest direction)
	for (size_t v_outgoing = 0; v_outgoing < L_nVels - 1; v_outgoing++) {

		// Identify site in this direction //
		dest_x = i+c[0][v_outgoing];
		dest_y = j+c[1][v_outgoing];
		dest_z = k+c[2][v_outgoing];

		// If this site is off-grid then cannot determine orientation of wall
		if (GridUtils::isOffGrid(dest_x,dest_y,dest_z,*this)) {

			continue;

		}

		// Only apply if destination is a fluid site or a solid site (i.e. an inlet site which feeds the domain)
		if (	LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eFluid || 
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eTransitionToCoarser ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eTransitionToFiner ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eSolid ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eRefinedSolid
			) {

			// Declarations
			double rho_wall;
			double f_plus = 0, f_zero = 0;
			std::vector<double> x_plus, x_minus, y_plus, y_minus, z_plus, z_minus;
			double u_normal;
			int normal_dir;

			// First 4 (2D) or 6 (3D) are normals so apply straight wall version of density
			if (v_outgoing < 2 * L_dims) {

				/* According to Latt & Chopard 2008 and the cited thesis by Latt 2007 we define the regularised
				 * boundary condition as folows:
				 *
				 * 1) Apply off-equilibrium bounceback to the unknown populations.
				 * 2) Compute off-equilibrium stress components L_PI^neq_ab = sum( c_ia c_ib f^neq_i ).
				 * 3) Substitute off-equilibrium definitions.
				 * 4) Compute regularised off-equilibrium part from f^neq = (w_i / (2*cs^4)) Q_iab * L_PI^neq_ab
				 * where Q_iab = dot(c_i, c_i) - cs^2 delta_ab.
				 * 5) Finally replace all populations on the inlet boundary node as f_i = f^eq_i + f^neq_i
				 * 
				 */				

				// Find nature of normal and assign normal velocity
				for (n = 0; n < L_dims; n++) {

					// Is it normal x, y or z vector?
					if (abs(c[n][v_outgoing]) == 1) {
						normal_dir = n;

						// Get appropriate normal velocity (+ve u_normal = incoming flow)
						switch (n) {
						case 0:
							u_normal = c[normal_dir][v_outgoing] * L_u_0x;
							break;
						case 1:
							u_normal = c[normal_dir][v_outgoing] * L_u_0y;
							break;
						case 2:
							u_normal = c[normal_dir][v_outgoing] * L_u_0z;
							break;
						}

						break;
					}
				}


				// Loop through direction vectors and compute the relevant momenta
				// from the known (outgoing) populations
				for (n = 0; n < L_nVels; n++) {

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
			u(i,j,k,0,M_lim,K_lim,L_dims) = L_u_0x;
			u(i,j,k,1,M_lim,K_lim,L_dims) = L_u_0y;
#if (L_dims == 3)
			u(i,j,k,2,M_lim,K_lim,L_dims) = L_u_0z;
#endif

			// Update feq to match the desired density and velocity and
			// apply off-equilibrium bounce-back to unknown populations
			for (n = 0; n < L_nVels; n++) {
				feq(i,j,k,n,M_lim,K_lim,L_nVels) = LBM_collide(i,j,k,n);


				// Different "if conditions" for corners/edges and normal walls
				if (v_outgoing < 2 * L_dims) {
					// Unknowns on normal wall share normal component
					if (c[normal_dir][n] == c[normal_dir][v_outgoing]) {
						ftmp[n] = ftmp[GridUtils::getOpposite(n)] - 
							(feq(i,j,k,GridUtils::getOpposite(n),M_lim,K_lim,L_nVels) - feq(i,j,k,n,M_lim,K_lim,L_nVels));
					}
				} else {
					// Unknowns on corners/edges are ones which do not share any component of outgoing vector
					if (c[0][n] != c[0][v_outgoing] && c[1][n] != c[1][v_outgoing]
#if (L_dims == 3)
					&& c[2][n] != c[2][v_outgoing]
#endif
					) {
						ftmp[n] = ftmp[GridUtils::getOpposite(n)] - 
							(feq(i,j,k,GridUtils::getOpposite(n),M_lim,K_lim,L_nVels) - feq(i,j,k,n,M_lim,K_lim,L_nVels));
					}
				}
			}

			// Get off-equilibrium stress components
			double Sxx = 0, Syy = 0, Sxy = 0;	// 2D & 3D
			double Szz = 0, Sxz = 0, Syz = 0;	// Just 3D
			for (n = 0; n < L_nVels; n++) {
				Sxx += c[0][n] * c[0][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,L_nVels));
				Syy += c[1][n] * c[1][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,L_nVels));
				Szz += c[2][n] * c[2][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,L_nVels));
				Sxy += c[0][n] * c[1][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,L_nVels));
				Sxz += c[0][n] * c[2][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,L_nVels));
				Syz += c[1][n] * c[2][n] * (ftmp[n] - feq(i,j,k,n,M_lim,K_lim,L_nVels));
			}

			// Compute regularised non-equilibrium components and add to feq to get new populations
			for (int v = 0; v < L_nVels; v++) {

				f(i,j,k,v,M_lim,K_lim,L_nVels) = 

					LBM_collide(i,j,k,v) +
						
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
		if (v_outgoing == L_nVels) {

			// Set macroscopic quantities to default values
			rho(i,j,k,M_lim,K_lim) = L_rho_in;
			u(i,j,k,0,M_lim,K_lim,L_dims) = L_u_0x;
			u(i,j,k,1,M_lim,K_lim,L_dims) = L_u_0y;
#if (L_dims == 3)
			u(i,j,k,2,M_lim,K_lim,L_dims) = L_u_0z;
#endif

			// Set f to default values
			for (int v = 0; v < L_nVels; v++) {
				feq(i,j,k,v,M_lim,K_lim,L_nVels) = LBM_collide(i,j,k,v);
				f(i,j,k,v,M_lim,K_lim,L_nVels) = feq(i,j,k,v,M_lim,K_lim,L_nVels);
			}
		}


	}

}

// ***************************************************************************************************
// Routine to apply Bfl boundary condition
void GridObj::bc_applyBfl(int i, int j, int k) {

	// Declarations
	MarkerData* m_data;
	ObjectManager* objMan = ObjectManager::getInstance();

	// For each even outgoing direction (saves BC being applied twice otherwise)
	for (int v_outgoing = 0; v_outgoing < L_nVels; v_outgoing+=2) {

		// Get oppposite direction
		int v_incoming = GridUtils::getOpposite(v_outgoing);

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

#ifdef L_BUILD_FOR_MPI
		// When using MPI, equivalent to off-grid is destination is if in periodic recv layer with 
		// periodic boundaries disabled.
		} else if (GridUtils::isOnRecvLayer(XPos[dest_i],YPos[dest_j],ZPos[dest_k]) 
			&& GridUtils::isOverlapPeriodic(dest_i,dest_j,dest_k,*this)) {

#if (!defined L_PERIODIC_BOUNDARIES)
			continue;	// Periodic boundaries disabled so do not try to apply BC
#endif

#endif	// L_BUILD_FOR_MPI

		}	// End of exclusions

		/* Not been filtered by above exclusions so try to apply boundary condition. */

		// Get marker ID
#ifdef L_BUILD_FOR_MPI
		// Convert locals to global for marker access
		std::vector<int> globals; GridUtils::local_to_global(i,j,k,this,globals);
		m_data = objMan->pBody[0].getMarkerData(globals[0],globals[1],globals[2]);
#else
		m_data = objMan->pBody[0].getMarkerData(i, j, k);
#endif

		/** Apply BC in pairs  -- BC 1 **/

		// Get value of q (outgoing direction, source store)
		double q = objMan->pBody[0].Q[v_outgoing][m_data->ID];

		// Choose which implementation is appropriate //     
		
		// Half-way Bounce Back
		if (q == 0) {

			f(i,j,k,v_outgoing,M_lim,K_lim,L_nVels) = 
				objMan->f_prestream(i,j,k,v_incoming,M_lim,K_lim,L_nVels);

		// q less than 0.5 BFL bounce back (post stream)
		} else if (q < 0.5 && q > 0) {

			f(i,j,k,v_outgoing,M_lim,K_lim,L_nVels) = 
				(1 - 2 * q) * 
				objMan->f_prestream(i - c[0][v_outgoing],j - c[1][v_outgoing],k - c[2][v_outgoing],v_incoming,M_lim,K_lim,L_nVels) +
				2 * q * 
				objMan->f_prestream(i,j,k,v_incoming,M_lim,K_lim,L_nVels);

		// q greater than or equal to 0.5 BFL bounce back (post stream)
		} else if (q >= 0.5 && q < 1) {

			f(i,j,k,v_outgoing,M_lim,K_lim,L_nVels) = 
				1 / (2 * q) *
				objMan->f_prestream(i,j,k,v_incoming,M_lim,K_lim,L_nVels) +
				(2 * q - 1) / (2 * q) * 
				objMan->f_prestream(i,j,k,v_outgoing,M_lim,K_lim,L_nVels);
                          
		}		


		/** Apply BC in pairs  -- BC 2 **/

		// Get value of q (incoming direction, destination store)
		q = objMan->pBody[0].Q[v_incoming + L_nVels][m_data->ID];

		// Half-way Bounce Back
		if (q == 0) {

			f(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,L_nVels) = 
				objMan->f_prestream(dest_i,dest_j,dest_k,v_outgoing,M_lim,K_lim,L_nVels);

		// q less than 0.5 BFL bounce back (post stream)
		} else if (q < 0.5 && q > 0) {

			f(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,L_nVels) = 
				(1 - 2 * q) * 
				objMan->f_prestream(dest_i - c[0][v_incoming],dest_j - c[1][v_incoming],dest_k - c[2][v_incoming],v_outgoing,M_lim,K_lim,L_nVels) +
				2 * q * 
				objMan->f_prestream(dest_i,dest_j,dest_k,v_outgoing,M_lim,K_lim,L_nVels);

		// q greater than or equal to 0.5 BFL bounce back (post stream)
		} else if (q >= 0.5 && q < 1) {

			f(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,L_nVels) = 
				1 / (2 * q) *
				objMan->f_prestream(dest_i,dest_j,dest_k,v_outgoing,M_lim,K_lim,L_nVels) +
				(2 * q - 1) / (2 * q) * 
				objMan->f_prestream(dest_i,dest_j,dest_k,v_incoming,M_lim,K_lim,L_nVels);
                          
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
void GridObj::bc_applySpecReflect(int label, int i, int j, int k) {

	int dest_x, dest_y, dest_z;

	// Check the 4 (2D) or 6 (3D) normals to get orientation of the wall
	for (size_t normal = 0; normal < (L_dims * 2); normal++) {

		// Identify site where the population will be streamed to //

		// Consider periodic stream (serial/non-MPI)
#if (defined L_PERIODIC_BOUNDARIES && !defined L_BUILD_FOR_MPI)

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

#ifdef L_BUILD_FOR_MPI
		// When using MPI, equivalent to off-grid is when destination is in 
		// periodic recv layer with periodic boundaries disabled.
		} else if ( GridUtils::isOnRecvLayer(XPos[dest_x],YPos[dest_y],ZPos[dest_z]) 
			&& GridUtils::isOverlapPeriodic(dest_x,dest_y,dest_z,*this) ) {

#if (!defined L_PERIODIC_BOUNDARIES)
			continue;	// Periodic boundaries disabled so do not try to apply BC
#endif

			/* If periodic boundaries are enabled then this is a valid stream 
			 * so allow flow through to next section */

#endif	// L_BUILD_FOR_MPI

		}	// End of exclusions

		/* Not been filtered by above exclusions so try to apply boundary condition. */

		// Only apply if destination is a fluid site
		if (	LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eFluid || 
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eTransitionToCoarser ||
				LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eTransitionToFiner
			) {

			// This normal is valid so apply specular reflection to sites getting 
			// correct direction depending on the normal
			for (size_t v_outgoing = 0; v_outgoing < L_nVels; v_outgoing++) {

				f(i,j,k,v_outgoing,M_lim,K_lim,L_nVels) = 
					f(dest_x,dest_y,dest_z,GridUtils::dir_reflect[normal][v_outgoing],M_lim,K_lim,L_nVels);
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

	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				
				// Reset solid site velocities to zero
				if (LatTyp(i,j,k,M_lim,K_lim) == eSolid || LatTyp(i,j,k,M_lim,K_lim) == eRefinedSolid) {
					
					u(i,j,k,0,M_lim,K_lim,L_dims) = 0.0;
					u(i,j,k,1,M_lim,K_lim,L_dims) = 0.0;
#if (L_dims == 3)
					u(i,j,k,2,M_lim,K_lim,L_dims) = 0.0;
#endif
				}

			}
		}
	}

}
// ***************************************************************************************************
// ***************************************************************************************************