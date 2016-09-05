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

/* This file holds all the code for the core LBM operations including collision,
streaming and macroscopic calulcation.
*/

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/definitions.h"
#include "../inc/globalvars.h"
#include "../inc/IVector.h"
#include "../inc/ObjectManager.h"
#include "../inc/MpiManager.h"

using namespace std;

// ***************************************************************************************************
// LBM multi-grid kernel applicable for both single and multi-grid (IBM assumed on coarse grid for now)
// IBM_flag dictates whether we need the predictor step or not
void GridObj::LBM_multi ( bool IBM_flag ) {

	// Start the clock to time the kernel
	clock_t secs, t_start = clock();


	///////////////////////////////
	// IBM pre-kernel processing //
	///////////////////////////////


	// Copy distributions prior to IBM predictive step
#ifdef L_IBM_ON
	// Local stores used to hold info prior to IBM predictive step
	IVector<double> f_ibm_initial, u_ibm_initial, rho_ibm_initial;

	// If IBM on and predictive loop flag true then store initial data and reset forces
	if (level == 0 && IBM_flag == true) { // Limit to level 0 immersed body for now

		*GridUtils::logfile << "Prediction step..." << std::endl;

		// Store lattice data
		f_ibm_initial = f;
		u_ibm_initial = u;
		rho_ibm_initial = rho;

		// Reset lattice and Cartesian force vectors at each site
		LBM_forcegrid(true);
	}
#else

	// If not using IBM then just reset the force vectors
	LBM_forcegrid(true);

#endif



	////////////////
	// LBM kernel //
	////////////////


	// Loop twice on refined levels as refinement ratio per level is 2
	int count = 1;
	do {

		// Apply boundary conditions (regularised must be applied before collision)
#if (defined L_INLET_ON && defined L_INLET_REGULARISED)
		LBM_boundary(2);
#endif

#ifdef L_MEGA_DEBUG
		/*DEBUG*/ io_lite((t+1)*100 + 0,"AFTER INLET BC");
#endif

		// Force lattice directions using current Cartesian force vector (adding gravity if necessary)
		LBM_forcegrid(false);

		// Collision on Lr
		LBM_collide();

#ifdef L_MEGA_DEBUG
		/*DEBUG*/ io_lite((t+1)*100 + 1,"AFTER COLLIDE");
#endif

		////////////////////
		// Refined levels //
		////////////////////

		// Check if lower level expected
		if (L_NumLev > level) {

			for (int reg = 0; reg < static_cast<int>(subGrid.size()); reg++) {

				// Explode
				LBM_explode(reg);

				// Call same routine for lower level
				subGrid[reg].LBM_multi(IBM_flag);

			}

			// Apply boundary conditions
#if (defined L_SOLID_BLOCK_ON || defined L_WALLS_ON || defined L_SOLID_FROM_FILE)
			LBM_boundary(1);	// Bounce-back (walls and solids)
#endif

#ifdef L_MEGA_DEBUG
			/*DEBUG*/ io_lite((t+1)*100 + 2,"AFTER SOLID BC");
#endif

#ifdef L_BFL_ON
			// Store the f values pre stream for BFL
			ObjectManager::getInstance()->f_prestream = f;
#endif

			// Stream
			LBM_stream();

#ifdef L_MEGA_DEBUG
			/*DEBUG*/ io_lite((t+1)*100 + 3,"AFTER STREAM");
#endif

			// Apply boundary conditions
#ifdef L_BFL_ON
			LBM_boundary(5);	// BFL boundary conditions
#endif
#ifdef L_MEGA_DEBUG
			/*DEBUG*/ io_lite((t+1)*100 + 4,"AFTER BFL");
#endif

			for (int reg = 0; reg < static_cast<int>(subGrid.size()); reg++) {

				// Coalesce
				LBM_coalesce(reg);

			}

#ifdef L_MEGA_DEBUG
			/*DEBUG*/ io_lite((t+1)*100 + 5,"AFTER COALESCE"); // Do not change this tag!
#endif


			///////////////////////
			// No refined levels //
			///////////////////////

		} else {

			// Apply boundary conditions
#if (defined L_SOLID_BLOCK_ON || defined L_WALLS_ON || defined L_SOLID_FROM_FILE)
			LBM_boundary(1);	// Bounce-back (walls and solids)
#endif

#ifdef L_MEGA_DEBUG
			/*DEBUG*/ io_lite((t+1)*100 + 2,"AFTER SOLID BC");
#endif

#ifdef L_BFL_ON
			// Store the f values pre stream for BFL
			ObjectManager::getInstance()->f_prestream = f;
#endif

			// Stream
			LBM_stream();

#ifdef L_MEGA_DEBUG
			/*DEBUG*/ io_lite((t+1)*100 + 3,"AFTER STREAM");
#endif


			// Apply boundary conditions
#ifdef L_BFL_ON
			LBM_boundary(5);	// BFL boundary conditions
#endif
#ifdef L_MEGA_DEBUG
			/*DEBUG*/ io_lite((t+1)*100 + 4,"AFTER BFL");
#endif

		}


		//////////////
		// Continue //
		//////////////

		// Apply boundary conditions
#ifdef L_OUTLET_ON
		LBM_boundary(3);	// Outlet
#endif

#ifdef L_MEGA_DEBUG
		/*DEBUG*/ io_lite((t+1)*100 + 6,"AFTER OUTLET BC");
#endif

		// Update macroscopic quantities (including time-averaged quantities)
		LBM_macro();

#ifdef L_MEGA_DEBUG
		/*DEBUG*/ io_lite((t+1)*100 + 7,"AFTER MACRO");
#endif

		// Check if on L0 and if so drop out as only need to loop once on coarsest level
		if (level == 0) {
			// Increment time step counter and break
			t++;
			break;
		}

		// Increment counters
		t++; count++;

	} while (count < 3);




	////////////////////////////////
	// IBM post-kernel processing //
	////////////////////////////////

	// Execute IBM procedure using newly computed predicted data
#ifdef L_IBM_ON
	if (level == 0 && IBM_flag == true) {

		// Reset force vectors on grid in preparation for spreading step
		LBM_forcegrid(true);

		// Calculate and apply IBM forcing to fluid
		ObjectManager::getInstance()->ibm_apply(*this);

		// Restore data to start of time step
		f = f_ibm_initial;
		u = u_ibm_initial;
		rho = rho_ibm_initial;

		// Relaunch kernel with IBM flag set to false (corrector step)
		// Corrector step does not reset force vectors but uses newly computed vector instead.
		*GridUtils::logfile << "Correction step..." << std::endl;
		t--;                // Predictor-corrector results in double time step (need to reset back 1)
		LBM_multi(false);

		// Move the body if necessary
		ObjectManager::getInstance()->ibm_move_bodies(*this);

	}
#endif

	// Get time of loop (includes sub-loops)
	secs = clock() - t_start;

	// Update average timestep time on this grid
	timeav_timestep *= (t-1);
	timeav_timestep += ((double)secs)/CLOCKS_PER_SEC;
	timeav_timestep /= t;

	if (t % L_out_every == 0) {
		// Performance data to logfile
		*GridUtils::logfile << "Grid " << level << ": Time stepping taking an average of " << timeav_timestep*1000 << "ms" << std::endl;
	}


	///////////////////////
	// MPI communication //
	///////////////////////

#ifdef L_BUILD_FOR_MPI
	/* Do MPI communication on this grid level before returning. */

	// Launch communication on this grid by passing its level and region number
	MpiManager::getInstance()->mpi_communicate(level, region_number);

#endif


}

// ***************************************************************************************************
// Takes Cartesian force vector and populates forces for each lattice direction.
// If reset_flag is true, resets the force vectors.
void GridObj::LBM_forcegrid(bool reset_flag) {

	/* This routine computes the forces applied along each direction on the lattice
	from Guo's 2002 scheme. The basic LBM must be modified in two ways: 1) the forces
	are added to the populations produced by the LBGK collision in the collision
	routine; 2) dt/2 * F is added to the momentum in the macroscopic calculation. This
	has been done already in the other routines.

	The forces along each direction are computed according to:

	F_i = (1 - 1/2*tau) * (w_i / cs^2) * ( (c_i - v) + (1/cs^2) * (c_i . v) * c_i ) . F

	where
	F_i = force applied to lattice in i-th direction
	tau = relaxation time
	w_i = weight in i-th direction
	cs = lattice sound speed
	c_i = vector of lattice speeds for i-th direction
	v = macroscopic velocity vector
	F = Cartesian force vector

	In the following, we assume
	lambda_i = (1 - 1/2*tau) * (w_i / cs^2)
	beta_i = (1/cs^2) * (c_i . v)

	The above in shorthand becomes summing over d dimensions:

	F_i = lambda_i * sum( F_d * (c_d_i (1+beta_i) - u_d )

	*/

	if (reset_flag) {

		// Reset lattice force vectors on every grid site
		std::fill(force_i.begin(), force_i.end(), 0.0);

		// Reset Cartesian force vector on every grid site
		std::fill(force_xyz.begin(), force_xyz.end(), 0.0);

	} else {
	// Else, compute forces

		// Get grid sizes
		size_t N_lim = XPos.size();
		size_t M_lim = YPos.size();
		size_t K_lim = ZPos.size();

		// Declarations
		double lambda_v;

		// Loop over grid and overwrite forces for each direction
		for (size_t i = 0; i < N_lim; i++) {
			for (size_t j = 0; j < M_lim; j++) {
				for (size_t k = 0; k < K_lim; k++) {

#ifdef GRAVITY_ON
					// Add gravity to any IBM forces currently stored
					force_xyz(i,j,k,L_grav_direction,M_lim,K_lim,L_dims) += rho(i,j,k,M_lim,K_lim) * L_grav_force;
#endif

					// Now compute force_i components from Cartesian force vector
					for (size_t v = 0; v < L_nVels; v++) {

						// Only apply to non-solid sites
						if (LatTyp(i,j,k,M_lim,K_lim) != eSolid) {

							// Reset beta_v
							double beta_v = 0.0;

							// Compute the lattice forces based on Guo's forcing scheme
							lambda_v = (1 - 0.5 * omega) * ( w[v] / (cs*cs) );

							// Dot product (sum over d dimensions)
							for (int d = 0; d < L_dims; d++) {
								beta_v +=  (c[d][v] * u(i,j,k,d,M_lim,K_lim,L_dims));
							}
							beta_v = beta_v * (1/(cs*cs));

							// Compute force using shorthand sum described above
							for (int d = 0; d < L_dims; d++) {
								force_i(i,j,k,v,M_lim,K_lim,L_nVels) += force_xyz(i,j,k,d,M_lim,K_lim,L_dims) *
									(c[d][v] * (1 + beta_v) - u(i,j,k,d,M_lim,K_lim,L_dims));
							}

							// Multiply by lambda_v
							force_i(i,j,k,v,M_lim,K_lim,L_nVels) = force_i(i,j,k,v,M_lim,K_lim,L_nVels) * lambda_v;

						}

					}

				}
			}
		}


#ifdef L_IBM_DEBUG
		// DEBUG -- write out force components
		std::ofstream testout;
		testout.open(GridUtils::path_str + "/force_i_LB.out", std::ios::app);
		testout << "\nNEW TIME STEP" << std::endl;
		for (size_t j = 1; j < M_lim - 1; j++) {
			for (size_t i = 0; i < N_lim; i++) {
				for (size_t v = 0; v < L_nVels; v++) {
					testout << force_i(i,j,0,v,M_lim,K_lim,L_nVels) << "\t";
				}
				testout << std::endl;
			}
			testout << std::endl;
		}
		testout.close();
#endif
	}

}


// ***************************************************************************************************
// Collision operator
void GridObj::LBM_collide( ) {

	/*
	Loop through the lattice points to compute the new distribution functions.
	Equilibrium based on:
	       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4)
	*/

	// Declarations and Grid size
	int N_lim = static_cast<int>(XPos.size());
	int M_lim = static_cast<int>(YPos.size());
	int K_lim = static_cast<int>(ZPos.size());

	// Create temporary lattice to prevent overwriting useful populations and initialise with same values as
	// pre-collision f grid. Initialise with current f values.
	IVector<double> f_new( f );


	// Loop over all lattice sites
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// Ignore refined sites and TL sites (based on Rohde refinement)
				if (LatTyp(i,j,k,M_lim,K_lim) == eRefined || LatTyp(i,j,k,M_lim,K_lim) == eTransitionToCoarser) {
					// Do nothing as taken care of on lower level grid

				} else {

					

#ifdef L_USE_KBC_COLLISION

					// KBC collision
					LBM_kbcCollide(i, j, k, M_lim, K_lim, f_new);

#else

					// Loop over directions and perform collision
					for (int v = 0; v < L_nVels; v++) {

						// Get feq value by calling overload of collision function
						feq(i,j,k,v,M_lim,K_lim,L_nVels) = LBM_collide(i, j, k, v, M_lim, K_lim);

						// LBGK collision
						f_new(i,j,k,v,M_lim,K_lim,L_nVels) = 
							f(i,j,k,v,M_lim,K_lim,L_nVels) - 
							omega * ( 
								f(i,j,k,v,M_lim,K_lim,L_nVels) - feq(i,j,k,v,M_lim,K_lim,L_nVels) 
									) +
							force_i(i,j,k,v,M_lim,K_lim,L_nVels);

					}

#endif

				}

			}
		}
	}

	// Update f from fnew
	f = f_new;

}

// ***************************************************************************************************
// Overload of collision function to allow calculation of feq only for initialisation
double GridObj::LBM_collide( int i, int j, int k, int v, int M_lim, int K_lim ) {

	/* LBGK equilibrium function is represented as:
		feq_i = rho * w_i * ( 1 + u_a c_ia / cs^2 + Q_iab u_a u_b / 2*cs^4 )
	where
		Q_iab = c_ia c_ib - cs^2 * delta_ab
	and
		delta_ab is the Kronecker delta.
	*/

	// Declare single feq value and intermediate values A and B
	double feq, A, B;

	// Compute the parts of the expansion for feq (we now have a dot product routine so could simplify this code)

#if (L_dims == 3)
		// Compute c_ia * u_a which is actually the dot product of c and u
		A = (c[0][v] * u(i,j,k,0,M_lim,K_lim,L_dims)) + (c[1][v] * u(i,j,k,1,M_lim,K_lim,L_dims)) + (c[2][v] * u(i,j,k,2,M_lim,K_lim,L_dims));

		/*
		Compute second term in the expansion
		Q_iab u_a u_b =
		(c_x^2 - cs^2)u_x^2 + (c_y^2 - cs^2)u_y^2 + (c_z^2 - cs^2)u_z^2
		+ 2c_x c_y u_x u_y + 2c_x c_z u_x u_z + 2c_y c_z u_y u_z
		*/

		B =	((c[0][v]*c[0][v]) - (cs*cs)) * (u(i,j,k,0,M_lim,K_lim,L_dims)*u(i,j,k,0,M_lim,K_lim,L_dims)) +
			((c[1][v]*c[1][v]) - (cs*cs)) * (u(i,j,k,1,M_lim,K_lim,L_dims)*u(i,j,k,1,M_lim,K_lim,L_dims)) +
			((c[2][v]*c[2][v]) - (cs*cs)) * (u(i,j,k,2,M_lim,K_lim,L_dims)*u(i,j,k,2,M_lim,K_lim,L_dims)) +
			2 * c[0][v]*c[1][v] * u(i,j,k,0,M_lim,K_lim,L_dims) * u(i,j,k,1,M_lim,K_lim,L_dims) +
			2 * c[0][v]*c[2][v] * u(i,j,k,0,M_lim,K_lim,L_dims) * u(i,j,k,2,M_lim,K_lim,L_dims) +
			2 * c[1][v]*c[2][v] * u(i,j,k,1,M_lim,K_lim,L_dims) * u(i,j,k,2,M_lim,K_lim,L_dims);
#else
		// 2D versions of the above
		A = (c[0][v] * u(i,j,k,0,M_lim,K_lim,L_dims)) + (c[1][v] * u(i,j,k,1,M_lim,K_lim,L_dims));

		B =	((c[0][v]*c[0][v]) - (cs*cs)) * (u(i,j,k,0,M_lim,K_lim,L_dims)*u(i,j,k,0,M_lim,K_lim,L_dims)) +
			((c[1][v]*c[1][v]) - (cs*cs)) * (u(i,j,k,1,M_lim,K_lim,L_dims)*u(i,j,k,1,M_lim,K_lim,L_dims)) +
			2 * c[0][v]*c[1][v] * u(i,j,k,0,M_lim,K_lim,L_dims) * u(i,j,k,1,M_lim,K_lim,L_dims);
#endif


	// Compute f^eq
	feq = rho(i,j,k,M_lim,K_lim) * w[v] * ( 1 + (A / (cs*cs)) + (B / (2*(cs*cs*cs*cs))) );

	return feq;

}

// ***************************************************************************************************
// KBC collision operator
void GridObj::LBM_kbcCollide( int i, int j, int k, int M_lim, int K_lim, IVector<double>& f_new ) {
	
	// Declarations
	double ds[L_nVels], dh[L_nVels], gamma;

	// Compute required moments and equilibrium moments
#if (L_dims == 3)
		
	// Most moments are required in 3D for the KBC-N4 model

	// Stress (second order)
	double M200 = 0.0, M200eq = 0.0;
	double M020 = 0.0, M020eq = 0.0;
	double M002 = 0.0, M002eq = 0.0;
	// Pis (second order)
	double M110 = 0.0, M110eq = 0.0;
	double M101 = 0.0, M101eq = 0.0;
	double M011 = 0.0, M011eq = 0.0;
	// Qs (third order)
	double M111 = 0.0, M111eq = 0.0;
	double M102 = 0.0, M102eq = 0.0;
	double M210 = 0.0, M210eq = 0.0;
	double M021 = 0.0, M021eq = 0.0;
	double M201 = 0.0, M201eq = 0.0;
	double M120 = 0.0, M120eq = 0.0;
	double M012 = 0.0, M012eq = 0.0;

	for (int v = 0; v < L_nVels; v++) {
		
		// Update feq
		feq(i,j,k,v,M_lim,K_lim,L_nVels) = LBM_collide(i, j, k, v, M_lim, K_lim);

		// These are actually rho * MXXX but no point in dividing to multiply later
		M200 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v]);
		M020 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[1][v]);
		M002 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[2][v] * c[2][v]);
		M110 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v]);
		M101 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[2][v]);
		M011 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[2][v]);
		M111 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v] * c[2][v]);
		M102 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[2][v] * c[2][v]);
		M210 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v] * c[1][v]);
		M021 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[1][v] * c[2][v]);
		M201 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v] * c[2][v]);
		M120 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v] * c[1][v]);
		M012 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[2][v] * c[2][v]);

		M200eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v]);
		M020eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[1][v]);
		M002eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[2][v] * c[2][v]);
		M110eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v]);
		M101eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[2][v]);
		M011eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[2][v]);
		M111eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v] * c[2][v]);
		M102eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[2][v] * c[2][v]);
		M210eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v] * c[1][v]);
		M021eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[1][v] * c[2][v]);
		M201eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v] * c[2][v]);
		M120eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v] * c[1][v]);
		M012eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[2][v] * c[2][v]);
	}

	// Compute ds
	for (int v = 0; v < L_nVels; v++) {

		// s part dictated by KBC model choice and directions
		if (c[0][v] == 0 && c[1][v] == 0 && c[2][v] == 0) {

			// First family
			ds[v] = ( -(M200 + M020 + M002) ) - 
					( -(M200eq + M020eq + M002eq) );

		} else if (c[0][v] != 0 && c[1][v] == 0 && c[2][v] == 0) {

			// Second family
			ds[v] =	( (2 * (M200 - M002) - (M020 - M002)) / 6 + (M200 + M020 + M002) / 6  - c[0][v] * 0.5 * (M120 + M102) ) - 
					( (2 * (M200eq - M002eq) - (M020eq - M002eq)) / 6 + (M200eq + M020eq + M002eq) / 6 - c[0][v] * 0.5 * (M120eq + M102eq));

		} else if (c[0][v] == 0 && c[1][v] != 0 && c[2][v] == 0) {

			// Third family
			ds[v] =	( (-(M200 - M002) + 2 * (M020 - M002)) / 6 + (M200 + M020 + M002) / 6 - c[1][v] * 0.5 * (M210 + M012) ) - 
					( (-(M200eq - M002eq) + 2 * (M020eq - M002eq)) / 6 + (M200eq + M020eq + M002eq) / 6 - c[1][v] * 0.5 * (M210eq + M012eq));

		} else if (c[0][v] == 0 && c[1][v] == 0 && c[2][v] != 0) {

			// Fourth family
			ds[v] =	( (-(M200 - M002) - (M020 - M002)) / 6 + (M200 + M020 + M002) / 6 - c[2][v] * 0.5 * (M201 + M021) ) - 
					( (-(M200eq - M002eq) + 2 * (M020eq - M002eq)) / 6 + (M200eq + M020eq + M002eq) / 6 - c[2][v] * 0.5 * (M201eq + M021eq) );

		} else if (c[0][v] != 0 && c[1][v] != 0 && c[2][v] == 0) {

			// Fifth family
			ds[v] =	( c[0][v] * c[1][v] * 0.25 * M110 + (c[1][v] * 0.25 * M210 + c[0][v] * 0.25 * M120) ) - 
					( c[0][v] * c[1][v] * 0.25 * M110eq + (c[1][v] * 0.25 * M210eq + c[0][v] * 0.25 * M120eq) );

		} else if (c[0][v] != 0 && c[1][v] == 0 && c[2][v] != 0) {

			// Sixth family
			ds[v] =	( c[0][v] * c[2][v] * 0.25 * M101 + (c[2][v] * 0.25 * M201 + c[0][v] * 0.25 * M102) ) - 
					( c[0][v] * c[2][v] * 0.25 * M101eq + (c[2][v] * 0.25 * M201eq + c[0][v] * 0.25 * M102eq) );

		} else if (c[0][v] == 0 && c[1][v] != 0 && c[2][v] != 0) {

			// Seventh family
			ds[v] =	( c[1][v] * c[2][v] * 0.25 * M011 + (c[2][v] * 0.25 * M021 + c[1][v] * 0.25 * M012) ) - 
					( c[1][v] * c[2][v] * 0.25 * M011eq + (c[2][v] * 0.25 * M021eq + c[1][v] * 0.25 * M012eq) );

		} else if (c[0][v] != 0 && c[1][v] != 0 && c[2][v] != 0) {

			// Eighth family
			ds[v] =	( c[0][v] * c[1][v] * c[2][v] * M111 / 8 ) - 
					( c[0][v] * c[1][v] * c[2][v] * M111eq / 8 );

		}


		// Compute dh
		dh[v] = f(i,j,k,v,M_lim,K_lim,L_nVels) - feq(i,j,k,v,M_lim,K_lim,L_nVels) - ds[v];

	}



#else

	// Only need to compute these 3 in 2D for KBC-D model
	double M20 = 0.0, M20eq = 0.0;
	double M02 = 0.0, M02eq = 0.0;
	double M11 = 0.0, M11eq = 0.0;

	for (int v = 0; v < L_nVels; v++) {
		
		// Update feq
		feq(i,j,k,v,M_lim,K_lim,L_nVels) = LBM_collide(i, j, k, v, M_lim, K_lim);

		// These are actually rho * MXX but no point in dividing to multiply later
		M20 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v]);
		M02 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[1][v]);
		M11 += f(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v]);

		M20eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[0][v]);
		M02eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[1][v] * c[1][v]);
		M11eq += feq(i,j,k,v,M_lim,K_lim,L_nVels) * (c[0][v] * c[1][v]);
	}

	// Compute ds
	for (int v = 0; v < L_nVels; v++) {

		// s part dictated by KBC model choice and directions
		if (c[0][v] == 0 && c[1][v] == 0) {

			// First family
			ds[v] = 0.0;

		} else if (c[0][v] != 0 && c[1][v] == 0) {

			// Second family
			ds[v] =	( 0.25 * (M20 - M02) ) - 
					( 0.25 * (M20eq - M02eq) );

		} else if (c[0][v] == 0 && c[1][v] != 0) {

			// Third family
			ds[v] =	( -0.25 * (M20 - M02) ) -
					( -0.25 * (M20eq - M02eq) );

		} else {

			// Fourth family
			ds[v] =	( 0.25 * c[0][v] * c[1][v] * M11 ) - 
					( 0.25 * c[0][v] * c[1][v] * M11eq );

		}


		// Compute dh
		dh[v] = f(i,j,k,v,M_lim,K_lim,L_nVels) - feq(i,j,k,v,M_lim,K_lim,L_nVels) - ds[v];

	}

#endif

	// Once all dh and ds have been computed, compute the products
	double top_prod = 0.0, bot_prod = 0.0;
	for (int v = 0; v < L_nVels; v++) {

		// Compute scalar products
		top_prod += ds[v] * dh[v] / feq(i,j,k,v,M_lim,K_lim,L_nVels);
		bot_prod += dh[v] * dh[v] / feq(i,j,k,v,M_lim,K_lim,L_nVels);

	}
	
	// Compute gamma
	if (bot_prod == 0.0) gamma = (2/omega);
	else gamma = (2/omega) - ( 2 - (2/omega) ) * (top_prod / bot_prod);

	// Finally perform collision
	for (int v = 0; v < L_nVels; v++) {

		// Perform collision
		f_new(i,j,k,v,M_lim,K_lim,L_nVels) = 
			f(i,j,k,v,M_lim,K_lim,L_nVels) - 
			(omega/2) * ( 2 * ds[v] + gamma * dh[v] ) +
			force_i(i,j,k,v,M_lim,K_lim,L_nVels);

	}


}


// ***************************************************************************************************
// Streaming operator
// Applies periodic BCs on level 0
void GridObj::LBM_stream( ) {

	/* This streaming operation obeys the following logic process:
	 *
	 *	1) Apply Source-based Exclusions (e.g. any fine sites does not need to be streamed)
	 *
	 * Then either:
	 *	2a) Apply Off-Grid Ops (e.g. retain incoming values if value streams off-grid or apply periodic BCs)
	 *
	 * Or
	 *	2b) Apply Destination-based Exclusions (e.g. do not stream to a do-nothing inlet)
	 *	2c) Stream value
	 *
	 * End
	 */

	// Declarations
	int N_lim = static_cast<int>(XPos.size());
	int M_lim = static_cast<int>(YPos.size());
	int K_lim = static_cast<int>(ZPos.size());
	int dest_x, dest_y, dest_z;
	int v_opp;

	// Create temporary lattice of zeros to prevent overwriting useful populations
	IVector<double> f_new( f.size(), 0.0 );	// Could just initialise to f to make the logic below simpler //


#ifdef L_DEBUG_STREAM
	/*DEBUG*/
	int count0 = 0, count1 = 0, count2 = 0, count3 = 0;
#endif

	// Stream one lattice site at a time
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				for (int v = 0; v < L_nVels; v++) {

					// Store opposite direction
					v_opp = GridUtils::getOpposite(v);

#ifdef L_DEBUG_STREAM
					/*DEBUG*/
					count0++;
#endif


					/////////////////////////////////
					// Streaming Source Exclusions //
					/////////////////////////////////

					/* This section prevents streaming operations given the type of site
					 * FROM which the streaming takes place regardless of the destination type.
					 */

					// Fine --> Any; do not stream in any direction
					if (LatTyp(i,j,k,M_lim,K_lim) == eRefined) {
						break;

					// Do-nothing-inlet --> Any; copy value to new grid (i.e. apply do-nothing inlet)
#if (defined L_INLET_ON && !defined L_INLET_REGULARISED && !defined L_INLET_NRBC)					
					} else if (LatTyp(i,j,k,M_lim,K_lim) == eInlet || LatTyp(i,j,k,M_lim,K_lim) == eRefinedInlet) {
						f_new(i,j,k,v,M_lim,K_lim,L_nVels) = f(i,j,k,v,M_lim,K_lim,L_nVels);
						// Carry on and stream
#endif
					}


					////////////////////////
					// Off-grid streaming //
					////////////////////////

					/* If destination off-grid then ask whether periodic boundaries in use
					 * and if so then only stream if Coarse --> Coarse.
					 * If not then retain the incoming value at the site as it will not receive
					 * an update from off-grid.
					 * If using MPI, periodic BCs are applied differently later.
					 */

					// Compute destination site (no periodicity)
					dest_x = i+c[0][v];
					dest_y = j+c[1][v];
					dest_z = k+c[2][v];


					// If off-grid
					if (	(dest_x >= N_lim || dest_x < 0) ||
							(dest_y >= M_lim || dest_y < 0)
#if (L_dims == 3)
							|| (dest_z >= K_lim || dest_z < 0)
#endif
						) {


#ifdef L_DEBUG_STREAM
							/*DEBUG*/
							count1++;
							*GridUtils::logfile << "Stream " << i << "," << j << "," << k << 
								" (" << XPos[i] << "," << YPos[j] << "," << ZPos[k] << ")" <<
								" to \t" << dest_x << "," << dest_y << "," << dest_z <<
								" : \toff-grid in " <<
								v << " direction. Count1 = " << count1 << ". Value is f = " 
								<< f(i,j,k,v,M_lim,K_lim,L_nVels) << std::endl;
#endif



							// Apply periodic boundary conditions (serial/non-MPI)
#if (defined L_PERIODIC_BOUNDARIES  && !defined L_BUILD_FOR_MPI)

							// Compute destination site indices using periodicity
							dest_x = (i+c[0][v] + N_lim) % N_lim;
							dest_y = (j+c[1][v] + M_lim) % M_lim;
							dest_z = (k+c[2][v] + K_lim) % K_lim;

							/* Only apply periodic BCs on coarsest level and if stream is:
							 * Coarse --> Coarse
							 * Solid --> Coarse (to allow half-way BB to be applied)
							 * Coarse --> Solid (for consistency with the above) */
							if (
								(level == 0) &&
								( 	(LatTyp(i,j,k,M_lim,K_lim) == eFluid && LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eFluid) ||
									(LatTyp(i,j,k,M_lim,K_lim) == eSolid && LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eFluid) ||
									(LatTyp(i,j,k,M_lim,K_lim) == eFluid && LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eSolid)	)
							) {
								// Stream periodically
								f_new(dest_x,dest_y,dest_z,v,M_lim,K_lim,L_nVels) = f(i,j,k,v,M_lim,K_lim,L_nVels);
								continue;	// Move on to next site
							}

#endif

							// If not using periodic boundary conditions then retain incoming value as 
							// will not receive an update from off-grid and continue
							f_new(i,j,k,v_opp,M_lim,K_lim,L_nVels) = f(i,j,k,v_opp,M_lim,K_lim,L_nVels);
							continue;



					} else {


						///////////////////////
						// On-grid streaming //
						///////////////////////


						/* If it is an on-grid stream and using MPI need some
						 * additional logic checking to make sure we prevent stream
						 * from a periodic recv layer site when periodic BCs are not in effect. */
						
#ifdef L_BUILD_FOR_MPI

						//////////////////////
						// MPI Periodic BCs //
						//////////////////////

						/* If source in recv layer and destination in sender layer 
						 * check if this is a periodic stream. If periodic BCs are 
						 * in use then allow stream. If not then do not stream from 
						 * recv layer to sender layer as this would constitute an 
						 * illegal periodic stream. */
						if (

						(

						// Condition 1: Source in a recv layer
						GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k])

						) && (

						// Condition 2: Destination on-grid (in a sender layer)
						GridUtils::isOnSenderLayer(XPos[dest_x],YPos[dest_y],ZPos[dest_z])
						

						) && (

						// Condition 3: Recv layer is linked to a periodic neighbour rank
						GridUtils::isOverlapPeriodic(i,j,k,*this)

						)

						) {


#ifdef L_DEBUG_STREAM
							/*DEBUG*/
							count2++;
							*GridUtils::logfile << "Stream " << i << "," << j << "," << k << 
								" (" << XPos[i] << "," << YPos[j] << "," << ZPos[k] << ")" <<
								" to \t" << dest_x << "," << dest_y << "," << dest_z <<
								" (" << XPos[dest_x] << "," << YPos[dest_y] << "," << ZPos[dest_z] << ")" << 
								" : \tperiodic stream " <<
									v << " direction. Count2 = " << count2 << ". Value is f = " 
									<< f(i,j,k,v,M_lim,K_lim,L_nVels) << std::endl;
#endif


						// Either apply periodic BCs or not on these sites (ones which stream from periodic recv site)

#ifdef L_PERIODIC_BOUNDARIES
							if (LatTyp(i,j,k,M_lim,K_lim) == eFluid && LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eFluid)
							{
								// Stream periodically
								f_new(dest_x,dest_y,dest_z,v,M_lim,K_lim,L_nVels) = f(i,j,k,v,M_lim,K_lim,L_nVels);
								continue;

							} else {
								// Incoming value at the destination site should be retained as no periodic BC to be applied
								f_new(dest_x,dest_y,dest_z,v,M_lim,K_lim,L_nVels) = f(dest_x,dest_y,dest_z,v,M_lim,K_lim,L_nVels);
								continue;

							}


						}
#else

							// If not using periodic BCs then incoming value at the destination site should be retained
							f_new(dest_x,dest_y,dest_z,v,M_lim,K_lim,L_nVels) = f(dest_x,dest_y,dest_z,v,M_lim,K_lim,L_nVels);
							continue;

						}

#endif	// L_PERIODIC_BOUNDARIES

#endif	// L_BUILD_FOR_MPI




						//////////////////////////////////////
						// Streaming Destination Exclusions //
						//////////////////////////////////////

						/* Filter out unwanted streaming operations by checking the
						 * source-destination pairings and retaining those values
						 * that you do not want to be overwritten by streaming.
						 * This section prevents streaming operations given the type of site
						 * TO which the streaming takes place regardless of the source type.
						 */

						// TL2lower --> TL2lower then ignore as done on lower grid stream
						if (
							(LatTyp(i,j,k,M_lim,K_lim) == eTransitionToFiner) &&
							(LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eTransitionToFiner)
						) {
							continue;

						// Any --> Do-nothing-inlet then ignore so as not to overwrite inlet site
#if (defined L_INLET_ON && !defined L_INLET_REGULARISED && !defined L_INLET_NRBC)
						} else if (LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eInlet || LatTyp(dest_x,dest_y,dest_z,M_lim,K_lim) == eRefinedInlet) {
							continue;
#endif
						}
						

#ifdef L_DEBUG_STREAM
						/*DEBUG*/
						count3++;
						*GridUtils::logfile << "Stream " << i << "," << j << "," << k << 
								" (" << XPos[i] << "," << YPos[j] << "," << ZPos[k] << ")" <<
								" to \t" << dest_x << "," << dest_y << "," << dest_z <<
								" (" << XPos[dest_x] << "," << YPos[dest_y] << "," << ZPos[dest_z] << ")" << 
								" : \ton-grid stream " <<
								v << " direction. Count3 = " << count3 << ". Value is f = " 
								<< f(i,j,k,v,M_lim,K_lim,L_nVels) << std::endl;
#endif


						// Stream population
						f_new(dest_x,dest_y,dest_z,v,M_lim,K_lim,L_nVels) = f(i,j,k,v,M_lim,K_lim,L_nVels);

					}


				}

			}
		}
	}


#ifdef L_DEBUG_STREAM
	/*DEBUG*/
	*GridUtils::logfile << "Counts were " << count0 << "," << count1 << "," << count2 << "," << count3 << std::endl;
#endif

	// Replace old grid with new grid
	f = f_new;

}


// ***************************************************************************************************
// Macroscopic quantity calculation and update of time-averaged quantities
void GridObj::LBM_macro( ) {

	// Declarations
	int N_lim = static_cast<int>(XPos.size());
	int M_lim = static_cast<int>(YPos.size());
	int K_lim = static_cast<int>(ZPos.size());
	double rho_temp = 0.0;
	double fux_temp = 0.0;
	double fuy_temp = 0.0;
	double fuz_temp = 0.0;
	double ta_temp = 0.0;


	// Loop over lattice
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				if (LatTyp(i,j,k,M_lim,K_lim) == eRefined) {

					// Refined site so set both density and velocity to zero
					rho(i,j,k,M_lim,K_lim) = 0.0;
					u(i,j,k,0,M_lim,K_lim,L_dims) = 0.0;
					u(i,j,k,1,M_lim,K_lim,L_dims) = 0.0;
#if (L_dims == 3)
					u(i,j,k,2,M_lim,K_lim,L_dims) = 0.0;
#endif

				} else if (LatTyp(i,j,k,M_lim,K_lim) == eSolid || LatTyp(i,j,k,M_lim,K_lim) == eRefinedSolid) {

					// Solid site so do not update density but set velocity to zero
					rho(i,j,k,M_lim,K_lim) = 1.0;
					u(i,j,k,0,M_lim,K_lim,L_dims) = 0.0;
					u(i,j,k,1,M_lim,K_lim,L_dims) = 0.0;
#if (L_dims == 3)
					u(i,j,k,2,M_lim,K_lim,L_dims) = 0.0;
#endif

				} else if ( LatTyp(i,j,k,M_lim,K_lim) == eSymmetry ||
					(LatTyp(i,j,k,M_lim,K_lim) == eInlet || LatTyp(i,j,k,M_lim,K_lim) == eRefinedInlet) ) {

					// Symmetry or Inlet BC which update themselves prior to collision
					continue;

				} else {

					// Any other of type of site compute both density and velocity from populations
					rho_temp = 0.0; fux_temp = 0.0; fuy_temp = 0.0; fuz_temp = 0.0;

					for (int v = 0; v < L_nVels; v++) {

						// Sum up to find mass flux
						fux_temp += (double)c[0][v] * f(i,j,k,v,M_lim,K_lim,L_nVels);
						fuy_temp += (double)c[1][v] * f(i,j,k,v,M_lim,K_lim,L_nVels);
						fuz_temp += (double)c[2][v] * f(i,j,k,v,M_lim,K_lim,L_nVels);

						// Sum up to find density
						rho_temp += f(i,j,k,v,M_lim,K_lim,L_nVels);

					}

					// Assign density
					rho(i,j,k,M_lim,K_lim) = rho_temp;

					// Add forces to momentum (rho * time step * 0.5 * force -- eqn 19 in Favier 2014)
					fux_temp += rho_temp * (1 / pow(2,level)) * 0.5 * force_xyz(i,j,k,0,M_lim,K_lim,L_dims);
					fuy_temp += rho_temp * (1 / pow(2,level)) * 0.5 * force_xyz(i,j,k,1,M_lim,K_lim,L_dims);
#if (L_dims == 3)
					fuz_temp += rho_temp * (1 / pow(2,level)) * 0.5 * force_xyz(i,j,k,2,M_lim,K_lim,L_dims);
#endif

					// Assign velocity
					u(i,j,k,0,M_lim,K_lim,L_dims) = fux_temp / rho_temp;
					u(i,j,k,1,M_lim,K_lim,L_dims) = fuy_temp / rho_temp;
#if (L_dims == 3)
					u(i,j,k,2,M_lim,K_lim,L_dims) = fuz_temp / rho_temp;
#endif

				}

				// Update time-averaged quantities thus...

				// Multiply current value by completed time steps to get sum
				ta_temp = rho_timeav(i,j,k,M_lim,K_lim) * (double)t;
				// Add new value
				ta_temp += rho(i,j,k,M_lim,K_lim);
				// Divide by completed time steps + 1 to get new average
				rho_timeav(i,j,k,M_lim,K_lim) = ta_temp / (double)(t+1);

				// Repeat for other quantities
				int pq_combo = 0;
				for (int p = 0; p < L_dims; p++) {
					ta_temp = ui_timeav(i,j,k,p,M_lim,K_lim,L_dims) * (double)t;
					ta_temp += u(i,j,k,p,M_lim,K_lim,L_dims);
					ui_timeav(i,j,k,p,M_lim,K_lim,L_dims) = ta_temp / (double)(t+1);
					// Do necessary products
					for (int q = p; q < L_dims; q++) {
						ta_temp = uiuj_timeav(i,j,k,pq_combo,M_lim,K_lim,(3*L_dims-3)) * (double)t;
						ta_temp += ( u(i,j,k,p,M_lim,K_lim,L_dims) * u(i,j,k,q,M_lim,K_lim,L_dims) );
						uiuj_timeav(i,j,k,pq_combo,M_lim,K_lim,(3*L_dims-3)) = ta_temp / (double)(t+1);
						pq_combo++;
					}
				}



			}
		}
	}

}


// ***************************************************************************************************
// Overload of macroscopic quantity calculation to allow it to be applied to a single site as used by
// the MPI unpacking routine to update the values for the next collision step. This routine does not
// update the time-averaged quantities.
void GridObj::LBM_macro( int i, int j, int k ) {

	// Declarations
	int N_lim = static_cast<int>(XPos.size());
	int M_lim = static_cast<int>(YPos.size());
	int K_lim = static_cast<int>(ZPos.size());
	double rho_temp = 0.0;
	double fux_temp = 0.0;
	double fuy_temp = 0.0;
	double fuz_temp = 0.0;


	if (LatTyp(i,j,k,M_lim,K_lim) == eRefined) {

		// Refined site so set both density and velocity to zero
		rho(i,j,k,M_lim,K_lim) = 0.0;
		u(i,j,k,0,M_lim,K_lim,L_dims) = 0.0;
		u(i,j,k,1,M_lim,K_lim,L_dims) = 0.0;
#if (L_dims == 3)
		u(i,j,k,2,M_lim,K_lim,L_dims) = 0.0;
#endif

	} else if (LatTyp(i,j,k,M_lim,K_lim) == eSolid || LatTyp(i,j,k,M_lim,K_lim) == eRefinedSolid) {

		// Solid site so do not update density but set velocity to zero
		rho(i,j,k,M_lim,K_lim) = 1.0;
		u(i,j,k,0,M_lim,K_lim,L_dims) = 0.0;
		u(i,j,k,1,M_lim,K_lim,L_dims) = 0.0;
#if (L_dims == 3)
		u(i,j,k,2,M_lim,K_lim,L_dims) = 0.0;
#endif

	} else if ( LatTyp(i,j,k,M_lim,K_lim) == eSymmetry ||
		(LatTyp(i,j,k,M_lim,K_lim) == eInlet || LatTyp(i,j,k,M_lim,K_lim) == eRefinedInlet) ) {

		// Symmetry or Inlet BC which update themselves prior to collision

	} else {

		// Any other of type of site compute both density and velocity from populations
		rho_temp = 0.0; fux_temp = 0.0; fuy_temp = 0.0; fuz_temp = 0.0;

		for (int v = 0; v < L_nVels; v++) {

			// Sum up to find mass flux
			fux_temp += (double)c[0][v] * f(i,j,k,v,M_lim,K_lim,L_nVels);
			fuy_temp += (double)c[1][v] * f(i,j,k,v,M_lim,K_lim,L_nVels);
			fuz_temp += (double)c[2][v] * f(i,j,k,v,M_lim,K_lim,L_nVels);

			// Sum up to find density
			rho_temp += f(i,j,k,v,M_lim,K_lim,L_nVels);

		}

		// Assign density
		rho(i,j,k,M_lim,K_lim) = rho_temp;

		// Add forces to momentum (rho * time step * 0.5 * force -- eqn 19 in Favier 2014)
		fux_temp += rho_temp * (1 / pow(2,level)) * 0.5 * force_xyz(i,j,k,0,M_lim,K_lim,L_dims);
		fuy_temp += rho_temp * (1 / pow(2,level)) * 0.5 * force_xyz(i,j,k,1,M_lim,K_lim,L_dims);
#if (L_dims == 3)
		fuz_temp += rho_temp * (1 / pow(2,level)) * 0.5 * force_xyz(i,j,k,2,M_lim,K_lim,L_dims);
#endif

		// Assign velocity
		u(i,j,k,0,M_lim,K_lim,L_dims) = fux_temp / rho_temp;
		u(i,j,k,1,M_lim,K_lim,L_dims) = fuy_temp / rho_temp;
#if (L_dims == 3)
		u(i,j,k,2,M_lim,K_lim,L_dims) = fuz_temp / rho_temp;
#endif

	}

}


// ***************************************************************************************************
// Explosion operation
void GridObj::LBM_explode( int RegionNumber ) {

	// Declarations
	GridObj* fGrid = NULL;
	GridUtils::getGrid(MpiManager::Grids,level+1,RegionNumber,fGrid);
	int y_start, x_start, z_start;
	int M_fine = static_cast<int>(fGrid->YPos.size());
	int M_coarse = static_cast<int>(YPos.size());
	int K_coarse = static_cast<int>(ZPos.size());
	int K_fine = static_cast<int>(fGrid->ZPos.size());

	// Loop over coarse grid (just region of interest)
	for (int i = fGrid->CoarseLimsX[0]; i <= fGrid->CoarseLimsX[1]; i++) {
		for (int j = fGrid->CoarseLimsY[0]; j <= fGrid->CoarseLimsY[1]; j++) {
			for (int k = fGrid->CoarseLimsZ[0]; k <= fGrid->CoarseLimsZ[1]; k++) {

				// If TL to lower level and point belongs to region then explosion required
				if (LatTyp(i,j,k,M_coarse,K_coarse) == eTransitionToFiner) {

					// Lookup indices for lower level
					x_start = fGrid->CoarseLimsX[0];
					y_start = fGrid->CoarseLimsY[0];
					z_start = fGrid->CoarseLimsZ[0];

					// Find indices of fine site
					vector<int> idx_fine = GridUtils::getFineIndices(i, x_start, j, y_start, k, z_start);
					int fi = idx_fine[0];
					int fj = idx_fine[1];
					int fk = idx_fine[2];

					// Only pass information to lower levels if lower level is expecting it 
					// (i.e. not a boundary site)
					if (fGrid->LatTyp(fi,fj,fk,M_fine,K_fine) == eTransitionToCoarser) {

						// Update fine grid values according to Rohde et al.
						for (int v = 0; v < L_nVels; v++) {

							// Get coarse site value
							double coarse_f = f(i,j,k,v,M_coarse,K_coarse,L_nVels);

#if (L_dims == 3)
							// 3D Case -- cube of 8 cells

							// Copy coarse to fine
							fGrid->f(fi,	fj,		fk,		v,M_fine,K_fine,L_nVels)	= coarse_f;
							fGrid->f(fi+1,	fj,		fk,		v,M_fine,K_fine,L_nVels)	= coarse_f;
							fGrid->f(fi,	fj+1,	fk,		v,M_fine,K_fine,L_nVels)	= coarse_f;
							fGrid->f(fi+1,	fj+1,	fk,		v,M_fine,K_fine,L_nVels)	= coarse_f;
							fGrid->f(fi,	fj,		fk+1,	v,M_fine,K_fine,L_nVels)	= coarse_f;
							fGrid->f(fi+1,	fj,		fk+1,	v,M_fine,K_fine,L_nVels)	= coarse_f;
							fGrid->f(fi,	fj+1,	fk+1,	v,M_fine,K_fine,L_nVels)	= coarse_f;
							fGrid->f(fi+1,	fj+1,	fk+1,	v,M_fine,K_fine,L_nVels)	= coarse_f;

#else

							// 2D Case -- square of 4 cells

							// Copy coarse to fine
							fGrid->f(fi,	fj,		v,M_fine,L_nVels)		= coarse_f;
							fGrid->f(fi+1,	fj,		v,M_fine,L_nVels)		= coarse_f;
							fGrid->f(fi,	fj+1,	v,M_fine,L_nVels)		= coarse_f;
							fGrid->f(fi+1,	fj+1,	v,M_fine,L_nVels)		= coarse_f;

#endif
						}
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
	GridObj* fGrid;
	GridUtils::getGrid(MpiManager::Grids,level+1,RegionNumber,fGrid);
	int y_start, x_start, z_start;
	int M_fine = static_cast<int>(fGrid->YPos.size());
	int M_coarse = static_cast<int>(YPos.size());
	int K_coarse = static_cast<int>(ZPos.size());
#if (L_dims == 3)
	int K_fine = static_cast<int>(fGrid->ZPos.size());
#else
	int K_fine = 0;
#endif


	// Loop over coarse grid (only region of interest)
	for (int i = fGrid->CoarseLimsX[0]; i <= fGrid->CoarseLimsX[1]; i++) {
		for (int j = fGrid->CoarseLimsY[0]; j <= fGrid->CoarseLimsY[1]; j++) {
			for (int k = fGrid->CoarseLimsZ[0]; k <= fGrid->CoarseLimsZ[1]; k++) {

				// If TL to lower level or BC is on a refined level then fetch values from lower level
				if (LatTyp(i,j,k,M_coarse,K_coarse) == eTransitionToFiner || 
					LatTyp(i,j,k,M_coarse,K_coarse) == eRefinedSolid ||
					LatTyp(i,j,k,M_coarse,K_coarse) == eRefinedSymmetry ||
					LatTyp(i,j,k,M_coarse,K_coarse) == eRefinedInlet) {

					// Lookup indices for lower level
					x_start = fGrid->CoarseLimsX[0];
					y_start = fGrid->CoarseLimsY[0];
					z_start = fGrid->CoarseLimsZ[0];

					// Find indices of fine site
					vector<int> idx_fine = GridUtils::getFineIndices(i, x_start, j, y_start, k, z_start);
					int fi = idx_fine[0];
					int fj = idx_fine[1];
					int fk = idx_fine[2];

					// Loop over directions
					for (int v = 0; v < L_nVels; v++) {

						// Check to see if f value is missing on coarse level
						if (f(i,j,k,v,M_coarse,K_coarse,L_nVels) == 0) {

#if (L_dims == 3)
							// 3D Case -- cube of 8 cells

							// Average the values
							f(i,j,k,v,M_coarse,K_coarse,L_nVels) = (
								fGrid->f(fi,	fj,		fk,		v,M_fine,K_fine,L_nVels) +
								fGrid->f(fi+1,	fj,		fk,		v,M_fine,K_fine,L_nVels) +
								fGrid->f(fi,	fj+1,	fk,		v,M_fine,K_fine,L_nVels) +
								fGrid->f(fi+1,	fj+1,	fk,		v,M_fine,K_fine,L_nVels) +
								fGrid->f(fi,	fj,		fk+1,	v,M_fine,K_fine,L_nVels) +
								fGrid->f(fi+1,	fj,		fk+1,	v,M_fine,K_fine,L_nVels) +
								fGrid->f(fi,	fj+1,	fk+1,	v,M_fine,K_fine,L_nVels) +
								fGrid->f(fi+1,	fj+1,	fk+1,	v,M_fine,K_fine,L_nVels)
								) / pow(2, L_dims);

#else

							// 2D Case -- square of 4 cells

							// Average the values
							f(i,j,k,v,M_coarse,K_coarse,L_nVels) = (
								fGrid->f(fi,	fj,		v,M_fine,L_nVels) +
								fGrid->f(fi+1,	fj,		v,M_fine,L_nVels) +
								fGrid->f(fi,	fj+1,	v,M_fine,L_nVels) +
								fGrid->f(fi+1,	fj+1,	v,M_fine,L_nVels)
								) / pow(2, L_dims);

#endif
						}

					}

				}

			}
		}
	}

}

// ***************************************************************************************************
