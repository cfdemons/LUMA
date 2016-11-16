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

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/ObjectManager.h"
#include "../inc/MpiManager.h"
#include "../inc/GridUtils.h"

// *****************************************************************************
/// \brief	Optimised LBM multi-grid kernel.
///
///			This kernel compresses the old kernel into a single loop in order to
///			make it more efficient. Capabilities are current limited with this 
///			kernel with incompatible options giving unpredictable results.
///			Use with caution.
///
///	\param	subcycle	sub-cycle to be performed if called from a subgrid.
void GridObj::LBM_multi_opt(int subcycle) {

	// Two iterations on sub-grid first
	if (subGrid.size()) {
		for (int i = 0; i < 2; ++i) {
			subGrid[0].LBM_multi_opt(i);
		}
	}

	// Start the clock to time this kernel
	clock_t secs, t_start = clock();

	// Loop over grid
	for (int i = 0; i < N_lim; ++i) {
		for (int j = 0; j < M_lim; ++j) {
			for (int k = 0; k < K_lim; ++k) {

				// Local index and type
				int id = k + j * K_lim + i * K_lim * M_lim;
				eType type_local = LatTyp[id];

				// Ignore Refined sites
				if (type_local == eRefined) continue;

				// STREAM //
				_LBM_stream_opt(i, j, k, id, subcycle);

				// MACROSCOPIC //
				_LBM_macro_opt(i, j, k, id, type_local);

				// FORCING //
#if (defined L_IBM_ON || defined L_GRAVITY_ON)
				if (type_local != eSolid) {	// Do not force solid sites
					_LBM_forceGrid_opt(id);
				}
				L_DACTION_WRITE_OUT_FORCES
#endif

					// COLLIDE //
					if (type_local != eTransitionToCoarser) { // Do not collide on UpperTL
						_LBM_collide_opt(id);
					}

			}
		}
	}

	// Swap distributions
	f.swap(fNew);

	// Increment internal loop counter
	++t;

	// Get time of loop
	secs = clock() - t_start;

	// Update average timestep time on this grid
	timeav_timestep *= (t - 1);
	timeav_timestep += ((double)secs) / CLOCKS_PER_SEC;
	timeav_timestep /= t;

	if (t % L_out_every == 0) {
		// Performance data to logfile
		*GridUtils::logfile << "Grid " << level << ": Time stepping taking an average of " << timeav_timestep * 1000 << "ms" << std::endl;
	}


	// MPI COMMUNICATION //
#ifdef L_BUILD_FOR_MPI

	// Launch communication on this grid by passing its level and region number
	MpiManager::getInstance()->mpi_communicate(level, region_number);

#endif

}

// *****************************************************************************
/// \brief	Optimised stream operation.
///
/// \param	i	x-index of current site.
/// \param	j	y-index of current site.
/// \param	k	z-index of current site.
///	\param	id	flattened ijk index.
///	\param	subcycle	number of sub-cycle being performed.
void GridObj::_LBM_stream_opt(int i, int j, int k, int id, int subcycle) {

	// Local value to save multiple loads
	eType src_type_local;

	for (int v = 0; v < L_nVels; ++v) {

		// Get indicies for source site (periodic by default)
		int src_x = (i - c_opt[v][0] + N_lim) % N_lim;
		int src_y = (j - c_opt[v][1] + M_lim) % M_lim;
		int src_z = (k - c_opt[v][2] + K_lim) % K_lim;

		// Source id and type
		int src_id = src_z + src_y * K_lim + src_x * K_lim * M_lim;
		src_type_local = LatTyp[src_id];

		// BOUNCEBACK
		if (src_type_local == eSolid) {
			// F value is its opposite (HWBB)
			fNew[v + id * L_nVels] =
				f[GridUtils::getOpposite(v) + id * L_nVels];
		}
		// VELOCITY BC
		else if (src_type_local == eInlet) {
			// Set f to equilibrium (forced equilibrium BC)
			fNew[v + id * L_nVels] = _LBM_equilibrium_opt(src_id, v);
		}
#if (L_NumLev > 0)
		// EXPLODE
		else if (src_type_local == eTransitionToCoarser &&
			subcycle == 0) {
			// Pull value from parent TL site
			_LBM_explode_opt(id, v, src_x, src_y, src_z);
		}
		// COALESCE
		else if (src_type_local == eRefined) {
			// Pull average value from child TL cluster to get value leaving fine grid
			_LBM_coalesce_opt(i, j, k, id, v);
		}
#endif
		// REGULAR STREAM
		else {
			// Pull population from source site
			fNew[v + id * L_nVels] = f[v + src_id * L_nVels];
		}
	}

}

// *****************************************************************************
/// \brief	Optimised coalesce operation.
///
/// \param	i	x-index of current site.
/// \param	j	y-index of current site.
/// \param	k	z-index of current site.
/// \param	v	lattice direction.
///	\param	region	region number from which values are to be coalesced.
void GridObj::_LBM_coalesce_opt(int i, int j, int k, int id, int v) {

	// Get pointer to child grid
	GridObj *childGrid = &subGrid[0];

	// Get sizes
	int cM_lim = childGrid->M_lim;
	int cK_lim = childGrid->K_lim;

	// Get indices of child site
	std::vector<int> cInd =
		GridUtils::getFineIndices(
		i, childGrid->CoarseLimsX[0],
		j, childGrid->CoarseLimsY[0],
		k, childGrid->CoarseLimsZ[0]);

	// Pull average value of f from child cluster
	double fNew_local = 0.0;
	for (int ii = 0; ii < 2; ++ii) {
		for (int jj = 0; jj < 2; ++jj) {
#if (L_dims == 3)
			for (int kk = 0; kk < 2; ++kk)
#else
			int kk = 0;
#endif
			{
				fNew_local +=
					childGrid->f[v +
					(cInd[2] + kk) * L_nVels +
					(cInd[1] + jj) * L_nVels * cK_lim +
					(cInd[0] + ii) * L_nVels * cK_lim * cM_lim];
			}
		}
	}

#if (L_dims == 3)
	fNew_local /= 8.0;
#else
	fNew_local /= 4.0;
#endif

	// Store back in memory
	fNew[v + id * L_nVels] = fNew_local;
}

// *****************************************************************************
/// \brief	Optimised explode operation.
///
/// \param	id	flattened ijk index.
///	\param	v	lattice direction.
///	\param	src_x	x-index of site where value is pulled from.
///	\param	src_y	y-index of site where value is pulled from.
///	\param	src_z	z-index of site where value is pulled from.
void GridObj::_LBM_explode_opt(int id, int v, int src_x, int src_y, int src_z) {

	// Get parent indices
	std::vector<int> pInd =
		GridUtils::getCoarseIndices(
		src_x, CoarseLimsX[0],
		src_y, CoarseLimsY[0],
		src_z, CoarseLimsZ[0]);

	// Pull value from parent
	fNew[v + id * L_nVels] =
		parentGrid->f[
			v +
				pInd[2] * L_nVels +
				pInd[1] * L_nVels * parentGrid->K_lim +
				pInd[0] * L_nVels * parentGrid->K_lim * parentGrid->M_lim
		];
}

// *****************************************************************************
/// \brief	Optimised equilibrium calculation.
///
///			Computes the equilibrium distribution in direction supplied at the 
///			given lattice site and returns the value.
///
/// \param id	flattened ijk index.
/// \param v	lattice direction.
/// \return		equilibrium function.
double GridObj::_LBM_equilibrium_opt(int id, int v) {

	// Declare intermediate values A and B
	double A, B;

	// Compute the parts of the expansion for feq

#if (L_dims == 3)
	A = (c_opt[v][0] * u[0 + id * L_dims]) +
		(c_opt[v][1] * u[1 + id * L_dims]) +
		(c_opt[v][2] * u[2 + id * L_dims]);

	B = ((c_opt[v][0] * c_opt[v][0]) - (cs*cs)) * (u[0 + id * L_dims] * u[0 + id * L_dims]) +
		((c_opt[v][1] * c_opt[v][1]) - (cs*cs)) * (u[1 + id * L_dims] * u[1 + id * L_dims]) +
		((c_opt[v][2] * c_opt[v][2]) - (cs*cs)) * (u[2 + id * L_dims] * u[2 + id * L_dims]) +
		2 * c_opt[v][0] * c_opt[v][1] * u[0 + id * L_dims] * u[1 + id * L_dims] +
		2 * c_opt[v][0] * c_opt[v][2] * u[0 + id * L_dims] * u[2 + id * L_dims] +
		2 * c_opt[v][1] * c_opt[v][2] * u[1 + id * L_dims] * u[2 + id * L_dims];
#else
	A = (c_opt[v][0] * u[0 + id * L_dims]) +
		(c_opt[v][1] * u[1 + id * L_dims]);

	B = ((c_opt[v][0] * c_opt[v][0]) - (cs*cs)) * (u[0 + id * L_dims] * u[0 + id * L_dims]) +
		((c_opt[v][1] * c_opt[v][1]) - (cs*cs)) * (u[1 + id * L_dims] * u[1 + id * L_dims]) +
		2 * c_opt[v][0] * c_opt[v][1] * u[0 + id * L_dims] * u[1 + id * L_dims];
#endif


	// Compute f^eq
	return rho[id] * w[v] * (1 + (A / (cs*cs)) + (B / (2 * (cs*cs*cs*cs))));

}

// *****************************************************************************
/// \brief	Optimised collision operation.
///
///			Uses either BGK or KBC depending on definitions. KBC implementation
///			is currently not optimised.
///
/// \param	id	flattened ijk index.
void GridObj::_LBM_collide_opt(int id) {

#ifdef L_USE_KBC_COLLISION
	LBM_kbcCollide();

#else
	// Perform collision operation
	for (int v = 0; v < L_nVels; ++v) {
		fNew[v + id * L_nVels] +=
			omega *	(
			_LBM_equilibrium_opt(id, v) -
			fNew[v + id * L_nVels]
			) +
			force_i[v + id * L_nVels];
	}
#endif

}

// *****************************************************************************
/// \brief	Optimised macroscopic operation.
///
/// \param	id	flattened ijk index.
void GridObj::_LBM_macro_opt(int i, int j, int k, int id, eType type_local) {

	// Only update fluid sites or TL to finer
	if (type_local == eFluid ||
		type_local == eTransitionToFiner) {

		// Reset
		double rho_temp = 0.0;
		double rhouX_temp = 0.0;
		double rhouY_temp = 0.0;
#if (L_dims == 3)
		double rhouZ_temp = 0.0;
#endif

		// Sum to find rho and momentum
		for (int v = 0; v < L_nVels; ++v) {
			rho_temp += fNew[v + id * L_nVels];
			rhouX_temp += c_opt[v][0] * fNew[v + id * L_nVels];
			rhouY_temp += c_opt[v][1] * fNew[v + id * L_nVels];
#if (L_dims == 3)
			rhouZ_temp += c_opt[v][2] * fNew[v + id * L_nVels];
#endif
		}

		// Add forces to momentum
#if (defined L_IBM_ON || defined L_GRAVITY_ON)
		rhouX_temp += 0.5 * force_xyz[0 + id * L_dims];
		rhouY_temp += 0.5 * force_xyz[1 + id * L_dims];
#if (L_dims == 3)
		rhouX_temp += 0.5 * force_xyz[2 + id * L_dims];
#endif
#endif

		// Divide by rho to get velocity
		u[0 + id * L_dims] = rhouX_temp / rho_temp;
		u[1 + id * L_dims] = rhouY_temp / rho_temp;
#if (L_dims == 3)
		u[2 + id * L_dims] = rhouZ_temp / rho_temp;
#endif

		// Assign density
		rho[id] = rho_temp;

	}

	// Update child TL sites for aethetic reasons only -- can be removed for performance
	if (type_local == eTransitionToFiner) {

		// Get child grid
		GridObj *childGrid = &subGrid[0];

		// Get indices
		std::vector<int> cInd =
			GridUtils::getFineIndices(
			i, childGrid->CoarseLimsX[0],
			j, childGrid->CoarseLimsY[0],
			k, childGrid->CoarseLimsZ[0]);

		// Get sizes
		int cM_lim = childGrid->M_lim;
		int cK_lim = childGrid->K_lim;

		for (int ii = 0; ii < 2; ++ii) {
			for (int jj = 0; jj < 2; ++jj) {
#if (L_dims == 3)
				for (int kk = 0; kk < 2; ++kk)
#else
				int kk = 0;
#endif
				{
					for (int d = 0; d < L_dims; ++d) {
						childGrid->u[
							d +
								(cInd[2] + kk) * L_dims +
								(cInd[1] + jj) * L_dims * cK_lim +
								(cInd[0] + ii) * L_dims * cK_lim * cM_lim
						] = u[d + id * L_dims];
					}

					subGrid[0].rho[
						(cInd[2] + kk) +
							(cInd[1] + jj) * cK_lim +
							(cInd[0] + ii) * cK_lim * cM_lim
					] = rho[id];
				}
			}
		}
	}

	// TIME-AVERAGED QUANTITIES //

#ifdef L_COMPUTE_TIME_AVERAGED_QUANTITIES
	// Multiply current value by completed time steps to get sum
	double ta_temp = rho_timeav[id] * (double)t;
	// Add new value
	ta_temp += rho[id];
	// Divide by completed time steps + 1 to get new average
	rho_timeav[id] = ta_temp / (double)(t + 1);

	// Repeat for other quantities
	int pq_combo = 0;
	for (int p = 0; p < L_dims; p++) {
		ta_temp = ui_timeav[p + id * L_dims] * (double)t;
		ta_temp += u[p + id * L_dims];
		ui_timeav[p + id * L_dims] = ta_temp / (double)(t + 1);
		// Do necessary products
		for (int q = p; q < L_dims; q++) {
			ta_temp = uiuj_timeav[pq_combo + id * (3 * L_dims - 3)] * (double)t;
			ta_temp += (u[p + id * L_dims] * u[q + id * L_dims]);
			uiuj_timeav[pq_combo + id * (3 * L_dims - 3)] = ta_temp / (double)(t + 1);
			pq_combo++;
		}
	}
#endif

}

// *****************************************************************************
/// \brief	Optimised body force calculator.
///
///			Takes Cartesian force vector and populates forces for each lattice 
///			direction. If reset_flag is true, then resets the force vectors to zero.
///
///	\param	id	flattened ijk index.
void GridObj::_LBM_forceGrid_opt(int id) {

	/* This routine computes the forces applied along each direction on the lattice
	from Guo's 2002 scheme. The basic LBM must be modified in two ways: 1) the forces
	are added to the populations produced by the collision in the collision
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

	F_i = lambda_i * sum( F_d * (c_d_i (1 +beta_i) - u_d )

	*/

	// Declarations
	double lambda_v, beta_v;

#ifdef L_GRAVITY_ON
	// Add gravity
	force_xyz[L_grav_direction + id * L_dims] =
		rho[id] * L_grav_force * refinement_ratio;
#endif

	// Now compute force_i components from Cartesian force vector
	for (size_t v = 0; v < L_nVels; v++) {

		// Reset beta_v
		beta_v = 0.0;

		// Compute the lattice forces based on Guo's forcing scheme
		lambda_v = (1 - 0.5 * omega) * (w[v] / (cs*cs));

		// Dot product (sum over d dimensions)
		for (int d = 0; d < L_dims; d++) {
			beta_v += (c_opt[v][d] * u[d + id * L_dims]);
		}
		beta_v = beta_v * (1 / (cs*cs));

		// Compute force using shorthand sum described above
		for (int d = 0; d < L_dims; d++) {
			force_i[v + id * L_nVels] += force_xyz[d + id * L_dims] *
				(c_opt[v][d] * (1 + beta_v) - u[d + id * L_dims]);
		}

		// Multiply by lambda_v
		force_i[v + id * L_nVels] *= lambda_v;

	}

}
// *****************************************************************************