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

#include "../inc/Matrix.h"

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
			subGrid[0]->LBM_multi_opt(i);
		}
	}

	// Start the clock to time this kernel
	clock_t secs, t_start = clock();

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

#ifdef L_LD_OUT
	// Reset object forces for momentum exchange force calculation
	objman->resetMomexBodyForces(this);
#endif

	// Loop over grid
	for (int i = 0; i < N_lim; ++i) {
		for (int j = 0; j < M_lim; ++j) {
			for (int k = 0; k < K_lim; ++k) {

				// Local index and type
				int id = k + j * K_lim + i * K_lim * M_lim;
				eType type_local = LatTyp[id];

				// MOMENTUM EXCHANGE //
#ifdef L_LD_OUT
				if (type_local == eSolid)
				{
					// Compute lift and drag contribution of this site
					objman->computeLiftDrag(i, j, k, this);
				}
#endif
				// IGNORE THESE SITES //
				if (type_local == eRefined || type_local == eSolid
#ifndef L_VELOCITY_REGULARISED
					|| type_local == eInlet
#endif
					) continue;

				// STREAM //
				_LBM_stream_opt(i, j, k, id, type_local, subcycle);

				// MACROSCOPIC //
				_LBM_macro_opt(i, j, k, id, type_local);

				// REGULARISED BCs //
#ifdef L_VELOCITY_REGULARISED
				if (type_local == eInlet)
					_LBM_regularised_velocity_opt(i, j, k, id);
#endif

				// FORCING //
#if (defined L_IBM_ON || defined L_GRAVITY_ON)
				// Do not force solid sites
				if (type_local != eSolid)
				{
					_LBM_forceGrid_opt(id);
				}
				L_DACTION_WRITE_OUT_FORCES
#endif
				// COLLIDE //
				if (type_local != eTransitionToCoarser) { // Do not collide on UpperTL

#ifdef L_USE_KBC_COLLISION
					_LBM_kbcCollide_opt(id);
#else
					_LBM_collide_opt(id);
#endif
				}

			}
		}
	}

	// Swap distributions
	f.swap(fNew);

#ifdef L_MOMEX_DEBUG
	if (level == objman->bbbOnGridLevel && region_number == objman->bbbOnGridReg)
	{
		// Close file for momentum exchange information (call before t increments)
		objman->toggleDebugStream(this);
	}
#endif

	// Increment internal loop counter
	++t;

	// Get time of loop
	secs = clock() - t_start;

	// Update average timestep time on this grid
	timeav_timestep *= (t - 1);
	timeav_timestep += ((double)secs) / CLOCKS_PER_SEC;
	timeav_timestep /= t;

	if (t % L_OUT_EVERY == 0) {
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
///	\param	type_local	type of current site.
///	\param	subcycle	number of sub-cycle being performed.
void GridObj::_LBM_stream_opt(int i, int j, int k, int id, eType type_local, int subcycle)
{

	// Local value to save multiple loads
	eType src_type_local;

	// Loop over velocities
	for (int v = 0; v < L_NUM_VELS; ++v)
	{

		// Forced equilibirum inlets and solid sites are simply reset to equilibrium
		if	(
			type_local == eSolid
#ifndef L_VELOCITY_REGULARISED
			|| type_local == eInlet
#endif
			)
		{
			fNew[v + id * L_NUM_VELS] = _LBM_equilibrium_opt(id, v);
			continue;
		}

		// Get indicies for source site (periodic by default)
		int src_x = (i - c_opt[v][0] + N_lim) % N_lim;
		int src_y = (j - c_opt[v][1] + M_lim) % M_lim;
		int src_z = (k - c_opt[v][2] + K_lim) % K_lim;

		// Source id and type
		int src_id = src_z + src_y * K_lim + src_x * K_lim * M_lim;
		src_type_local = LatTyp[src_id];

		// BFL BOUNCEBACK
		if (type_local == eBFL || src_type_local == eBFL)
		{
			// Try to apply BFL BC on streaming link
			if (_LBM_applyBFL_opt(id, src_id, v, i, j, k, src_x, src_y, src_z)) continue;
		}

		// BOUNCEBACK
		if (src_type_local == eSolid)
		{
			// F value is its opposite (HWBB)
			fNew[v + id * L_NUM_VELS] =
				f[GridUtils::getOpposite(v) + id * L_NUM_VELS];
		}

		// VELOCITY BC (forced equilbirium)
#ifndef L_VELOCITY_REGULARISED
		else if (src_type_local == eInlet)
		{
			// Set f to equilibrium (forced equilibrium BC)
			fNew[v + id * L_NUM_VELS] = _LBM_equilibrium_opt(src_id, v);
		}
#endif

#if (L_NUM_LEVELS > 0)	// Only need to check these options when using refinement
		// EXPLODE
		else if (src_type_local == eTransitionToCoarser && subcycle == 0)
		{
			// Pull value from parent TL site
			_LBM_explode_opt(id, v, src_x, src_y, src_z);
		}

		// COALESCE
		else if (src_type_local == eRefined && type_local == eTransitionToFiner)
		{
			// Pull average value from child TL cluster to get value leaving fine grid
			_LBM_coalesce_opt(i, j, k, id, v);
		}
#endif

		// REGULAR STREAM
		else
		{
			// Pull population from source site
			fNew[v + id * L_NUM_VELS] = f[v + src_id * L_NUM_VELS];
		}

	}

}

// *****************************************************************************
/// \brief	Optimised application of regularised velocity BC
///
///			This implementation assumes the velocity is normal to the wall as is
///			the case in the Latt et al. 2008 paper on which the conditions are
///			based. https://doi.org/10.1103/PhysRevE.77.056703
///
/// \param	i			x-index of current site.
/// \param	j			y-index of current site.
/// \param	k			z-index of current site.
///	\param	id			flattened ijk index.
void GridObj::_LBM_regularised_velocity_opt(int i, int j, int k, int id)
{
	// Declarations
	double normalVelocity;
	int normalVector[3] = { 0, 0, 0 };
	eCartesianDirection normalDirection;
	double f_plus = 0.0, f_zero = 0.0;
	double Sxx = 0, Syy = 0, Sxy = 0;	// 2D & 3D
	double fneq;
	double Szz = 0, Sxz = 0, Syz = 0;	// Just 3D
	unsigned int edgeCount = 0;

	// Find which wall we are on and proceed to assign key variables //

	// Left wall
	if (XPos[i] > 0.0 && XPos[i] < dh)
	{
		// Assign wall velocity
		normalVelocity = GridUnits::ud2ulbm(L_UX0, this);

		// Select the wall normal velocity direction
		normalDirection = eXDirection;
		normalVector[normalDirection] = 1;
		edgeCount++;
	}

	// Right wall
	if (XPos[i] < L_BX && XPos[i] > L_BX - dh)
	{
		normalVelocity = -GridUnits::ud2ulbm(L_UX0, this);
		normalDirection = eXDirection;
		normalVector[normalDirection] = -1;
		edgeCount++;
	}

	// Bottom wall
	else if (YPos[j] > 0.0 && YPos[j] < dh)
	{
		normalVelocity = GridUnits::ud2ulbm(L_UY0, this);
		normalDirection = eYDirection;
		normalVector[normalDirection] = 1;
		edgeCount++;
	}
	
	// Top wall
	else if (YPos[j] < L_BY && YPos[j] > L_BY - dh)
	{
		normalVelocity = -GridUnits::ud2ulbm(L_UY0, this);
		normalDirection = eYDirection;
		normalVector[normalDirection] = -1;
		edgeCount++;
	}

#if (L_DIMS == 3)

	// Front wall
	else if (ZPos[k] > 0.0 && ZPos[k] < dh)
	{
		normalVelocity = GridUnits::ud2ulbm(L_UZ0, this);
		normalDirection = eZDirection;
		normalVector[normalDirection] = 1;
		edgeCount++;
	}
	
	// Back wall
	else if (ZPos[k] < L_BZ && ZPos[k] > L_BZ - dh)
	{
		normalVelocity = -GridUnits::ud2ulbm(L_UZ0, this);
		normalDirection = eZDirection;
		normalVector[normalDirection] = -1;
		edgeCount++;
	}
#endif

	/* If it is a corner or an edge then get wall density by extrapolating in 
	 * the direction of the normal. Since this is non-local, when using MPI
	 * cannot be performed correctly if site on receiver layer */
	if (edgeCount > 1)
	{
		// Quadratic extrapolation (non-local)
		int i1 = i + normalVector[eXDirection];
		int i2 = i1 + normalVector[eXDirection];
		int j1 = j + normalVector[eYDirection];
		int j2 = j1 + normalVector[eYDirection];
		int k1 = k + normalVector[eZDirection];
		int k2 = k1 + normalVector[eZDirection];

#ifdef L_BUILD_FOR_MPI
		if (!GridUtils::isOnRecvLayer(XPos[i], YPos[j], ZPos[k]) &&
			(GridUtils::isOffGrid(i1, j1, k1, this) || GridUtils::isOffGrid(i2, j2, k2, this)))
		{
			std::string msg;
			msg += "Regularised BC edges/corners perform quadratic extrapolation.";
			msg += "The inlets on this block do not have enough points next to the inlet sites to perform this procedure. Exiting.";
			L_ERROR(msg, GridUtils::logfile);
		}
#endif
		rho[id] = 2.0 * rho(i1, j1, k1, M_lim, K_lim) - rho(i2, j2, k2, M_lim, K_lim);

	}

	// Only in one edge region so not a corner or an edge
	else
	{

		// Compute f_plus and f_0
		for (int v = 0; v < L_NUM_VELS; ++v)
		{
			// If has opposite normal direction component then part of f_plus
			if (c_opt[v][normalDirection] == -normalVector[normalDirection])
			{
				// Add to known momentum leaving the domain
				f_plus += fNew[v + id * L_NUM_VELS];

			}
			// If it is perpendicular to wall part of f_zero
			else if (c_opt[v][normalDirection] == 0)
			{
				f_zero += fNew[v + id * L_NUM_VELS];
			}

		}

		// Update density
		rho[id] = (1.0 / (1.0 - normalVelocity)) * (2.0 * f_plus + f_zero);
	}

	// Update velocity
	u[eXDirection + id * L_DIMS] = GridUnits::ud2ulbm(L_UX0, this);
	u[eYDirection + id * L_DIMS] = GridUnits::ud2ulbm(L_UY0, this);
#if (L_DIMS == 3)
	u[eZDirection + id * L_DIMS] = GridUnits::ud2ulbm(L_UZ0, this);
#endif

	// Loop over directions now macroscopic are up-to-date
	for (int v = 0; v < L_NUM_VELS; ++v)
	{

		// Apply off-equilibrium BB to unknown components //

		// Unknowns for a normal case share the normal vector components
		if (edgeCount == 1 && c_opt[v][normalDirection] == normalVector[normalDirection])
		{
			fNew[v + id * L_NUM_VELS] = _LBM_equilibrium_opt(id, v) +
				(fNew[GridUtils::getOpposite(v) + id * L_NUM_VELS] - _LBM_equilibrium_opt(id, GridUtils::getOpposite(v)));
		}

		// Unknown in edge cases are ones who share at least one of the normal components
		else if (edgeCount > 1 && (
			c_opt[v][eXDirection] == normalVector[eXDirection] ||
			c_opt[v][eYDirection] == normalVector[eYDirection]
#if (L_DIMS == 3)
			|| c_opt[v][eZDirection] == normalVector[eZDirection]
#endif
			)
			)
		{
			fNew[v + id * L_NUM_VELS] = _LBM_equilibrium_opt(id, v) +
				(fNew[GridUtils::getOpposite(v) + id * L_NUM_VELS] - _LBM_equilibrium_opt(id, GridUtils::getOpposite(v)));
		}

		// Store off-equilibrium and update stress components
		fneq = fNew[v + id * L_NUM_VELS] - _LBM_equilibrium_opt(id, v);

		// Compute off-equilibrium stress components
		Sxx += c_opt[v][eXDirection] * c_opt[v][eXDirection] * fneq;
		Syy += c_opt[v][eYDirection] * c_opt[v][eYDirection] * fneq;
		Sxy += c_opt[v][eXDirection] * c_opt[v][eYDirection] * fneq;
#if (L_DIMS == 3)
		Szz += c_opt[v][eZDirection] * c_opt[v][eZDirection] * fneq;
		Sxz += c_opt[v][eXDirection] * c_opt[v][eZDirection] * fneq;
		Syz += c_opt[v][eYDirection] * c_opt[v][eZDirection] * fneq;
#endif

	}

	// Compute regularised non-equilibrium components and add to feq to get new populations
	for (int v = 0; v < L_NUM_VELS; v++)
	{
		fNew[v + id * L_NUM_VELS] = _LBM_equilibrium_opt(id, v) +
			(w[v] / (2.0 * SQ(cs) * SQ(cs))) *
			(
			((c_opt[v][eXDirection] * c_opt[v][eXDirection] - SQ(cs)) * Sxx) +
			((c_opt[v][eYDirection] * c_opt[v][eYDirection] - SQ(cs)) * Syy) +
			((c_opt[v][eZDirection] * c_opt[v][eZDirection] - SQ(cs)) * Szz) +
			(2.0 * c_opt[v][eXDirection] * c_opt[v][eYDirection] * Sxy) +
			(2.0 * c_opt[v][eXDirection] * c_opt[v][eZDirection] * Sxz) +
			(2.0 * c_opt[v][eYDirection] * c_opt[v][eZDirection] * Syz)
			);
	}



}

// *****************************************************************************
/// \brief	Optimised coalesce operation.
///
/// \param	i	x-index of current site.
/// \param	j	y-index of current site.
/// \param	k	z-index of current site.
///	\param	id	flattened ijk index.
/// \param	v	lattice direction.
void GridObj::_LBM_coalesce_opt(int i, int j, int k, int id, int v) {

	// Get pointer to child grid
	GridObj *childGrid = subGrid[0];

	// Get sizes
	int cM_lim = childGrid->M_lim;
	int cK_lim = childGrid->K_lim;

	// Get indices of child site
	std::vector<int> cInd =
		GridUtils::getFineIndices(
		i, childGrid->CoarseLimsX[eMinimum],
		j, childGrid->CoarseLimsY[eMinimum],
		k, childGrid->CoarseLimsZ[eMinimum]);

	// Pull average value of f from child cluster
	double fNew_local = 0.0;
	for (int ii = 0; ii < 2; ++ii) {
		for (int jj = 0; jj < 2; ++jj) {
#if (L_DIMS == 3)
			for (int kk = 0; kk < 2; ++kk)
#else
			int kk = 0;
#endif
			{
				fNew_local +=
					childGrid->f[v +
					(cInd[2] + kk) * L_NUM_VELS +
					(cInd[1] + jj) * L_NUM_VELS * cK_lim +
					(cInd[0] + ii) * L_NUM_VELS * cK_lim * cM_lim];
			}
		}
	}

#if (L_DIMS == 3)
	fNew_local /= 8.0;
#else
	fNew_local /= 4.0;
#endif

	// Store back in memory
	fNew[v + id * L_NUM_VELS] = fNew_local;

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
		src_x, CoarseLimsX[eMinimum],
		src_y, CoarseLimsY[eMinimum],
		src_z, CoarseLimsZ[eMinimum]);

	// Pull value from parent
	fNew[v + id * L_NUM_VELS] =
		parentGrid->f[
			v +
				pInd[2] * L_NUM_VELS +
				pInd[1] * L_NUM_VELS * parentGrid->K_lim +
				pInd[0] * L_NUM_VELS * parentGrid->K_lim * parentGrid->M_lim
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

#if (L_DIMS == 3)
	A = (c_opt[v][0] * u[0 + id * L_DIMS]) +
		(c_opt[v][1] * u[1 + id * L_DIMS]) +
		(c_opt[v][2] * u[2 + id * L_DIMS]);

	B = (SQ(c_opt[v][0]) - SQ(cs)) * SQ(u[0 + id * L_DIMS]) +
		(SQ(c_opt[v][1]) - SQ(cs)) * SQ(u[1 + id * L_DIMS]) +
		(SQ(c_opt[v][2]) - SQ(cs)) * SQ(u[2 + id * L_DIMS]) +
		2 * c_opt[v][0] * c_opt[v][1] * u[0 + id * L_DIMS] * u[1 + id * L_DIMS] +
		2 * c_opt[v][0] * c_opt[v][2] * u[0 + id * L_DIMS] * u[2 + id * L_DIMS] +
		2 * c_opt[v][1] * c_opt[v][2] * u[1 + id * L_DIMS] * u[2 + id * L_DIMS];
#else
	A = (c_opt[v][0] * u[0 + id * L_DIMS]) +
		(c_opt[v][1] * u[1 + id * L_DIMS]);

	B = (SQ(c_opt[v][0]) - SQ(cs)) * SQ(u[0 + id * L_DIMS]) +
		(SQ(c_opt[v][1]) - SQ(cs)) * SQ(u[1 + id * L_DIMS]) +
		2 * c_opt[v][0] * c_opt[v][1] * u[0 + id * L_DIMS] * u[1 + id * L_DIMS];
#endif


	// Compute f^eq
	return rho[id] * w[v] * ( 1.0 + (A / SQ(cs)) + (B / (2.0 * SQ(cs)*SQ(cs)) ) );

}

// *****************************************************************************
/// \brief	Compute Smagorinksy-modified relaxation
/// 
///			Model taken from "DNS and LES of decaying isotropic turbulence with 
///			and without frame rotation using lattice Boltzmann method" by Yu, 
///			Huidan Girimaji, Sharath S. Luo, Li Shi  [2005]
///
///	\param	id 		flattened ijk index. 
/// \param 	omega 	Relaxation frequency. 
/// \return 		Smagorinsky-modified omega value
double GridObj::_LBM_smag(int id, double omega)
{
	// Calculate the non equilibrium stress tensor
	Matrix2D<double> nonEquiStress(3,3);
	double fneq[L_NUM_VELS];
 
	// Compute non-equilibrium values
	for (int v = 0; v < L_NUM_VELS; ++v)
		fneq[v] = fNew[v + id*L_NUM_VELS] - _LBM_equilibrium_opt(id, v);

	// Calculate diagonal and upper diagonal of the non equilibrium stress tensor
	for (int i = 0; i < L_DIMS; ++i)
	{
		for (int j = i; j < L_DIMS; ++j)
		{
			nonEquiStress[i][j] = 0.0;
			for (int v = 0; v < L_NUM_VELS; ++v)
			{
				nonEquiStress[i][j] += c_opt[v][i] * c_opt[v][j] * fneq[v];
			}
		}
	}

	// Fill in the lower diagonal of the non equilibrium stress tensor (since this tensor is symmetric)
	for (int i = 1; i < L_DIMS; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			nonEquiStress[i][j] = nonEquiStress[j][i];
		}
	}

	// Inner product
	double Q = sqrt(2.0 * (nonEquiStress % nonEquiStress));

	// Compute tau correction
	double tau = 1.0 / omega;
	double tau_t = 0.5 * (sqrt(SQ(tau) + 2.0 * L_SQRT2 * SQ(L_CSMAG) * L_RHOIN * SQ(cs) * SQ(cs) * Q ) - tau);  
	return ( 1.0 / (tau + tau_t) );
}

// *****************************************************************************
/// \brief	Optimised collision operation.
///
///			BGK collision operator. If Smagnorinksy turned on, will modify the 
///			value of omega locally.
///
/// \param	id	flattened ijk index.
void GridObj::_LBM_collide_opt(int id)
{

#ifdef L_USE_BGKSMAG
	// Compute Smagorinksy-modified relaxation
	double omega_s = _LBM_smag(id, omega);
#else
	double omega_s = omega;
#endif

	// Perform collision operation (using omega_s -- modified if using Smagorinksy)
	for (int v = 0; v < L_NUM_VELS; ++v)
	{
		fNew[v + id * L_NUM_VELS] +=
			omega_s *	(
			_LBM_equilibrium_opt(id, v) -
			fNew[v + id * L_NUM_VELS]
			) +
			force_i[v + id * L_NUM_VELS];
	}

}

// *****************************************************************************
/// \brief	Optimised macroscopic operation.
///
/// \param	i	x-index of current site.
/// \param	j	y-index of current site.
/// \param	k	z-index of current site.
/// \param	id	flattened ijk index.
///	\param	type_local	type of site under consideration
void GridObj::_LBM_macro_opt(int i, int j, int k, int id, eType type_local) {

	// Only update fluid sites (including BFL) or TL to finer
	if (type_local == eFluid || type_local == eBFL ||
		type_local == eTransitionToFiner) {

		// Reset
		double rho_temp = 0.0;
		double rhouX_temp = 0.0;
		double rhouY_temp = 0.0;
#if (L_DIMS == 3)
		double rhouZ_temp = 0.0;
#endif

		// Sum to find rho and momentum
		for (int v = 0; v < L_NUM_VELS; ++v) {
			rho_temp += fNew[v + id * L_NUM_VELS];
			rhouX_temp += c_opt[v][0] * fNew[v + id * L_NUM_VELS];
			rhouY_temp += c_opt[v][1] * fNew[v + id * L_NUM_VELS];
#if (L_DIMS == 3)
			rhouZ_temp += c_opt[v][2] * fNew[v + id * L_NUM_VELS];
#endif
		}

		// Add forces to momentum
#if (defined L_IBM_ON || defined L_GRAVITY_ON)
		rhouX_temp += 0.5 * force_xyz[0 + id * L_DIMS];
		rhouY_temp += 0.5 * force_xyz[1 + id * L_DIMS];
#if (L_DIMS == 3)
		rhouX_temp += 0.5 * force_xyz[2 + id * L_DIMS];
#endif
#endif

		// Divide by rho to get velocity
		u[0 + id * L_DIMS] = rhouX_temp / rho_temp;
		u[1 + id * L_DIMS] = rhouY_temp / rho_temp;
#if (L_DIMS == 3)
		u[2 + id * L_DIMS] = rhouZ_temp / rho_temp;
#endif

		// Assign density
		rho[id] = rho_temp;

	}

	// Update child TL sites for aethetic reasons only -- can be removed for performance
	if (type_local == eTransitionToFiner) {

		// Get child grid
		GridObj *childGrid = subGrid[0];

		// Get indices
		std::vector<int> cInd =
			GridUtils::getFineIndices(
			i, childGrid->CoarseLimsX[eMinimum],
			j, childGrid->CoarseLimsY[eMinimum],
			k, childGrid->CoarseLimsZ[eMinimum]);

		// Get sizes
		int cM_lim = childGrid->M_lim;
		int cK_lim = childGrid->K_lim;

		for (int ii = 0; ii < 2; ++ii) {
			for (int jj = 0; jj < 2; ++jj) {
#if (L_DIMS == 3)
				for (int kk = 0; kk < 2; ++kk)
#else
				int kk = 0;
#endif
				{
					for (int d = 0; d < L_DIMS; ++d) {
						childGrid->u[
							d +
								(cInd[2] + kk) * L_DIMS +
								(cInd[1] + jj) * L_DIMS * cK_lim +
								(cInd[0] + ii) * L_DIMS * cK_lim * cM_lim
						] = u[d + id * L_DIMS];
					}

					subGrid[0]->rho[
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
	for (int p = 0; p < L_DIMS; p++) {
		ta_temp = ui_timeav[p + id * L_DIMS] * (double)t;
		ta_temp += u[p + id * L_DIMS];
		ui_timeav[p + id * L_DIMS] = ta_temp / (double)(t + 1);
		// Do necessary products
		for (int q = p; q < L_DIMS; q++) {
			ta_temp = uiuj_timeav[pq_combo + id * (3 * L_DIMS - 3)] * (double)t;
			ta_temp += (u[p + id * L_DIMS] * u[q + id * L_DIMS]);
			uiuj_timeav[pq_combo + id * (3 * L_DIMS - 3)] = ta_temp / (double)(t + 1);
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
	// Add gravity and reset the lattice forces
	force_xyz[L_GRAVITY_DIRECTION + id * L_DIMS] =
		rho[id] * gravity * refinement_ratio;
	memset(&force_i[id * L_NUM_VELS], 0, sizeof(double) * L_NUM_VELS);
#endif

	// Now compute force_i components from Cartesian force vector
	for (size_t v = 0; v < L_NUM_VELS; v++) {

		// Reset beta_v
		beta_v = 0.0;

		// Compute the lattice forces based on Guo's forcing scheme
		lambda_v = (1 - 0.5 * omega) * (w[v] / (cs*cs));

		// Dot product (sum over d dimensions)
		for (int d = 0; d < L_DIMS; d++) {
			beta_v += (c_opt[v][d] * u[d + id * L_DIMS]);
		}
		beta_v = beta_v * (1 / (cs*cs));

		// Compute force using shorthand sum described above
		for (int d = 0; d < L_DIMS; d++) {
			force_i[v + id * L_NUM_VELS] += force_xyz[d + id * L_DIMS] *
				(c_opt[v][d] * (1 + beta_v) - u[d + id * L_DIMS]);
		}

		// Multiply by lambda_v
		force_i[v + id * L_NUM_VELS] *= lambda_v;

	}

}

// *****************************************************************************
/// \brief	Optimised BFL application.
///
///			Applies the BFL boundary condition at the present site. BFL site may
///			be either source or current site and so this method will decide
///			appropriate action based on each case as well as wall location.
///			Currently only first BFL body considered.
///
/// \param	id		flattened ijk index.
///	\param	src_id	flatted ijk index of the source site.
///	\param	v		lattice direction.
///	\param	i		x-index of current site.
///	\param	j		y-index of current site.
///	\param	k		z-index of current site.
///	\param	src_x	x-index of site where value is pulled from.
///	\param	src_y	y-index of site where value is pulled from.
///	\param	src_z	z-index of site where value is pulled from.
///	\return	boolen a indicator as to whether BFL was applied. If not, regular streaming is continued.
bool GridObj::_LBM_applyBFL_opt(int id, int src_id, int v, int i, int j, int k, int src_x, int src_y, int src_z)
{
	/* BFL can be applied easily if only one of the two sites are labelled as BFL.
	 * If both are labelled BFL, the algorithm here checks to see if there is any 
	 * intersecting wall assuming only one wall per voxel. If there are two 
	 * intersecting walls, then the BC favours the nearest. */

	// Initiate marker data store pointer on stack and retrieve Q
	MarkerData *m_data;
	double q_link = -1;		// Set to invalid value by default
	bool bCurrentSiteBflSite = true;
	int markerID;

	// Check whether current site is BFL site and get Q value
	if (LatTyp(i, j, k, M_lim, K_lim) == eBFL)
	{
		m_data = ObjectManager::getInstance()->pBody[0].getMarkerData(XPos[i], YPos[j], ZPos[k]);
		markerID = m_data->ID;
		delete m_data;
		q_link = ObjectManager::getInstance()->pBody[0].Q[GridUtils::getOpposite(v) + L_NUM_VELS * markerID];
		
	}

	/* If q value is valid then current site is a BFL site with link-intersecting 
	 * wall. If not, then we can check to see if the source site is a BFL site 
	 * and has a link-intersecting wall. */
	if (q_link == -1)
	{
		m_data = ObjectManager::getInstance()->pBody[0].getMarkerData(XPos[src_x], YPos[src_y], ZPos[src_z]);
		if (m_data->isValid())
		{
			bCurrentSiteBflSite = false;
			markerID = m_data->ID;
			q_link = ObjectManager::getInstance()->pBody[0].Q[v + L_NUM_VELS * markerID];
		}
		delete m_data;
	}
		
	/* BFL BC must only be applied if the pull link intersects the wall. Wall may
	 * or may not intersect the pull link hence this method must handle the case 
	 * where either the source site or the current site is BFL. If 
	 * q_link == -1, there is no wall intersection between the current
	 * site and the source site on that link. */

	// No intersection //
	if (q_link == -1)
	{
		// Regular stream as no intersection
		return false;
	}

	// Intersection and wall must be nearer current site than source site //
	else if (bCurrentSiteBflSite)
	{
		/* Here, the wall can be considered to be closer to BFL site but we need 
		 * to perform interpolation on pre-stream values pointing towards the
		 * wall from the BFL site and one site further away from the wall. This
		 * stencil site can be found from the pull direction unit vector. */
		int stencil_i = i + c_opt[v][0];
		int stencil_j = j + c_opt[v][1];
		int stencil_k = k + c_opt[v][2];
		int stencil_id = k + j * K_lim + i * K_lim * M_lim;

		/* Only apply if the stencil is valid on this rank. Otherwise, doesn't 
		 * matter as it is likely a halo site which will get overwritten anyway. */
		if (stencil_i >= 0 && stencil_i < N_lim &&
			stencil_j >= 0 && stencil_j < M_lim &&
			stencil_k >= 0 && stencil_k < K_lim)
		{
			// Interpolate pre-stream value then perform bounceback stream
			fNew[v + id * L_NUM_VELS] =
				(1 - 2 * q_link) *
				(f[GridUtils::getOpposite(v) + stencil_id * L_NUM_VELS] - f[GridUtils::getOpposite(v) + id * L_NUM_VELS])
				+ f[GridUtils::getOpposite(v) + id * L_NUM_VELS];

			// Momentum exchange -- don't include forces computed on halo sites to avoid duplicates
#ifdef L_LD_OUT
			if (!GridUtils::isOnRecvLayer(XPos[i], YPos[j], ZPos[k]))
				ObjectManager::getInstance()->computeLiftDrag(v, id, this, markerID);
#endif
		}
	}

	// Intersection and wall must be nearer source site than current site //
	else
	{
		/* Wall must be nearer the source site than the current site. We can 
		 * compute bounced value at current site from post-stream interpolated
		 * values pointing away from the wall. */
		fNew[v + id * L_NUM_VELS] =
			(1 - 2 * q_link) *
			((f[v + id * L_NUM_VELS] - f[GridUtils::getOpposite(v) + id * L_NUM_VELS]) / (2 - 2 * q_link))
			+ f[GridUtils::getOpposite(v) + id * L_NUM_VELS];

		// Momentum exchange -- don't include forces computed on halo sites to avoid duplicates
#ifdef L_LD_OUT
		if (!GridUtils::isOnRecvLayer(XPos[i], YPos[j], ZPos[k]))
			ObjectManager::getInstance()->computeLiftDrag(v, id, this, markerID);
#endif
	}

	return true;

}

// *****************************************************************************
/// \brief	Optimised KBC collision operator.
///
///			Applies KBC collision operator using the KBC-N4 and KBC-D models in 
///			3D and 2D, respectively.
///
/// \param id		flattened index of the lattice site.
void GridObj::_LBM_kbcCollide_opt(int id)
{

	// Declarations
	double ds[L_NUM_VELS];
	double dh[L_NUM_VELS];
	double fneq[L_NUM_VELS];
	double gamma;
	std::vector<double> Mneq;
	std::vector<int> C;

	// Compute required moments and equilibrium moments //
#if (L_DIMS == 3)
	int numMoments = 13;
#else
	int numMoments = 3;
#endif
	Mneq.resize(numMoments, 0.0);
	C.resize(numMoments * L_NUM_VELS, 1);

	for (int v = 0; v < L_NUM_VELS; v++)
	{

		// Update feq and store fneq
		feq[v + id * L_NUM_VELS] = _LBM_equilibrium_opt(id, v);
		fneq[v] = f[v + id * L_NUM_VELS] - feq[v + id * L_NUM_VELS];

		// 2-index and 3-index non-equilibrium moments
		int idx = 0;
		for (int sig = 0; sig < L_DIMS; ++sig)
		{
			for (int gam = sig; gam < L_DIMS; ++gam)
			{
				C[idx + v * numMoments] = c_opt[v][sig] * c_opt[v][gam];
				Mneq[idx] += fneq[v] * C[idx + v * numMoments];
				idx++;

#if (L_DIMS == 3)
				// Only need this inner loop for 3D moments
				for (int del = gam; del < L_DIMS; ++del)
				{
					// Don't include if all the same index
					if (sig != gam || gam != del || sig != del)
					{
						C[idx + v * numMoments] = c_opt[v][sig] * c_opt[v][gam] * c_opt[v][del];
						Mneq[idx] += fneq[v] * C[idx + v * numMoments];
						idx++;
					}
				}
#endif
			}
		}		
	}

	// Compute ds
	for (int v = 0; v < L_NUM_VELS; v++)
	{
		/* s part dictated by KBC model choice and directions.
		 * These are different in 2D and 3D so compile different options. */

#if (L_DIMS == 3)

		if (c_opt[v][0] == 0)
		{
			if (c_opt[v][1] == 0)
			{
				if (c_opt[v][2] == 0)
				{
					// First family
					ds[v] = (-(Mneq[0] + Mneq[8] + Mneq[12]));
				}
				else
				{
					// Fourth family
					ds[v] = ((-(Mneq[0] - Mneq[12]) - (Mneq[8] - Mneq[12])) / 6.0 + (Mneq[0] + Mneq[8] + Mneq[12]) / 6.0 - c_opt[v][2] * 0.5 * (Mneq[2] + Mneq[9]));
				}
			}
			else
			{
				if (c_opt[v][2] == 0)
				{
					// Third family
					ds[v] = ((-(Mneq[0] - Mneq[12]) + 2.0 * (Mneq[8] - Mneq[12])) / 6.0 + (Mneq[0] + Mneq[8] + Mneq[12]) / 6.0 - c_opt[v][1] * 0.5 * (Mneq[1] + Mneq[11]));
				}
				else
				{	// Seventh family
					ds[v] = (C[10 + v * numMoments] * 0.25 * Mneq[10] + (c_opt[v][2] * 0.25 * Mneq[9] + c_opt[v][1] * 0.25 * Mneq[11]));
				}
			}
		}
		else
		{
			if (c_opt[v][1] == 0)
			{
				if (c_opt[v][2] == 0)
				{
					// Second family
					ds[v] = ((2.0 * (Mneq[0] - Mneq[12]) - (Mneq[8] - Mneq[12])) / 6.0 + (Mneq[0] + Mneq[8] + Mneq[12]) / 6.0 - c_opt[v][0] * 0.5 * (Mneq[4] + Mneq[7]));
				}
				else
				{
					// Sixth family
					ds[v] = (C[6 + v * numMoments] * 0.25 * Mneq[6] + (c_opt[v][2] * 0.25 * Mneq[2] + c_opt[v][0] * 0.25 * Mneq[7]));
				}
			}
			else
			{
				if (c_opt[v][2] == 0)
				{
					// Fifth family
					ds[v] = (C[3 + v * numMoments] * 0.25 * Mneq[3] + (c_opt[v][1] * 0.25 * Mneq[1] + c_opt[v][0] * 0.25 * Mneq[4]));
				}
				else
				{
					// Eighth family
					ds[v] = (C[5 + v * numMoments] * Mneq[5] / 8.0);
				}
			}
		}
#else

		if (c_opt[v][0] == 0)
		{
			if (c_opt[v][1] == 0)
			{
				// First family
				ds[v] = 0.0;
			}
			else
			{
				// Third family
				ds[v] = -0.25 * (Mneq[0] - Mneq[2]);
			}
		}
		else
		{
			if (c_opt[v][1] == 0)
			{
				// Second family
				ds[v] = 0.25 * (Mneq[0] - Mneq[2]);
			}
			else
			{
				// Fourth family
				ds[v] = 0.25 * C[1 + v * numMoments] * Mneq[1];
			}
		}

#endif
		// Compute dh
		dh[v] = fneq[v] - ds[v];
	}

	// Once all dh and ds have been computed, compute the products
	double top_prod = 0.0, bot_prod = 0.0;
	for (int v = 0; v < L_NUM_VELS; v++)
	{
		// Compute scalar products
		top_prod += ds[v] * dh[v] / feq[v + id * L_NUM_VELS];
		bot_prod += dh[v] * dh[v] / feq[v + id * L_NUM_VELS];
	}

	// Compute 1/beta
	double beta_m1 = 2.0 / omega;

	// Compute gamma
	if (bot_prod == 0.0) gamma = 2.0;	// Regularised?
	else gamma = beta_m1 - (2.0 - beta_m1) * (top_prod / bot_prod);

	// Finally perform collision
	for (int v = 0; v < L_NUM_VELS; v++)
	{
		// Perform collision
		fNew[v + id * L_NUM_VELS] =
			f[v + id * L_NUM_VELS] -
			(1.0 / beta_m1) * (2.0 * ds[v] + gamma * dh[v]) +
			force_i[v + id * L_NUM_VELS];
	}

}
// *****************************************************************************