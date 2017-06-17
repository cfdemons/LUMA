/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) The University of Manchester 2017
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * further distribution commericially or otherwise without written consent.
 *
 */

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/ObjectManager.h"



// *****************************************************************************
/// Perform IBM procedure.
void ObjectManager::ibm_apply2(int level) {

//	// Loop through all bodies and find support
//	for (int ib = 0; ib < iBody.size(); ib++) {
//
//		// If a flexible body recalculate support and delta values
//		if (iBody[ib].isFlexible) {
//
//			// Recompute support for new marker positions
//			for (int m = 0; m < static_cast<int>(iBody[ib].markers.size()); m++) {
//
//				// Erase support vectors
//				iBody[ib].markers[m].supp_i.clear();
//				iBody[ib].markers[m].supp_j.clear();
//				iBody[ib].markers[m].supp_k.clear();
//				iBody[ib].markers[m].supp_x.clear();
//				iBody[ib].markers[m].supp_y.clear();
//				iBody[ib].markers[m].supp_z.clear();
//				iBody[ib].markers[m].deltaval.clear();
//				iBody[ib].markers[m].support_rank.clear();
//
//				// Get the marker position
//				double x = iBody[ib].markers[m].position[eXDirection];
//				double y = iBody[ib].markers[m].position[eYDirection];
//				double z = iBody[ib].markers[m].position[eZDirection];
//
//				// Get ijk of enclosing voxel
//				std::vector<int> ijk;
//				GridUtils::getEnclosingVoxel(x, y, z, iBody[ib]._Owner, &ijk);
//				iBody[ib].markers[m].supp_i.push_back(ijk[eXDirection]);
//				iBody[ib].markers[m].supp_j.push_back(ijk[eYDirection]);
//				iBody[ib].markers[m].supp_k.push_back(ijk[eZDirection]);
//
//				// Get the x-position of the support marker
//				iBody[ib].markers[m].supp_x.push_back(iBody[ib]._Owner->XPos[ijk[eXDirection]]);
//				iBody[ib].markers[m].supp_y.push_back(iBody[ib]._Owner->YPos[ijk[eYDirection]]);
//				iBody[ib].markers[m].supp_z.push_back(iBody[ib]._Owner->ZPos[ijk[eZDirection]]);
//
//				// Set rank of first support marker and marker ID in body
//				iBody[ib].markers[m].support_rank.push_back(GridUtils::safeGetRank());
//
//				// Recompute support
//				ibm_findSupport(ib, m);
//			}
//		}
//	}
//
//	// Calculate epsilon
//	if (hasMovingBodies)
//		ibm_findEpsilon();

	// Interpolate the velocity onto the markers
	ibm_interpolate(level);

	// Compute force
	ibm_computeForce(level);


#ifdef L_BUILD_FOR_MPI
	MpiManager *mpim = MpiManager::getInstance();
	MPI_Barrier(mpim->world_comm);
#endif
	exit(0);
}


// *****************************************************************************
/// Perform IBM procedure.
void ObjectManager::ibm_apply() {

//#ifdef L_BUILD_FOR_MPI
//	// Get MPI Manager Instance
//	MpiManager *mpim = MpiManager::getInstance();
//#endif
//	int rank = GridUtils::safeGetRank();
//
//	// Loop over array of IB_bodies and perform IB operations
//	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {
//
//#ifdef L_IBM_DEBUG
//		// DEBUG -- write out support coordinates
//		std::ofstream suppout;
//		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
//			suppout.open(GridUtils::path_str + "/Supp_" + std::to_string(iBody[ib].id) + "_" + std::to_string(m) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
//			suppout << "\nNEW TIME STEP" << std::endl;
//			suppout << "x\ty\tz" << std::endl;
//			for (size_t i = 0; i < iBody[ib].markers[m].supp_i.size(); i++) {
//				if (L_DIMS == 3) {
//					suppout << iBody[ib]._Owner->XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << iBody[ib]._Owner->YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << iBody[ib]._Owner->ZPos[iBody[ib].markers[m].supp_k[i]] << std::endl;
//				} else {
//					suppout << iBody[ib]._Owner->XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << iBody[ib]._Owner->YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << 0.0 << std::endl;
//				}
//			}
//			suppout.close();
//		}
//#endif
//
//#ifdef L_IBM_DEBUG
//		// DEBUG -- write out epsilon values
//		std::ofstream epout;
//		epout.open(GridUtils::path_str + "/Epsilon_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
//		epout << "\nNEW TIME STEP" << std::endl;
//		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
//			epout << iBody[ib].markers[m].epsilon << std::endl;
//		}
//		epout.close();
//#endif
//
//			// Interpolate velocity
//			ibm_interpol(ib);
//
//#ifdef L_IBM_DEBUG
//		// DEBUG -- write out interpolate velocity values
//		std::ofstream predout;
//		predout.open(GridUtils::path_str + "/interpVel_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
//		predout << "\nNEW TIME STEP" << std::endl;
//		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
//			predout << iBody[ib].markers[m].interpVel[0] << "\t" << iBody[ib].markers[m].interpVel[1] << std::endl;
//		}
//		predout.close();
//#endif
//
//			// Compute restorative force
//			ibm_computeForce(ib);
//
//#ifdef L_IBM_DEBUG
//		// DEBUG -- write out Lagrange force values
//		std::ofstream forceout;
//		forceout.open(GridUtils::path_str + "/force_xyz_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
//		forceout << "\nNEW TIME STEP" << std::endl;
//		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
//			forceout << iBody[ib].markers[m].force_xyz[0] << "\t" << iBody[ib].markers[m].force_xyz[1] << std::endl;
//		}
//		forceout.close();
//#endif
//
//			// Spread force back to lattice (Cartesian vector)
//			ibm_spread(ib);
//
//		}
}

// *****************************************************************************
/// \brief	Moves iBodies after applying IBM.
///
///			Wrapper for relocating markers of an iBody be calling appropriate
///			positional update routine.
///
void ObjectManager::ibm_moveBodies() {

//	// Loop over bodies launching positional update if movable to compute new locations of markers
//	*GridUtils::logfile << "Relocating markers as required..." << std::endl;
//	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {
//
//		// If body is movable it needs a positional update
//		if (iBody[ib].isMovable) {
//
//			// Call structural or forced positional update and recompute support
//			ibm_positionUpdate(ib);
//
//#ifndef L_STOP_EPSILON_RECOMPUTE
//			// Recompute epsilon
//			ibm_findEpsilon();
//#endif
//
//		}
//	}
//
//#if defined L_INSERT_FILARRAY
//	// Special bit for filament-based plates where flexible centreline is used to update position of others in group
//	*GridUtils::logfile << "Filament-based plate positional update..." << std::endl;
//	ibm_positionUpdateGroup(999);
//#endif

}


// *****************************************************************************
/// \brief	Initialise the array of iBodies.
///
///			Computes support and epsilon values.
///
void ObjectManager::ibm_initialise() {

	// Get rank safely
	int rank = GridUtils::safeGetRank();

	// Loop over the number of bodies in the iBody array
	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {

		// Write out marker positions for debugging
#ifdef L_IBM_DEBUG
		ibm_debug_markerPosition(ib);
#endif

		// Compute extended support for each marker in the body (or body portion)
		for (int m = 0; m < static_cast<int>(iBody[ib].markers.size()); m++)
			ibm_findSupport(ib, m);
	}

	// Get MPI Manager Instance
#ifdef L_BUILD_FOR_MPI
	ibm_buildMPIComms();
#endif

	// Find epsilon for the body
	ibm_findEpsilon();

	// Write out epsilon
#ifdef L_IBM_DEBUG
	for (int ib = 0; ib < iBody.size(); ib++)
		ibm_debug_epsilon(ib);
#endif
}

// *****************************************************************************
/// \brief	Method to evaluate delta kernel at supplied location.
///
///			Radius and dilation must be in the same units.
///
/// \param	radius		location at which kernel should be evaluated.
/// \param	dilation	width of kernel function.
/// \return	value of kernel function.
double ObjectManager::ibm_deltaKernel(double radius, double dilation) {

	double mag_r, value;

	// Absolute value of radius
	mag_r = fabs(radius) / dilation;

	// Piecemeal function evaluation
	if (mag_r > 1.5) {
		value = 0.0;
	} else if (mag_r > 0.5) {
		value = (5.0 - (3.0 * mag_r) - sqrt(-3.0 * pow( 1.0 - mag_r, 2) + 1.0)) / 6.0;
	} else {
		value = (1.0 + sqrt(1.0 - 3.0 * pow(mag_r, 2) ) ) / 3.0;
	}

	return value;
}

// *****************************************************************************
/// \brief	Finds support points for iBody.
///
///			Support for given marker in given body is sought on the owning grid.
///
/// \param	ib	body under consideration.
/// \param	m	marker whose support is to be found.
void ObjectManager::ibm_findSupport(int ib, int m) {

	// Get the rank
	int rank = GridUtils::safeGetRank();

	// Declarations
	int inear, jnear;							// Nearest node indices (global indices)
	std::vector<double> nearpos;				// Position of the nearest marker
	int knear;
	int s = 0;										// The sth support marker which has been found

#ifdef L_BUILD_FOR_MPI
	MpiManager *mpim = MpiManager::getInstance();
	int estimated_rank_offset[3] = { 0, 0, 0 };
#endif

	// Find closest support node
	std::vector<int> nearijk;
	GridUtils::getEnclosingVoxel(
		iBody[ib].markers[m].position[eXDirection],
		iBody[ib].markers[m].position[eYDirection],
		iBody[ib].markers[m].position[eZDirection],
		iBody[ib]._Owner,
		&nearijk
		);

	// Indices of closest support
	inear = nearijk[eXDirection];
	jnear = nearijk[eYDirection];
	knear = nearijk[eZDirection];

	// Get position of nearest node
	nearpos.push_back(iBody[ib]._Owner->XPos[inear]);
	nearpos.push_back(iBody[ib]._Owner->YPos[jnear]);
	nearpos.push_back(iBody[ib]._Owner->ZPos[knear]);

	// Side length of support region defined as 3 x dilation parameter
	iBody[ib].markers[m].dilation = 1.0;

	// Store normalised area of support region
	iBody[ib].markers[m].local_area = 1.0;

	// Vector to store estimated positions of the possible support points
	double estimated_position[3] = { 0, 0, 0 };

	// Get the position for the first support point (which already exists from the addMarker method)
	estimated_position[eXDirection] = nearpos[eXDirection];
	estimated_position[eYDirection] = nearpos[eYDirection];
#if (L_DIMS == 3)
	estimated_position[eZDirection] = nearpos[eZDirection];
#endif

	// Get the deltaval for the first support point
	ibm_initialiseSupport(ib, m, s, estimated_position);

	// Write out info for support site
#ifdef L_IBM_DEBUG
	ibm_debug_supportInfo(ib, m, s);
#endif

	// Loop over surrounding 5 lattice sites and check if within support region
	for (int i = inear - 5; i <= inear + 5; i++) {
		for (int j = jnear - 5; j <= jnear + 5; j++) {
#if (L_DIMS == 3)
			for (int k = knear - 5; k <= knear + 5; k++)
#else
			int k = 0;
#endif
			{
				/* Estimate position of support point rather than read from the 
				 * grid in case the point is outside the rank. 
				 * Estimate only works since LBM lattice uniformly spaced. */				
				estimated_position[eXDirection] = nearpos[eXDirection] + (i - inear) * iBody[ib]._Owner->dh;
				estimated_position[eYDirection] = nearpos[eYDirection] + (j - jnear) * iBody[ib]._Owner->dh;
				estimated_position[eZDirection] = nearpos[eZDirection] + (k - knear) * iBody[ib]._Owner->dh;

				/* Find distance between Lagrange marker and proposed support point and
				 * Check if inside the cage (convert to lattice units) */
				if	(
					(fabs(iBody[ib].markers[m].position[eXDirection] - estimated_position[eXDirection]) / iBody[ib]._Owner->dh
					< 1.5 * iBody[ib].markers[m].dilation)
					&&
					(fabs(iBody[ib].markers[m].position[eYDirection] - estimated_position[eYDirection]) / iBody[ib]._Owner->dh
					< 1.5 * iBody[ib].markers[m].dilation)
#if (L_DIMS == 3)
					&&
					(fabs(iBody[ib].markers[m].position[eZDirection] - estimated_position[eZDirection]) / iBody[ib]._Owner->dh 
					< 1.5 * iBody[ib].markers[m].dilation)
#endif
					)
				{

					// Skip the nearest as already added when marker constructed
					if (i == inear && j == jnear
#if (L_DIMS == 3)
						&& k == knear
#endif			
						)
					{
						continue;
					}
					else
					{

						// Increment the support point counter
						s++;

						// Lies within support region so add support point data
						iBody[ib].markers[m].supp_i.push_back(i);
						iBody[ib].markers[m].supp_j.push_back(j);
						iBody[ib].markers[m].supp_k.push_back(k);

						iBody[ib].markers[m].supp_x.push_back(estimated_position[eXDirection]);
						iBody[ib].markers[m].supp_y.push_back(estimated_position[eYDirection]);
						iBody[ib].markers[m].supp_z.push_back(estimated_position[eZDirection]);

						// Add owning rank as this one for now
						iBody[ib].markers[m].support_rank.push_back(rank);

#ifdef L_BUILD_FOR_MPI
						/* Estimate which rank this point belongs to by seeing which
						 * edge of the grid it is off. Use estimated rather than
						 * actual positions. Compare to sender layer edges as if
						 * it is on the recv layer it is belongs to the neighbour. */
						if (estimated_position[eXDirection] < mpim->sender_layer_pos.X[eLeftMin])
							estimated_rank_offset[eXDirection] = -1;
						if (estimated_position[eXDirection] > mpim->sender_layer_pos.X[eRightMax])
							estimated_rank_offset[eXDirection] = 1;
						if (estimated_position[eYDirection] < mpim->sender_layer_pos.Y[eLeftMin])
							estimated_rank_offset[eYDirection] = -1;
						if (estimated_position[eYDirection] > mpim->sender_layer_pos.Y[eRightMax])
							estimated_rank_offset[eYDirection] = 1;
#if (L_DIMS == 3)
						if (estimated_position[eZDirection] < mpim->sender_layer_pos.Z[eLeftMin])
							estimated_rank_offset[eZDirection] = -1;
						if (estimated_position[eZDirection] > mpim->sender_layer_pos.Z[eRightMax])
							estimated_rank_offset[eZDirection] = 1;
#endif

						// Get MPI direction of the neighbour that owns this point
						int owner_direction = GridUtils::getMpiDirection(estimated_rank_offset);
						if (owner_direction != -1) {

							// Owned by a neighbour so correct the support rank
							iBody[ib].markers[m].support_rank.back() = mpim->neighbour_rank[owner_direction];
						}

						// Reset estimated rank offset
						estimated_rank_offset[eXDirection] = 0;
						estimated_rank_offset[eYDirection] = 0;
#if (L_DIMS == 3)
						estimated_rank_offset[eZDirection] = 0;
#endif
#endif
						/* Initialise delta information for the set of support points
							 * including those not on this rank using estimated positions. */
						ibm_initialiseSupport(ib, m, s, estimated_position);

						// Write out info for support site
#ifdef L_IBM_DEBUG
						ibm_debug_supportInfo(ib, m, s);
#endif
					}
				}
			}
		}
	}
}

// *****************************************************************************
/// \brief	Initialise data associated with support points found.
///
///			Finds and stores the delta values of the support points.
///
/// \param	ib	iBody being operated on.
/// \param	m	marker of interest.
/// \param	s	support point of interest.
///	\param	estimated_position	vector containing the estimated position of the support point.
void ObjectManager::ibm_initialiseSupport(int ib, int m, int s, double estimated_position[])
{

	// Declarations
	double dist_x, dist_y, delta_x, delta_y;	// Distances and deltas
#if (L_DIMS == 3)
	double dist_z, delta_z;
#endif

	// Distance between Lagrange marker and support node in lattice units
	dist_x = (estimated_position[eXDirection] - iBody[ib].markers[m].position[eXDirection]) / iBody[ib]._Owner->dh;
	dist_y = (estimated_position[eYDirection] - iBody[ib].markers[m].position[eYDirection]) / iBody[ib]._Owner->dh;
#if (L_DIMS == 3)
	dist_z = (estimated_position[eZDirection] - iBody[ib].markers[m].position[eZDirection]) / iBody[ib]._Owner->dh;
#endif

	// Store delta function value
	delta_x = ibm_deltaKernel(dist_x, iBody[ib].markers[m].dilation);
	delta_y = ibm_deltaKernel(dist_y, iBody[ib].markers[m].dilation);
#if (L_DIMS == 3)
	delta_z = ibm_deltaKernel(dist_z, iBody[ib].markers[m].dilation);
#endif

	iBody[ib].markers[m].deltaval.push_back(delta_x * delta_y
#if (L_DIMS == 3)
		* delta_z
#endif	
		);
}

// *****************************************************************************
/// \brief	Interpolate velocity field onto markers
void ObjectManager::ibm_interpolate(int level) {

	// Get rank
	int rank = GridUtils::safeGetRank();

	// Loop through all bodies
	for (int ib = 0; ib < iBody.size(); ib++) {

		// Only interpolate the bodies that exist on this grid level
		if (iBody[ib]._Owner->level == level) {

			// Get grid sizes
			size_t M_lim = iBody[ib]._Owner->M_lim;
#if (L_DIMS == 3)
			size_t K_lim = iBody[ib]._Owner->K_lim;
#endif

			// For each marker
			for (size_t m = 0; m < iBody[ib].markers.size(); m++) {

				// Reset the values of interpolated velocity
				std::fill(iBody[ib].markers[m].interpVel.begin(), iBody[ib].markers[m].interpVel.end(), 0.0);

				// Loop over support nodes
				for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {

					// Only interpolate over data this rank actually owns at the moment
					if (rank == iBody[ib].markers[m].support_rank[i]) {

						// Loop over directions x y z
						for (int dir = 0; dir < L_DIMS; dir++) {

							// Read given velocity component from support node, multiply by delta function
							// for that support node and sum to get interpolated velocity.
#if (L_DIMS == 3)
							iBody[ib].markers[m].interpVel[dir] += iBody[ib]._Owner->u(
									iBody[ib].markers[m].supp_i[i],
									iBody[ib].markers[m].supp_j[i],
									iBody[ib].markers[m].supp_k[i],
									dir, M_lim, K_lim, L_DIMS) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#else
							iBody[ib].markers[m].interpVel[dir] += iBody[ib]._Owner->u(
									iBody[ib].markers[m].supp_i[i],
									iBody[ib].markers[m].supp_j[i],
									dir, M_lim, L_DIMS) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#endif
						}
					}
				}
			}
		}
	}



	// Pass the necessary values between ranks
#ifdef L_BUILD_FOR_MPI
	ibm_gatherOffRankVels(level);
#endif

	// Write out interpolate velocity
#ifdef L_IBM_DEBUG
	for (int ib = 0; ib < iBody.size(); ib++) {
		if (iBody[ib]._Owner->level == level)
			ibm_debug_interpVel(ib);
	}
#endif
}

// *****************************************************************************
/// \brief	Compute restorative force at each marker in a body.
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_computeForce(int level) {

	// Loop over markers
	for (int ib = 0; ib < iBody.size(); ib++) {

		// Only do if this body is on this grid level
		if (iBody[ib]._Owner->level == level) {
			for (int m = 0; m < iBody[ib].markers.size(); m++) {
				for (int dir = 0; dir < L_DIMS; dir++) {

					// Compute restorative force (in lattice units)
					iBody[ib].markers[m].force_xyz[dir] = (iBody[ib].markers[m].interpVel[dir] - iBody[ib].markers[m].markerVel[dir]) / 1.0;
				}
			}

		// Write out force on markers
#ifdef L_IBM_DEBUG
		ibm_debug_markerForce(ib);
#endif
		}
	}
}

// *****************************************************************************
/// \brief	Spread restorative force back onto marker support.
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_spread(int ib) {
//
//	int rank = GridUtils::safeGetRank();
//
//	// For each marker
//	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
//		// Loop over support nodes
//		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
//
//			// Get size of grid
//			size_t M_lim = iBody[ib]._Owner->M_lim;
//			size_t K_lim = iBody[ib]._Owner->K_lim;
//
//			for (size_t dir = 0; dir < L_DIMS; dir++) {
//				// Add contribution of current marker force to support node Cartesian force vector using delta values computed when support was computed
//				iBody[ib]._Owner->force_xyz(
//					iBody[ib].markers[m].supp_i[i],
//					iBody[ib].markers[m].supp_j[i],
//					iBody[ib].markers[m].supp_k[i],
//					dir, M_lim, K_lim, L_DIMS) +=
//					iBody[ib].markers[m].deltaval[i] *
//					iBody[ib].markers[m].force_xyz[dir] *
//					iBody[ib].markers[m].epsilon *
//					iBody[ib].spacing/iBody[ib]._Owner->dh;
//			}
//		}
//	}
//
//#ifdef L_IBM_DEBUG
//
//#ifdef L_BUILD_FOR_MPI
//	// Get MPI Manager Instance
//	MpiManager *mpim = MpiManager::getInstance();
//#endif
//
//	// DEBUG -- write out res vector
//	std::ofstream testout;
//	testout.open(GridUtils::path_str + "/force_xyz_supp" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out", std::ios::app);
//	testout << "\nNEW TIME STEP" << std::endl;
//	// Get size of grid
//	size_t M_lim = iBody[ib]._Owner->M_lim;
//	size_t K_lim = iBody[ib]._Owner->K_lim;
//	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
//		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
//			testout << iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 0, M_lim, K_lim, L_DIMS) << "\t"
//					<< iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 1, M_lim, K_lim, L_DIMS) << std::endl;
//		}
//		testout << std::endl;
//	}
//	testout.close();
//#endif

}

// *****************************************************************************
/// \brief	Compute epsilon for a given iBody.
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_findEpsilon() {

	// Dummy iBody vector for gathering information into for epsilon calculation
	std::vector<IBBody> iBodyTemp;

	// Build dummy iBody vector
#ifndef L_BUILD_FOR_MPI

	// TODO Need to think of a better way to do this without copying iBody -> use pointers
	// Populate dummy iBody vector
	for (int ib = 0; ib < iBody.size(); ib++)
		iBodyTemp.push_back(iBody[ib]);

#else

	// Get MPI manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Pass marker positions and support site information and fill recvBuffer
	std::vector<std::vector<double>> recvBuffer;
	mpim->mpi_epsilonCommGather(recvBuffer);

	// Now build temporary iBody objects with custom constructor
	for (int ib = 0; ib < iBody.size(); ib++) {
		if (mpim->my_rank == iBody[ib].owningRank)
			iBodyTemp.emplace_back(iBody[ib], recvBuffer);
	}
#endif

	// Size epsilon values
	std::vector<std::vector<double>> epsilon(bodyIDToIdx.size(), std::vector<double>(0));

	// Loop through all iBodys this rank owns
	for (int ib = 0; ib < iBodyTemp.size(); ib++) {

		/* The Reproducing Kernel Particle Method (see Pinelli et al. 2010, JCP) requires suitable weighting
		to be computed to ensure conservation while using the interpolation functions. Epsilon is this weighting.
		We can use built-in libraries to solve the ensuing linear system in future. */

		// Declarations
		double Delta_I, Delta_J;

		//////////////////////////////////
		//	Build coefficient matrix A	//
		//		with a_ij values.		//
		//////////////////////////////////

		// Initialise 2D std vector with zeros
		std::vector< std::vector<double> > A(iBodyTemp[ib].markers.size(), std::vector<double>(iBodyTemp[ib].markers.size(), 0.0));


		// Loop over support of marker I and integrate delta value multiplied by delta value of marker J.
		for (size_t I = 0; I < iBodyTemp[ib].markers.size(); I++) {

			// Loop over markers J
			for (size_t J = 0; J < iBodyTemp[ib].markers.size(); J++) {

				// Sum delta values evaluated for each support of I
				for (size_t s = 0; s < iBodyTemp[ib].markers[I].deltaval.size(); s++) {

					Delta_I = iBodyTemp[ib].markers[I].deltaval[s];
#if (L_DIMS == 3)
					Delta_J =
						ibm_deltaKernel(
						(iBodyTemp[ib].markers[J].position[eXDirection] - iBodyTemp[ib].markers[I].supp_x[s]) / iBodyTemp[ib]._Owner->dh,
						iBodyTemp[ib].markers[J].dilation
						) *
						ibm_deltaKernel(
						(iBodyTemp[ib].markers[J].position[eYDirection] - iBodyTemp[ib].markers[I].supp_y[s]) / iBodyTemp[ib]._Owner->dh,
						iBodyTemp[ib].markers[J].dilation
						) *
						ibm_deltaKernel(
						(iBodyTemp[ib].markers[J].position[eZDirection] - iBodyTemp[ib].markers[I].supp_z[s]) / iBodyTemp[ib]._Owner->dh,
						iBodyTemp[ib].markers[J].dilation
						);
#else
					Delta_J =
						ibm_deltaKernel(
						(iBodyTemp[ib].markers[J].position[eXDirection] - iBodyTemp[ib].markers[I].supp_x[s]) / iBodyTemp[ib]._Owner->dh,
						iBodyTemp[ib].markers[J].dilation
						) *
						ibm_deltaKernel(
						(iBodyTemp[ib].markers[J].position[eYDirection] - iBodyTemp[ib].markers[I].supp_y[s]) / iBodyTemp[ib]._Owner->dh,
						iBodyTemp[ib].markers[J].dilation
						);
#endif
					// Multiply by local area (or volume in 3D)
					A[I][J] += Delta_I * Delta_J * iBodyTemp[ib].markers[I].local_area;
				}

				// Multiply by arc length between markers in lattice units
				A[I][J] = A[I][J] * (iBodyTemp[ib].spacing / iBodyTemp[ib]._Owner->dh);
			}
		}

		// Create vectors
		epsilon[iBodyTemp[ib].id].resize(iBodyTemp[ib].markers.size(), 0.0);
		std::vector<double> bVector (iBodyTemp[ib].markers.size(), 1.0);

		//////////////////
		// Solve system //
		//////////////////

		// Settings
		double tolerance = 1.0e-5;
		int maxiterations = 2500;

		// Biconjugate gradient stabilised method for solving asymmetric linear systems
		double minimum_residual_achieved = ibm_bicgstab(A, bVector, epsilon[iBodyTemp[ib].id], tolerance, maxiterations);
	}

	// Distribute the epsilon values back to the iBody
#ifndef L_BUILD_FOR_MPI

	// Loop through all iBodies
	for (int ib = 0; ib < iBody.size(); ib++) {
		for (int m = 0; m < iBody[ib].markers.size(); m++) {
			iBody[ib].markers[m].epsilon = epsilon[ib][m];
		}
	}

#else

	// Perform MPI communication and insert correct epsilon values
	mpim->mpi_epsilonCommScatter(epsilon);
#endif
}

// *****************************************************************************
/// \brief	Biconjugate gradient method.
///
///			Biconjugate gradient stabilised method of solving a linear system 
///			Ax = b. Solution is performed iteratively.
///
/// \param	Amatrix			the A matrix in the linear system.
/// \param	bVector			the b vector in the linear system.
/// \param	epsilon			epsilon paramters for each marker.
/// \param	tolerance		tolerance of solution.
/// \param	maxiterations	maximum number of iterations.
/// \returns the minimum residual achieved by the solver.
double ObjectManager::ibm_bicgstab(
	std::vector< std::vector<double> >& Amatrix, std::vector<double>& bVector,
	std::vector<double>& epsilon,
	double tolerance, int maxiterations) {

	// Declarations //

	// Scalars
    double bic_alpha, bic_omega, bic_beta, res_current;
	double res_min = 100.0;		// Arbitrary selection of initial minimum residual -- deliberately big so that first loop is going to generate a epsilon vector with residual better than this.
	// Number of markers
    size_t nls = epsilon.size();
	// Vectors
    std::vector<double> bic_rho (2, 0.0); // Need both i and i-1 instances at same time so need to declare as a 2x1 vector
	std::vector<double> bic_s (nls, 0.0);
	std::vector<double> bic_t (nls, 0.0);
	std::vector<double> epsilon_best (nls, 0.0);
	std::vector<double> bic_r, bic_v, bic_p, bic_rhat;

	// Step 1: Use initial guess to compute r vector
	std::vector<double> bic_Ax = GridUtils::matrix_multiply(Amatrix,epsilon);
	for (size_t i = 0; i < nls; i++) {
		bic_r.push_back(bVector[i] - bic_Ax[i]);
	}

	// Step 2: Choose arbitrary vector r_hat
	bic_rhat = bic_r;

	// Step 3: rho0 = alpha = omega0 = 1
	bic_alpha = 1.0;
	bic_omega = bic_alpha;
    bic_rho[0] = bic_alpha;

	// Step 4: v0 = p0 = 0
	for (size_t i = 0; i < nls; i++) {
        bic_v.push_back(0.0);
		bic_p.push_back(0.0);
    }

	// Step 5: Iterate
	for (int i = 1; i < maxiterations; i++) {

		// Step 5a: Compute new rho
		bic_rho[1] = GridUtils::dotprod(bic_rhat, bic_r);

		// Step 5b: Compute beta
		bic_beta = (bic_rho[1] / bic_rho[0]) * (bic_alpha / bic_omega);

		// Step 5c: Compute new p vector
		for (size_t j = 0; j < nls; j++) {
			bic_p[j] = bic_r[j] + bic_beta * ( bic_p[j] - bic_omega * bic_v[j] );
		}

		// Step 5d: Compute new v vector
		bic_v = GridUtils::matrix_multiply(Amatrix,bic_p);

		// Step 5e: Compute alpha
		bic_alpha = bic_rho[1] / GridUtils::dotprod(bic_rhat, bic_v);
		// bic_rho is not used again so copy last element back ready for next iteration
		bic_rho[0] = bic_rho[1];

		// Step 5f: Compute new s vector
		for (size_t j = 0; j < nls; j++) {
			bic_s[j] = bic_r[j] - (bic_alpha * bic_v[j]);
		}

		// Step 5g: Compute t vector
		bic_t = GridUtils::matrix_multiply(Amatrix,bic_s);

		// Step 5h: Compute new omega
		bic_omega = GridUtils::dotprod(bic_t, bic_s) / GridUtils::dotprod(bic_t, bic_t);

		// Step 5i: Update epsilon
		for (size_t j = 0; j < nls; j++) {
			epsilon[j] += bic_alpha * bic_p[j] + bic_omega * bic_s[j];
		}

		// Step 5j: Compute residual
		for (size_t j = 0; j < nls; j++) {
			bic_r[j] = bic_s[j] - bic_omega * bic_t[j];
		}

		// Step 5k: Check residual and store epsilon if best so far
		res_current = sqrt(GridUtils::dotprod(bic_r, bic_r));
		if ( res_current < res_min) {
			res_min = res_current;	// Note best residual
			for (size_t j = 0; j < nls; j++) {
				epsilon_best[j] = epsilon[j];
			}
			// Check tolerance of best residual
			if ( res_min <= tolerance ) {
				// Before exiting, update epsilon to best epsilon
				for (size_t j = 0; j < nls; j++) {
					epsilon[j] = epsilon_best[j];
				}
				return res_min;
			}
		}

		if (i == maxiterations-1) {
			// Warn that max iterations hit
			*GridUtils::logfile << "Max iterations hit -- values of epsilon might not be converged. Try adjusting the number of Lagrange markers or the grid resolution to adjust support overlap. Setting Epsilon to 2." << std::endl;
			for (size_t j = 0; j < nls; j++) {
				epsilon[j] = 2;
			}
		}
	}


	// Before exiting, update epsilon to best epsilon
	for (size_t j = 0; j < nls; j++) {
		epsilon[j] = epsilon_best[j];
	}
	return res_min;
}


// *****************************************************************************
/// \brief	Write out marker positions for debugging.
///
///			Writes out the marker positions for debugging
///
/// \param	ib			id of body being written out.
///
void ObjectManager::ibm_debug_markerPosition(int ib) {

	// Only do if there is data to write
	if (iBody[ib].markers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream bodyout;
		bodyout.open(GridUtils::path_str + "/IBbody_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out");
		bodyout << "Marker\tx\ty\tz" << std::endl;

		// Loop through markers
		for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
			bodyout << iBody[ib].markers[i].id << "\t" << iBody[ib].markers[i].position[eXDirection] << "\t" << iBody[ib].markers[i].position[eYDirection] << "\t" << iBody[ib].markers[i].position[eZDirection] << std::endl;
		}
		bodyout.close();
	}
}

// *****************************************************************************
/// \brief	Write out support site information for debugging.
///
///			Writes out support site information for debugging
///
/// \param	ib			id of body being written out.
/// \param	m			id of marker being written out.
void ObjectManager::ibm_debug_supportInfo(int ib, int m, int s) {

	// Only do if there is data to write
	if (iBody[ib].markers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream supportout;
		supportout.open(GridUtils::path_str + "/IBsupport_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out", std::ios::app);

		// Write out the first (nearest) support marker
		if (m == 0 && s == 0)
			supportout << "Marker\tRank\ti\tj\tk\tX\tY\tZ\tdeltaVal" << std::endl;

		// Write out info
		supportout
			<< iBody[ib].markers[m].id << "\t"
			<< iBody[ib].markers[m].support_rank.back() << "\t"
			<< iBody[ib].markers[m].supp_i.back() << "\t"
			<< iBody[ib].markers[m].supp_j.back() << "\t"
			<< iBody[ib].markers[m].supp_k.back() << "\t"
			<< iBody[ib].markers[m].supp_x.back() << "\t"
			<< iBody[ib].markers[m].supp_y.back() << "\t"
			<< iBody[ib].markers[m].supp_z.back() << "\t"
			<< iBody[ib].markers[m].deltaval.back() << std::endl;

		supportout.close();
	}
}

// *****************************************************************************
/// \brief	Write out epsilon for debugging.
///
///			Writes out the epsilon values for debugging
///
/// \param	ib			id of body being written out.
///
void ObjectManager::ibm_debug_epsilon(int ib) {

	// Only do if there is data to write
	if (iBody[ib].markers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream epout;
		epout.open(GridUtils::path_str + "/Epsilon_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
		epout << "NEW TIME STEP" << std::endl;
		epout << "Marker\tEpsilon" << std::endl;

		// Loop through markers
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			epout << iBody[ib].markers[m].id << "\t" << iBody[ib].markers[m].epsilon << std::endl;
		}
		epout << std::endl;
		epout.close();
	}
}




// *****************************************************************************
/// \brief	Write out interpolated velocity for debugging.
///
///			Writes out the interpolated velocity values for debugging
///
/// \param	ib			id of body being written out.
///
void ObjectManager::ibm_debug_interpVel(int ib) {

	// Only do if there is data to write
	if (iBody[ib].markers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream predout;
		predout.open(GridUtils::path_str + "/interpVel_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
		predout << "NEW TIME STEP" << std::endl;
		predout << "Marker\tVelX\tVelY\tVelZ" << std::endl;

		// Loop through markers
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			predout << iBody[ib].markers[m].id << "\t" << iBody[ib].markers[m].interpVel[eXDirection] << "\t" << iBody[ib].markers[m].interpVel[eYDirection] << "\t" << iBody[ib].markers[m].interpVel[eZDirection] << std::endl;
		}
		predout << std::endl;
		predout.close();
	}
}

// *****************************************************************************
/// \brief	Write out force on markers for debugging.
///
///			Writes out force on markers for debugging
///
/// \param	ib			id of body being written out.
///
void ObjectManager::ibm_debug_markerForce(int ib) {

	// Only do if there is data to write
	if (iBody[ib].markers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream forceout;
		forceout.open(GridUtils::path_str + "/force_xyz_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
		forceout << "NEW TIME STEP" << std::endl;
		forceout << "Marker\tFx\tFy\tFz" << std::endl;

		// Loop through markers
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			forceout << iBody[ib].markers[m].id << "\t" << iBody[ib].markers[m].force_xyz[eXDirection] << "\t" << iBody[ib].markers[m].force_xyz[eYDirection] << "\t" << iBody[ib].markers[m].force_xyz[eZDirection] << std::endl;
		}
		forceout << std::endl;
		forceout.close();
	}
}







// *****************************************************************************
/// \brief	Moves iBodies after applying IBM.
///
///			Wrapper for relocating markers of an iBody be calling appropriate
///			positional update routine.
///
void ObjectManager::ibm_debug_supportVel(int ib) {

	// Get rank safely
	int rank = GridUtils::safeGetRank();

	// Open file and write header
	std::ofstream testout;
	testout.open(GridUtils::path_str + "/velSupp" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out", std::ios::app);
	testout << "\nNEW TIME STEP" << std::endl;

	// Loop through markers and support sites
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
			testout << iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], 0, iBody[ib]._Owner->M_lim, L_DIMS) << "\t"
					<< iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], 1, iBody[ib]._Owner->M_lim, L_DIMS) << std::endl;
		}
		testout << std::endl;
	}
	testout.close();
}


// *****************************************************************************
/// \brief	Moves iBodies after applying IBM.
///
///			Wrapper for relocating markers of an iBody be calling appropriate
///			positional update routine.
///
void ObjectManager::ibm_debug_supportForce(int ib) {

	// Get rank safely
	int rank = GridUtils::safeGetRank();

	// Open file and write header
	std::ofstream testout;
	testout.open(GridUtils::path_str + "/force_xyz_supp" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out", std::ios::app);
	testout << "\nNEW TIME STEP" << std::endl;

	// Get size of grid
	size_t M_lim = iBody[ib]._Owner->M_lim;
	size_t K_lim = iBody[ib]._Owner->K_lim;

	// Loop through markers and support sites
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
			testout << iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 0, M_lim, K_lim, L_DIMS) << "\t"
					<< iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 1, M_lim, K_lim, L_DIMS) << std::endl;
		}
		testout << std::endl;
	}
	testout.close();
}

// ****************************************************************************
