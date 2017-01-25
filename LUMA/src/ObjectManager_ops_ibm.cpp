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
/// Perform IBM procedure.
void ObjectManager::ibm_apply() {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	// Loop over array of IB_bodies and perform IB operations
	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {

#ifdef L_IBM_DEBUG
		// DEBUG -- write out support coordinates
		std::ofstream suppout;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			suppout.open(GridUtils::path_str + "/Supp_" + std::to_string(ib) + "_" + std::to_string(m) + "_rank" + std::to_string(mpim->my_rank) + ".out",std::ios::app);
			suppout << "\nNEW TIME STEP" << std::endl;
			suppout << "x\ty\tz" << std::endl;
			for (size_t i = 0; i < iBody[ib].markers[m].supp_i.size(); i++) {
				if (L_DIMS == 3) {
					suppout << iBody[ib]._Owner->XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << iBody[ib]._Owner->YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << iBody[ib]._Owner->ZPos[iBody[ib].markers[m].supp_k[i]] << std::endl;
				} else {
					suppout << iBody[ib]._Owner->XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << iBody[ib]._Owner->YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << 0.0 << std::endl;
				}
			}
			suppout.close();
		}
#endif

#ifdef L_IBM_DEBUG
		// DEBUG -- write out epsilon values
		std::ofstream epout;
		epout.open(GridUtils::path_str + "/Epsilon_" + std::to_string(ib) + "_rank" + std::to_string(mpim->my_rank) + ".out",std::ios::app);
		epout << "\nNEW TIME STEP" << std::endl;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			epout << iBody[ib].markers[m].epsilon << std::endl;
		}
		epout.close();
#endif

			// Interpolate velocity
			ibm_interpol(ib);

#ifdef L_IBM_DEBUG
		// DEBUG -- write out interpolate velocity values
		std::ofstream predout;
		predout.open(GridUtils::path_str + "/interpVel_" + std::to_string(ib) + "_rank" + std::to_string(mpim->my_rank) + ".out",std::ios::app);
		predout << "\nNEW TIME STEP" << std::endl;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			predout << iBody[ib].markers[m].fluid_vel[0] << "\t" << iBody[ib].markers[m].fluid_vel[1] << std::endl;
		}
		predout.close();
#endif

			// Compute restorative force
			ibm_computeForce(ib);

#ifdef L_IBM_DEBUG
		// DEBUG -- write out Lagrange force values
		std::ofstream forceout;
		forceout.open(GridUtils::path_str + "/force_xyz_" + std::to_string(ib) + "_rank" + std::to_string(mpim->my_rank) + ".out",std::ios::app);
		forceout << "\nNEW TIME STEP" << std::endl;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			forceout << iBody[ib].markers[m].force_xyz[0] << "\t" << iBody[ib].markers[m].force_xyz[1] << std::endl;
		}
		forceout.close();
#endif

			// Spread force back to lattice (Cartesian vector)
			ibm_spread(ib);

		}
}

// *****************************************************************************
/// \brief	Moves iBodies after applying IBM.
///
///			Wrapper for relocating markers of an iBody be calling appropriate
///			positional update routine.
///
void ObjectManager::ibm_moveBodies() {

	// Loop over bodies launching positional update if deformable to compute new locations of markers
	*GridUtils::logfile << "Relocating markers as required..." << std::endl;
	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {

		// If body is deformable it needs a positional update
		if (iBody[ib].deformable) {

			// Call structural or forced positional update and recompute support
			ibm_positionUpdate(ib);

#ifndef L_STOP_EPSILON_RECOMPUTE
			// Recompute epsilon
			ibm_findEpsilon(ib);
#endif

		}
	}

#if defined L_INSERT_FILARRAY
	// Special bit for filament-based plates where flexible centreline is used to update position of others in group
	*GridUtils::logfile << "Filament-based plate positional update..." << std::endl;
	ibm_positionUpdateGroup(999);
#endif

}


// *****************************************************************************
/// \brief	Initialise the array of iBodies.
///
///			Computes support and epsilon values.
///
void ObjectManager::ibm_initialise() {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	// Loop over the number of bodies in the iBody array
	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {

#ifdef L_IBM_DEBUG
		// DEBUG -- write out marker coordinates
		std::ofstream bodyout;
		bodyout.open(GridUtils::path_str + "/IBbody_" + std::to_string(ib) + "_rank" + std::to_string(mpim->my_rank) + ".out");
		bodyout << "x\ty\tz" << std::endl;
		for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
			bodyout << iBody[ib].markers[i].position[0] << "\t" << iBody[ib].markers[i].position[1] << "\t" << iBody[ib].markers[i].position[2] << std::endl;
		}
		bodyout.close();
#endif

		// Compute extended support for each marker in the body (or body portion)
		for (int m = 0; m < static_cast<int>(iBody[ib].markers.size()); m++)
		{
			ibm_findSupport(ib, m);		// Pass body ID and marker ID
		}

		// Find epsilon for the body
		ibm_findEpsilon(ib);

	}
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

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	// Declarations
	int inear, jnear;							// Nearest node indices (global indices)
	std::vector<double> nearpos;				// Position of the nearest marker
	int knear;
	int s = 0;										// The sth support marker which has been found

#ifdef L_BUILD_FOR_MPI
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
	inear = nearijk[eXDirection];
	jnear = nearijk[eYDirection];
	knear = nearijk[eZDirection];

	// Get position of nearest node
	nearpos.push_back(iBody[ib]._Owner->XPos[nearijk[eXDirection]]);
	nearpos.push_back(iBody[ib]._Owner->YPos[nearijk[eYDirection]]);
	nearpos.push_back(iBody[ib]._Owner->ZPos[nearijk[eZDirection]]);

	// Side length of support region defined as 3 x dilation parameter
	iBody[ib].markers[m].dilation = 1.0;

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

#ifdef L_IBM_DEBUG
	// DEBUG -- Write out support as they are found
	std::ofstream supportout;
	supportout.open(GridUtils::path_str
		+ "/IBsupport_" + std::to_string(ib)
		+ "_rank" + std::to_string(mpim->my_rank)
		+ ".out", std::ios::app);

		// Write out the first (nearest) support marker
	if (m == 0)
		supportout << "Marker\tRank\tX\tY\tZ" << std::endl;
	supportout
		<< m << "\t"
		<< iBody[ib].markers[m].support_rank.back() << "\t"
		<< iBody[ib]._Owner->XPos[iBody[ib].markers[m].supp_i.back()] << "\t"
		<< iBody[ib]._Owner->YPos[iBody[ib].markers[m].supp_j.back()] << "\t"
		<< iBody[ib]._Owner->ZPos[iBody[ib].markers[m].supp_k.back()] << std::endl;

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

						// Add owning rank as this one for now
						iBody[ib].markers[m].support_rank.push_back(mpim->my_rank);


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
						if (owner_direction != -1)
						{
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

#ifdef L_IBM_DEBUG
						// DEBUG -- add this support point to the list
						supportout
							<< m << "\t"
							<< iBody[ib].markers[m].support_rank.back() << "\t"
							<< estimated_position[eXDirection] << "\t"
							<< estimated_position[eYDirection] << "\t"
							<< estimated_position[eZDirection] << std::endl;
#endif
					}
				}
			}
		}
	}

#ifdef L_IBM_DEBUG
	supportout.close();
#endif

	// Store normalised area of support region
	iBody[ib].markers[m].local_area = 1; // Area remains constant at the local grid level
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
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_interpol(int ib) {

	// Get grid sizes
	size_t M_lim = iBody[ib]._Owner->M_lim;
#if (L_DIMS == 3)
	size_t K_lim = iBody[ib]._Owner->K_lim;
#endif


	// For each marker
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {

		// Reset the values of interpolated velocity
		std::fill(iBody[ib].markers[m].fluid_vel.begin(), iBody[ib].markers[m].fluid_vel.end(), 0.0);

		// Loop over support nodes
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {

			// Loop over directions x y z
			for (int dir = 0; dir < L_DIMS; dir++) {

				// Read given velocity component from support node, multiply by delta function
				// for that support node and sum to get interpolated velocity.

#if (L_DIMS == 3)
			iBody[ib].markers[m].fluid_vel[dir] += 
				iBody[ib]._Owner->u(
				iBody[ib].markers[m].supp_i[i],
				iBody[ib].markers[m].supp_j[i],
				iBody[ib].markers[m].supp_k[i],
				dir,
				M_lim, K_lim, L_DIMS
				) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#else
			iBody[ib].markers[m].fluid_vel[dir] += 
				iBody[ib]._Owner->u(
				iBody[ib].markers[m].supp_i[i],
				iBody[ib].markers[m].supp_j[i],
				dir,
				M_lim, L_DIMS
				) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#endif
			}
		}
	}

#ifdef L_IBM_DEBUG
	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	// DEBUG -- write out res vector
	std::ofstream testout;
	testout.open(GridUtils::path_str + "/velSupp" + std::to_string(ib) + "_rank" + std::to_string(mpim->my_rank) + ".out", std::ios::app);
	testout << "\nNEW TIME STEP" << std::endl;
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
			testout << iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], 0, M_lim, L_DIMS) << "\t"
					<< iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], 1, M_lim, L_DIMS) << std::endl;
		}
		testout << std::endl;
	}
	testout.close();
#endif

}

// *****************************************************************************
/// \brief	Compute restorative force at each marker in a body.
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_computeForce(int ib) {

	// Loop over markers
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		for (int dir = 0; dir < L_DIMS; dir++) {
			// Compute restorative force (in lattice units)
			iBody[ib].markers[m].force_xyz[dir] = 
				(iBody[ib].markers[m].desired_vel[dir] - iBody[ib].markers[m].fluid_vel[dir]) /	1.0;	// Time step in grid-normalised lattice units
		}
	}
}

// *****************************************************************************
/// \brief	Spread restorative force back onto marker support.
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_spread(int ib) {

	// For each marker
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		// Loop over support nodes
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {

			// Get size of grid
			size_t M_lim = iBody[ib]._Owner->M_lim;
			size_t K_lim = iBody[ib]._Owner->K_lim;

			for (size_t dir = 0; dir < L_DIMS; dir++) {
				// Add contribution of current marker force to support node Cartesian force vector using delta values computed when support was computed
				iBody[ib]._Owner->force_xyz(
					iBody[ib].markers[m].supp_i[i],
					iBody[ib].markers[m].supp_j[i],
					iBody[ib].markers[m].supp_k[i],
					dir, M_lim, K_lim, L_DIMS) +=
					iBody[ib].markers[m].deltaval[i] * 
					iBody[ib].markers[m].force_xyz[dir] * 
					iBody[ib].markers[m].epsilon * 
					iBody[ib].spacing/iBody[ib]._Owner->dh;
			}
		}
	}

#ifdef L_IBM_DEBUG

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	// DEBUG -- write out res vector
	std::ofstream testout;
	testout.open(GridUtils::path_str + "/force_xyz_supp" + std::to_string(ib) + "_rank" + std::to_string(mpim->my_rank) + ".out", std::ios::app);
	testout << "\nNEW TIME STEP" << std::endl;
	// Get size of grid
	size_t M_lim = iBody[ib]._Owner->M_lim;
	size_t K_lim = iBody[ib]._Owner->K_lim;
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
			testout << iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 0, M_lim, K_lim, L_DIMS) << "\t"
					<< iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 1, M_lim, K_lim, L_DIMS) << std::endl;
		}
		testout << std::endl;
	}
	testout.close();
#endif

}

// *****************************************************************************
/// \brief	Compute epsilon for a given iBody.
/// \param	ib	iBody being operated on.
double ObjectManager::ibm_findEpsilon(int ib) {

	/***** TODO: In order to compute epsilon using this linear system a single process needs 
	to know all the delta values of all the support points in the body to build the whole system. 
	We could simply nominate a rank to do this based on the body ID but we will have to think about
	how to gather the data -- this is one of the problems of not using the other approach I guess. *****/

#ifdef L_BUILD_FOR_MPI

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	// Create the vector of structs that needs to be sent
	std::vector<IBInfo> msg_deltas;
	for (int m = 0; m < iBody[ib].markers.size(); m++)
		msg_deltas.emplace_back(&iBody[ib], m, eIBDeltaSum);


	if (mpim->my_rank == 1) {
		std::cout << "Rank " << mpim->my_rank << std::endl;
		std::cout << "Size of vector is " << msg_deltas.size() << std::endl;
		for (int i = 0; i < msg_deltas.size(); i++) {
			std::cout << msg_deltas[i].markerX << " " << msg_deltas[i].markerY << " " <<  msg_deltas[i].markerZ << std::endl;
		}
	}

	// Map the IBInfo to an MPI_Type_struct object (just use the first element of vector to get the mapping)
	MPI_Datatype mpi_struct_type;
	msg_deltas[0].mapToMpiStruct(eIBDeltaSum, &mpi_struct_type);


	if (mpim->my_rank == 1) {
		int vecSizeSend = msg_deltas.size();
		MPI_Send(&vecSizeSend, 1, MPI_INT, 0, 99, MPI_COMM_WORLD);

		MPI_Send(&msg_deltas.front(), vecSizeSend, mpi_struct_type, 0, 98, MPI_COMM_WORLD);
	}

	std::vector<IBInfo> msg_deltasRecv;
	if (mpim->my_rank == 0) {
		int vecSizeRecv;
		MPI_Recv(&vecSizeRecv, 1, MPI_INT, 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		msg_deltasRecv.resize(vecSizeRecv);
//		std::cout << "Size of vector is " << vecSizeRecv << std::endl;
//		std::cout << "Size of vector is " << msg_deltasRecv.size() << std::endl;
		MPI_Recv(&msg_deltasRecv.front(), vecSizeRecv, mpi_struct_type, 1, 98, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		std::cout << "Rank " << mpim->my_rank << std::endl;
		std::cout << "Size of vector is " << msg_deltasRecv.size() << std::endl;
		for (int i = 0; i < msg_deltasRecv.size(); i++) {
			std::cout << msg_deltasRecv[i].markerX << " " << msg_deltasRecv[i].markerY << " " << msg_deltasRecv[i].markerZ << std::endl;
		}
	}


	// Send to managing rank


	// If managing rank then unpack data and build A matrix


	// Solve for epsilon


	// Create new message with all epsilon values for every point
	IBInfo *msg_epsilon = new IBInfo();

	// Send messages to all ranks


	// Each rank receive return message and copy epsilon value to matching markers


	// Destroy objects once we have finished
	//delete msg_deltas;
	delete msg_epsilon;

#endif


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
	std::vector< std::vector<double> > A (iBody[ib].markers.size(), std::vector<double>(iBody[ib].markers.size(), 0.0) );


	// Loop over support of marker I and integrate delta value multiplied by delta value of marker J.
	for (size_t I = 0; I < iBody[ib].markers.size(); I++) {

		// Loop over markers J
		for (size_t J = 0; J < iBody[ib].markers.size(); J++) {

			// Sum delta values evaluated for each support of I
			for (size_t s = 0; s < iBody[ib].markers[I].supp_i.size(); s++) {

				Delta_I = iBody[ib].markers[I].deltaval[s];
#if (L_DIMS == 3)
				Delta_J =
					ibm_deltaKernel(
					(iBody[ib].markers[J].position[0] - iBody[ib]._Owner->XPos[iBody[ib].markers[I].supp_i[s]]) / iBody[ib]._Owner->dh, 
					iBody[ib].markers[J].dilation
					) *
					ibm_deltaKernel(
					(iBody[ib].markers[J].position[1] - iBody[ib]._Owner->YPos[iBody[ib].markers[I].supp_j[s]]) / iBody[ib]._Owner->dh,
					iBody[ib].markers[J].dilation
					) *
					ibm_deltaKernel(
					(iBody[ib].markers[J].position[2] - iBody[ib]._Owner->ZPos[iBody[ib].markers[I].supp_k[s]]) / iBody[ib]._Owner->dh,
					iBody[ib].markers[J].dilation
					);
#else
				Delta_J =
					ibm_deltaKernel(
					(iBody[ib].markers[J].position[0] - iBody[ib]._Owner->XPos[iBody[ib].markers[I].supp_i[s]]) / iBody[ib]._Owner->dh,
					iBody[ib].markers[J].dilation
					) *
					ibm_deltaKernel(
					(iBody[ib].markers[J].position[1] - iBody[ib]._Owner->YPos[iBody[ib].markers[I].supp_j[s]]) / iBody[ib]._Owner->dh,
					iBody[ib].markers[J].dilation
					);
#endif
				// Multiply by local area (or volume in 3D)
				A[I][J] += Delta_I * Delta_J * iBody[ib].markers[I].local_area;
			}

			// Multiply by arc length between markers in lattice units
			A[I][J] = A[I][J] * (iBody[ib].spacing / iBody[ib]._Owner->dh);

		}

	}


#ifdef L_IBM_DEBUG

	// DEBUG -- write out A
	std::ofstream Aout;
	Aout.open(GridUtils::path_str + "/Amatrix_" + std::to_string(ib) + "_rank" + std::to_string(mpim->my_rank) + ".out");
	for (size_t i = 0; i < A.size(); i++) {
		Aout << "\n";
		for (size_t j = 0; j < A.size(); j++) {
			Aout << A[i][j] << "\t";
		}
	}
	Aout.close();
#endif


	// Create vectors
	std::vector<double> epsilon (iBody[ib].markers.size(), 0.0);
	std::vector<double> bVector (iBody[ib].markers.size(), 1.0);


	//////////////////
	// Solve system //
	//////////////////

	// Settings
    double tolerance = 1.0e-5;
	int maxiterations = 2500;
	double minimum_residual_achieved;

    // Biconjugate gradient stabilised method for solving asymmetric linear systems
    minimum_residual_achieved = ibm_bicgstab(A, bVector, epsilon, tolerance, maxiterations);

	// Now assign epsilon to the markers
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {

		iBody[ib].markers[m].epsilon = epsilon[m];
	}

	return minimum_residual_achieved;

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

// ****************************************************************************
