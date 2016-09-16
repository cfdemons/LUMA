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

	// Loop over array of IB_bodies and perform IB operations
	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {

#ifdef L_IBM_DEBUG
		// DEBUG -- write out support coordinates
		std::ofstream suppout;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			suppout.open(GridUtils::path_str + "/Supp_" + std::to_string(ib) + "_" + std::to_string(m) + "_rank" + std::to_string(MpiManager::my_rank) + ".out",std::ios::app);
			suppout << "\nNEW TIME STEP" << std::endl;
			suppout << "x\ty\tz" << std::endl;
			for (size_t i = 0; i < iBody[ib].markers[m].supp_i.size(); i++) {
				if (L_dims == 3) {
					suppout << g.XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << g.YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << g.ZPos[iBody[ib].markers[m].supp_k[i]] << std::endl;
				} else {
					suppout << g.XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << g.YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << 0.0 << std::endl;
				}
			}
			suppout.close();
		}
#endif

#ifdef L_IBM_DEBUG
		// DEBUG -- write out epsilon values
		std::ofstream epout;
		epout.open(GridUtils::path_str + "/Epsilon_" + std::to_string(ib) + "_rank" + std::to_string(MpiManager::my_rank) + ".out",std::ios::app);
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
		predout.open(GridUtils::path_str + "/interpVel_" + std::to_string(ib) + "_rank" + std::to_string(MpiManager::my_rank) + ".out",std::ios::app);
		predout << "\nNEW TIME STEP" << std::endl;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			predout << iBody[ib].markers[m].fluid_vel[0] << "\t" << iBody[ib].markers[m].fluid_vel[1] << std::endl;
		}
		predout.close();
#endif

			// Compute restorative force
			ibm_computeforce(ib);

#ifdef L_IBM_DEBUG
		// DEBUG -- write out Lagrange force values
		std::ofstream forceout;
		forceout.open(GridUtils::path_str + "/force_xyz_" + std::to_string(ib) + "_rank" + std::to_string(MpiManager::my_rank) + ".out",std::ios::app);
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
void ObjectManager::ibm_move_bodies() {

	// Loop over bodies launching positional update if deformable to compute new locations of markers
	*GridUtils::logfile << "Relocating markers as required..." << std::endl;
	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {

		// If body is deformable it needs a positional update
		if (iBody[ib].deformable) {

			// Call structural or forced positional update and recompute support
			ibm_position_update(ib);

#ifndef L_STOP_EPSILON_RECOMPUTE
			// Recompute epsilon
			ibm_findepsilon(ib);
#endif

		}
	}

#if defined L_INSERT_FILARRAY
	// Special bit for filament-based plates where flexible centreline is used to update position of others in group
	*GridUtils::logfile << "Filament-based plate positional update..." << std::endl;
	ibm_position_update_grp(999);
#endif

}


// *****************************************************************************
/// \brief	Initialise the array of iBodies.
///
///			Computes support and epsilon values.
///
void ObjectManager::ibm_initialise() {

	// Loop over the number of bodies in the iBody array
	for (int ib = 0; ib < static_cast<int>(iBody.size()); ib++) {

#ifdef L_IBM_DEBUG
		// DEBUG -- write out marker coordinates
		std::ofstream bodyout;
		bodyout.open(GridUtils::path_str + "/IBbody_" + std::to_string(ib) + "_rank" + std::to_string(MpiManager::my_rank) + ".out");
		bodyout << "x\ty\tz" << std::endl;
		for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
			bodyout << iBody[ib].markers[i].position[0] << "\t" << iBody[ib].markers[i].position[1] << "\t" << iBody[ib].markers[i].position[2] << std::endl;
		}
		bodyout.close();
#endif

		// Compute support for each marker
		for (int m = 0; m < static_cast<int>(iBody[ib].markers.size()); m++) {
			ibm_findsupport(ib, m);	// Pass body ID and marker ID
		}

		// Find epsilon for the body
		ibm_findepsilon(ib);

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
double ObjectManager::ibm_deltakernel(double radius, double dilation) {

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
void ObjectManager::ibm_findsupport(int ib, int m) {

	// Declarations
	int inear, jnear;					// Nearest node indices
	double dist_x, dist_y, delta_x, delta_y;	// Distances and deltas
#if (L_dims == 3)
	// Extras for 3D
	double dist_z, delta_z;
	int knear;
#endif


	// Find closest support node (simulate std::round as not availble on MSVC2012)
	inear = (int)std::floor(
		(iBody[ib].markers[m].position[0] - (iBody[ib]._Owner->XPos[0] - iBody[ib]._Owner->dx / 2.0)) 
		/ iBody[ib]._Owner->dx
		);
	jnear = (int)std::floor(
		(iBody[ib].markers[m].position[1] - (iBody[ib]._Owner->YPos[0] - iBody[ib]._Owner->dy / 2.0)) 
		/ iBody[ib]._Owner->dy
		);

#if (L_dims == 3)
	knear = (int)std::floor(
		(iBody[ib].markers[m].position[2] - (iBody[ib]._Owner->ZPos[0] - iBody[ib]._Owner->dz / 2.0)) 
		/ iBody[ib]._Owner->dz
		);
#endif


	// Define limits of support region (only have to do one since lattice is uniformly spaced with dx = dy = dz)
	// Following protocol for arbitrary grid spacing each marker should have at least 3 support nodes in each direction
	double h_plus = std::max(
		std::abs((iBody[ib]._Owner->XPos[inear + 1] - iBody[ib]._Owner->XPos[inear]) / iBody[ib]._Owner->dx), 
		std::abs((iBody[ib]._Owner->XPos[inear] - iBody[ib]._Owner->XPos[inear - 1]) / iBody[ib]._Owner->dx)
		);
	double h_minus = std::min(
		std::abs((iBody[ib]._Owner->XPos[inear + 1] - iBody[ib]._Owner->XPos[inear]) / iBody[ib]._Owner->dx), 
		std::abs((iBody[ib]._Owner->XPos[inear] - iBody[ib]._Owner->XPos[inear - 1]) / iBody[ib]._Owner->dx)
		);

	// Side length of support region defined as 3 x dilation paramter which is found from:
	iBody[ib].markers[m].dilation = (5.0/6.0) * h_plus + (1.0/6.0) * h_minus;	// TODO Find out why this has such a drastic effect on everything
		//+ ( (1.0/9.0)); // * (1 / pow(2,iBody[ib]._Owner->level)) );	// This last term is a small fraction of the local grid spacing in lattice units


	// Test to see if required support nodes are available
#if (L_dims == 3)

	if ( inear - 5 < 0 || static_cast<size_t>(inear + 5) >= iBody[ib]._Owner->N_lim ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body " << std::to_string(ib) << " is too near the X boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		exit(LUMA_FAILED);

	} else if ( jnear - 5 < 0 || static_cast<size_t>(jnear + 5) >= iBody[ib]._Owner->M_lim ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body " << std::to_string(ib) << " is too near the Y boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		exit(LUMA_FAILED);

	} else if ( knear - 5 < 0 || static_cast<size_t>(knear + 5) >= iBody[ib]._Owner->K_lim ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body " << std::to_string(ib) << " is too near the Z boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		exit(LUMA_FAILED);
	}


	// Loop over surrounding 5 nodes to find support nodes
	for (int i = inear - 5; i <= inear + 5; i++) {
		for (int j = jnear - 5; j <= jnear + 5; j++) {
			for (int k = knear - 5; k <= knear + 5; k++) {

				// Find distance between Lagrange marker and possible support node and decide whether in cage or not
				if (
					( fabs(iBody[ib]._Owner->XPos[inear] - iBody[ib]._Owner->XPos[i])/iBody[ib]._Owner->dx < 1.5*iBody[ib].markers[m].dilation ) &&
					( fabs(iBody[ib]._Owner->YPos[jnear] - iBody[ib]._Owner->YPos[j])/iBody[ib]._Owner->dx < 1.5*iBody[ib].markers[m].dilation ) &&
					( fabs(iBody[ib]._Owner->ZPos[knear] - iBody[ib]._Owner->ZPos[k])/iBody[ib]._Owner->dx < 1.5*iBody[ib].markers[m].dilation )
					) {

						// Lies within support region so store information
						iBody[ib].markers[m].supp_i.push_back(i);
						iBody[ib].markers[m].supp_j.push_back(j);
						iBody[ib].markers[m].supp_k.push_back(k);

						// Store normalised area of support region (actually a volume) computed
						// from the local grid spacing in lattice units (dx = 1 / 2^level = 1 on L0)
						iBody[ib].markers[m].local_area = pow( 1 / pow(2,iBody[ib]._Owner->level) ,3) ;

						// Distance between Lagrange marker and support node in lattice units
						dist_x = (iBody[ib]._Owner->XPos[i]-iBody[ib].markers[m].position[0]) / iBody[ib]._Owner->dx;
						dist_y = (iBody[ib]._Owner->YPos[j]-iBody[ib].markers[m].position[1]) / iBody[ib]._Owner->dy;
						dist_z = (iBody[ib]._Owner->ZPos[k]-iBody[ib].markers[m].position[2]) / iBody[ib]._Owner->dz;

						// Store delta function value
						delta_x = ibm_deltakernel(dist_x, iBody[ib].markers[m].dilation);
						delta_y = ibm_deltakernel(dist_y, iBody[ib].markers[m].dilation);
						delta_z = ibm_deltakernel(dist_z, iBody[ib].markers[m].dilation);

						iBody[ib].markers[m].deltaval.push_back( delta_x * delta_y * delta_z );

				}
			}
		}
	}

#else

	// 2D check support region
	if ( inear - 5 < 0 || static_cast<size_t>(inear + 5) >= iBody[ib]._Owner->N_lim ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body " << std::to_string(ib) << " is too near the X boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		exit(LUMA_FAILED);

	} else if ( jnear - 5 < 0 || static_cast<size_t>(jnear + 5) >= iBody[ib]._Owner->M_lim ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body " << std::to_string(ib) << " is too near the Y boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		exit(LUMA_FAILED);

	}

	// 2D version to find support nodes
	for (int i = inear - 5; i <= inear + 5; i++) {
		for (int j = jnear - 5; j <= jnear + 5; j++) {
			int k = 0;

			// Find distance between Lagrange marker and possible support node and decide whether in cage or not
			if (	( fabs(iBody[ib].markers[m].position[0] - iBody[ib]._Owner->XPos[i])/iBody[ib]._Owner->dx < 1.5*iBody[ib].markers[m].dilation ) &&
					( fabs(iBody[ib].markers[m].position[1] - iBody[ib]._Owner->YPos[j])/iBody[ib]._Owner->dx < 1.5*iBody[ib].markers[m].dilation )
				) {

					// Lies within support region so store information
					iBody[ib].markers[m].supp_i.push_back(i);
					iBody[ib].markers[m].supp_j.push_back(j);
					iBody[ib].markers[m].supp_k.push_back(k);

					// Store normalised area of support region = dx^2 = (1/2^level) ^ 2.
					iBody[ib].markers[m].local_area = 1; // Area remains constant at the local grid level

					//  Distance between Lagrange marker and support node in lattice units
					dist_x = (iBody[ib]._Owner->XPos[i]-iBody[ib].markers[m].position[0]) / iBody[ib]._Owner->dx;
					dist_y = (iBody[ib]._Owner->YPos[j]-iBody[ib].markers[m].position[1]) / iBody[ib]._Owner->dy;

					// Store delta function value
					delta_x = ibm_deltakernel(dist_x, iBody[ib].markers[m].dilation);
					delta_y = ibm_deltakernel(dist_y, iBody[ib].markers[m].dilation);

					iBody[ib].markers[m].deltaval.push_back( delta_x * delta_y );

			}
		}
	}
#endif
}

// *****************************************************************************
/// \brief	Interpolate velocity field onto markers
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_interpol(int ib) {

	// Get grid sizes
	size_t M_lim = iBody[ib]._Owner->M_lim;
#if (L_dims == 3)
	size_t K_lim = iBody[ib]._Owner->K_lim;
#endif


	// For each marker
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {

		// Reset the values of interpolated velocity
		std::fill(iBody[ib].markers[m].fluid_vel.begin(), iBody[ib].markers[m].fluid_vel.end(), 0.0);

		// Loop over support nodes
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {

			// Loop over directions x y z
			for (int dir = 0; dir < L_dims; dir++) {

				// Read given velocity component from support node, multiply by delta function
				// for that support node and sum to get interpolated velocity.

#if (L_dims == 3)
			iBody[ib].markers[m].fluid_vel[dir] += iBody[ib]._Owner->u(	iBody[ib].markers[m].supp_i[i],
														iBody[ib].markers[m].supp_j[i],
														iBody[ib].markers[m].supp_k[i],
														dir,
														M_lim, K_lim, L_dims
														) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#else
			iBody[ib].markers[m].fluid_vel[dir] += iBody[ib]._Owner->u(	iBody[ib].markers[m].supp_i[i],
														iBody[ib].markers[m].supp_j[i],
														dir,
														M_lim, L_dims
														) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#endif
			}
		}
	}

#ifdef L_IBM_DEBUG
		// DEBUG -- write out res vector
		std::ofstream testout;
		testout.open(GridUtils::path_str + "/velSupp" + std::to_string(ib) + "_rank" + std::to_string(MpiManager::my_rank) + ".out", std::ios::app);
		testout << "\nNEW TIME STEP" << std::endl;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
				testout << iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], 0, M_lim, L_dims) << "\t"
									  << iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], 1, M_lim, L_dims) << std::endl;
			}
			testout << std::endl;
		}
		testout.close();
#endif

}

// *****************************************************************************
/// \brief	Compute restorative force at each marker in a body.
/// \param	ib	iBody being operated on.
void ObjectManager::ibm_computeforce(int ib) {

	// Loop over markers
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		for (int dir = 0; dir < L_dims; dir++) {
			// Compute restorative force (in lattice units)
			iBody[ib].markers[m].force_xyz[dir] = (iBody[ib].markers[m].desired_vel[dir] - iBody[ib].markers[m].fluid_vel[dir]) /
				1.0;	// Time step in grid-normalised lattice units
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

			for (size_t dir = 0; dir < L_dims; dir++) {
				// Add contribution of current marker force to support node Cartesian force vector using delta values computed when support was computed
				iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], dir, M_lim, K_lim, L_dims) +=
					iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].force_xyz[dir] * iBody[ib].markers[m].epsilon * iBody[ib].spacing/iBody[ib]._Owner->dx;
			}
		}
	}

#ifdef L_IBM_DEBUG
		// DEBUG -- write out res vector
		std::ofstream testout;
		testout.open(GridUtils::path_str + "/force_xyz_supp" + std::to_string(ib) + "_rank" + std::to_string(MpiManager::my_rank) + ".out", std::ios::app);
		testout << "\nNEW TIME STEP" << std::endl;
		// Get size of grid
		size_t M_lim = iBody[ib]._Owner->M_lim;
		size_t K_lim = iBody[ib]._Owner->K_lim;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
				testout << iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 0, M_lim, K_lim, L_dims) << "\t"
									  << iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], 1, M_lim, K_lim, L_dims) << std::endl;
			}
			testout << std::endl;
		}
		testout.close();
#endif

}

// *****************************************************************************
/// \brief	Compute epsilon for a given iBody.
/// \param	ib	iBody being operated on.
double ObjectManager::ibm_findepsilon(int ib) {

	/* The Reproducing Kernel Particle Method (see Pinelli et al. 2010, JCP) requires suitable weighting
	to be computed to ensure conservation while using the interpolation functions. Epsilon is this weighting.
	We can use built-in libraries to solve the ensuing linear system in future. */

	// Declarations
	double Delta_I, Delta_J;

	///////////////////////////////////
	//	Build coefficient matrix A	//
	//		with a_ij values.		//
	///////////////////////////////////

	// Initialise 2D std vector with zeros
	std::vector< std::vector<double> > A (iBody[ib].markers.size(), std::vector<double>(iBody[ib].markers.size(), 0.0) );


	// Loop over support of marker I and integrate delta value multiplied by delta value of marker J.
	for (size_t I = 0; I < iBody[ib].markers.size(); I++) {

		// Loop over markers J
		for (size_t J = 0; J < iBody[ib].markers.size(); J++) {

			// Sum delta values evaluated for each support of I
			for (size_t s = 0; s < iBody[ib].markers[I].supp_i.size(); s++) {

				Delta_I = iBody[ib].markers[I].deltaval[s];
#if (L_dims == 3)
				Delta_J =
					ibm_deltakernel(
					(iBody[ib].markers[J].position[0] - iBody[ib]._Owner->XPos[iBody[ib].markers[I].supp_i[s]]) / iBody[ib]._Owner->dx, 
					iBody[ib].markers[J].dilation
					) *
					ibm_deltakernel(
					(iBody[ib].markers[J].position[1] - iBody[ib]._Owner->YPos[iBody[ib].markers[I].supp_j[s]]) / iBody[ib]._Owner->dx,
					iBody[ib].markers[J].dilation
					) *
					ibm_deltakernel(
					(iBody[ib].markers[J].position[2] - iBody[ib]._Owner->ZPos[iBody[ib].markers[I].supp_k[s]]) / iBody[ib]._Owner->dx,
					iBody[ib].markers[J].dilation
					);
#else
				Delta_J =
					ibm_deltakernel(
					(iBody[ib].markers[J].position[0] - iBody[ib]._Owner->XPos[iBody[ib].markers[I].supp_i[s]]) / iBody[ib]._Owner->dx,
					iBody[ib].markers[J].dilation
					) *
					ibm_deltakernel(
					(iBody[ib].markers[J].position[1] - iBody[ib]._Owner->YPos[iBody[ib].markers[I].supp_j[s]]) / iBody[ib]._Owner->dx,
					iBody[ib].markers[J].dilation
					);
#endif
				// Multiply by local area (or volume in 3D)
				A[I][J] += Delta_I * Delta_J * iBody[ib].markers[I].local_area;
			}

			// Multiply by arc length between markers in lattice units
			A[I][J] = A[I][J] * (iBody[ib].spacing / iBody[ib]._Owner->dx);

		}

	}


#ifdef L_IBM_DEBUG
	// DEBUG -- write out A
	std::ofstream Aout;
	Aout.open(GridUtils::path_str + "/Amatrix_" + std::to_string(ib) + "_rank" + std::to_string(MpiManager::my_rank) + ".out");
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


	///////////////////
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
