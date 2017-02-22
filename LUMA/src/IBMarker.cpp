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
#include "../inc/IBMarker.h"

// ***************************************************************************************************
/// \brief Custom constructor with position.
/// \param xPos			x-position of marker.
/// \param yPos			y-position of marker.
/// \param zPos			z-position of marker.
///	\param body_owner	Grid on which primary support point is to be found
/// \param isFlexible	flag to indicate whether marker is movable or not.
IBMarker::IBMarker(double xPos, double yPos, double zPos, GridObj const * const body_owner, bool isFlexible) {

	// Assign position and type
	this->isFlexible = isFlexible;
	this->position.push_back(xPos);
	this->position.push_back(yPos);
	this->position.push_back(zPos);

	// Stationary point
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);

	// Resize vectors
	this->fluid_vel.resize(L_DIMS);
	this->force_xyz.resize(L_DIMS);

	std::vector<int> ijk;
	GridUtils::getEnclosingVoxel(xPos, yPos, zPos, body_owner, &ijk);
	supp_i.push_back(ijk[eXDirection]);
	supp_j.push_back(ijk[eYDirection]);
	supp_k.push_back(ijk[eZDirection]);

#ifdef L_BUILD_FOR_MPI
	support_rank.push_back(MpiManager::getInstance()->my_rank);
#else
	support_rank.push_back(0);
#endif

}


// ***************************************************************************************************
