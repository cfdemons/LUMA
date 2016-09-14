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
/// \param xPos x-position of marker.
/// \param yPos y-position of marker.
/// \param zPos z-position of marker.
/// \param flex_rigid flag to indicate whether marker is movable or not.
IBMarker::IBMarker(double xPos, double yPos, double zPos, bool flex_rigid) {

	// Assign position and type
	this->flex_rigid = flex_rigid;
	this->position.push_back(xPos);
	this->position.push_back(yPos);
	this->position.push_back(zPos);

	// Stationary point
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);

	// Resize vectors
	this->fluid_vel.resize(L_dims);
	this->force_xyz.resize(L_dims);

}


// ***************************************************************************************************
