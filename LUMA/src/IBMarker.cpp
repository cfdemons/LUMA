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
/// \brief Custom constructor for IBMarker (calls the non-default base constructor).
/// \param xPos			x-position of marker.
/// \param yPos			y-position of marker.
/// \param zPos			z-position of marker.
///	\param body_owner	Grid on which primary support point is to be found
IBMarker::IBMarker(double xPos, double yPos, double zPos, GridObj const * const body_owner) : Marker(xPos, yPos, zPos, body_owner) {

	// Initialise all values to zero
	this->epsilon = 0.0;
	this->local_area = 0.0;
	this->dilation = 0.0;

	// Stationary point
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);

	// Resize vectors
	this->fluid_vel.resize(L_DIMS);
	this->force_xyz.resize(L_DIMS);

	// Set old position to initial position
	this->position_old.resize(this->position.size());
	this->position_old[eXDirection] = this->position[eXDirection];
	this->position_old[eYDirection] = this->position[eYDirection];
	this->position_old[eZDirection] = this->position[eZDirection];
}


// ***************************************************************************************************
