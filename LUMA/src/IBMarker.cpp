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
#include "../inc/IBMarker.h"

// ***************************************************************************************************
/// \brief Custom constructor for IBMarker (calls the non-default base constructor).
/// \param xPos			x-position of marker.
/// \param yPos			y-position of marker.
/// \param zPos			z-position of marker.
/// \param markerID		ID of marker within body
///	\param body_owner	Grid on which primary support point is to be found
IBMarker::IBMarker(double xPos, double yPos, double zPos, int markerID, GridObj const * const body_owner) : Marker(xPos, yPos, zPos, markerID, body_owner) {

	// Initialise all values to zero
	this->epsilon = 0.0;
	this->local_area = 1.0;
	this->dilation = 1.0;
	this->interpRho = 0.0;

	// Stationary point
	this->markerVel.push_back(0.0);
	this->markerVel.push_back(0.0);
	this->markerVel.push_back(0.0);

	// Resize vectors
	this->interpMom.resize(L_DIMS);
	this->force_xyz.resize(L_DIMS);

	// Set old position to initial position
	this->position0.resize(this->position.size());
	this->position0[eXDirection] = this->position[eXDirection];
	this->position0[eYDirection] = this->position[eYDirection];
	this->position0[eZDirection] = this->position[eZDirection];

	// Get rank which owns area where this marker is
	this->owningRank = GridUtils::getRankfromPosition(position);
}


// ***************************************************************************************************
