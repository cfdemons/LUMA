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
#include "../inc/BFLMarker.h"

// Implementation for the BFLMarker class //


/// Default constructor
BFLMarker::BFLMarker(void)
{
}

/// Default destructor
BFLMarker::~BFLMarker(void)
{
}

/// \brief Custom constructor with position.
/// \param x			x-position of marker
/// \param y			y-position of marker
/// \param z			z-position of marker
/// \param markerID		ID of marker within body
///	\param body_owner	Grid on which primary support is to be found.
BFLMarker::BFLMarker(double x, double y, double z, int markerID, GridObj const * const body_owner) : Marker(x, y, z, markerID, body_owner)
{
}
