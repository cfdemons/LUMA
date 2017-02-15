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

#ifndef BFLMARKER_H
#define BFLMARKER_H

#include "stdafx.h"
#include "Marker.h"

/// \brief	BFL marker.
///
///			This class declaration is for a BFL Lagrange point.
///			A collection of these points form BFL body.
class BFLMarker :
	public Marker
{

	// Make ObjectManager a friend class so it can access the protected data of BFLMarker objects
	friend class ObjectManager;

	// Allow BFLbody to access the marker positions etc.
	friend class BFLBody;

public:
	// Default constructor and destructor
	BFLMarker(void);
	~BFLMarker(void);

	// Custom constructor when positions are passed
	BFLMarker(double x, double y, double z, int markerID, GridObj const * const body_owner);

};

#endif
