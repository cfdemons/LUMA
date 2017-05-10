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

protected:
	// Per marker force storage
	double forceX = 0.0;	///< Instantaneous X-direction force on marker
	double forceY = 0.0;	///< Instantaneous Y-direction force on marker
	double forceZ = 0.0;	///< Instantaneous Z-direction force on marker

};

#endif
