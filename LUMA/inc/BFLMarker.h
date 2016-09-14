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

	// Allow BFLbody to access the marker positions etc.
	friend class BFLBody;

public:
	// Default constructor and destructor
	BFLMarker(void);
	~BFLMarker(void);

	// Constructor
	BFLMarker(double x, double y, double z);

};

#endif