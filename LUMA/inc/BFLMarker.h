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

#pragma once
#include "Marker.h"

// Definition for the BFLMarker class which represents a marker which defines a BFL object //
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

protected:


	/*
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/





	/*
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/




};

