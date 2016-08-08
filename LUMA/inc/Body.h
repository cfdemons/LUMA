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
#include <vector>

class GridObj;

// Markers may be of different types depending on the body
template <typename MarkerType>

/** Represents a general body **/
class Body
{

	// Default Constructor / Destructor //
public:
	Body(void)
	{
	};
	~Body(void)
	{
	};

	Body(GridObj* g)
	{
		this->_Owner = g;
	};

	// Protected Members //
	
protected:
	double spacing;						// Spacing of the markers in physical units
	std::vector<MarkerType> markers;	// Array of markers which make up the body
	bool closed_surface;				// Flag to specify whether or not it is a closed surface (for output)
	GridObj* _Owner;					// Pointer to owning grid
	

	// Methods //

	// Add marker to the body
	void addMarker(double x, double y, double z)
	{
		markers.emplace_back(x, y, z);	// Add a new marker object to the array
	}

};
