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
#ifndef PCPTS_H
#define PCPTS_H

#include "stdafx.h"
/// \brief	Class to hold point cloud data.
///
///			A container class for hold the X, Y and Z positions of 
///			points in a point cloud.
class PCpts {

public:

	/// Default constructor
	PCpts(void) {};

	/// Default destructor
	~PCpts(void) {};

	std::vector<double> x;		///< Vector of X positions
	std::vector<double> y;		///< Vector of Y positions
	std::vector<double> z;		///< Vector of Z positions
};

#endif