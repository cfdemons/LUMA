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
#include "../inc/FEMNode.h"

// *****************************************************************************
///	\brief	Default constructor for finite element node class
FEMNode::FEMNode () {

	// Set members to default values
	ID = 0;
	angles0 = 0.0;
	angles = 0.0;
}


// *****************************************************************************
///	\brief	Default destructor for finite element node class
FEMNode::~FEMNode () {
}


// *****************************************************************************
///	\brief	Custom constructor for building FEM node from input parameters
///
///	\param	idx				marker ID
///	\param	x				x-position
///	\param	y				y-position
///	\param	z				z-position
///	\param	inputAngles		node angle
FEMNode::FEMNode (int idx, double x, double y, double z, std::vector<double> &inputAngles) {

	// Set values
	ID = idx;

	// Set positions
	position0.push_back(x);
	position0.push_back(y);
	position0.push_back(z);

	// Resize vectors
	position.resize(3);

	// Set values
	position = position0;
	angles0 = inputAngles[0] * L_PI / 180.0;
	angles = angles0;
}
