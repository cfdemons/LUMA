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

/* This file defines the constructors and methods for the immersed boundary body object.
*/

#include "../inc/stdafx.h"
#include "../inc/FEMNode.h"

// ***************************************************************************************************
/// \brief Default constructor for FEM body.
FEMNode::FEMNode () {

	// Set members to default values
	ID = 0;
}


// ***************************************************************************************************
/// \brief Default destructor for FEM node.
FEMNode::~FEMNode () {
}


// ***************************************************************************************************
/// \brief Custom constructor for building FEM node from input parameters
FEMNode::FEMNode (int idx, double x, double y, double z, std::vector<double> &inputAngles) {

	// Set values
	ID = idx;

	// Resize vectors
	position.resize(3);
	angles0.resize(2);
	angles.resize(2);

	// Set positions
	position0.push_back(x);
	position0.push_back(y);
	position0.push_back(z);

	// Set values
	position = position0;
	angles0 = inputAngles;
	angles = angles0;
}
