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
#include "../inc/FEMElement.h"

// ***************************************************************************************************
/// \brief Default constructor for FEM element.
FEMElement::FEMElement () {

	// Set members to default values
	ID = 0;
	length0 = 0;
	length = 0;
	area = 0;
	E = 0;
	density = 0;
}

// ***************************************************************************************************
/// \brief Default constructor for FEM element.
FEMElement::~FEMElement () {
}


// ***************************************************************************************************
/// \brief Custom constructor to build FEM element from inputs.
FEMElement::FEMElement (int i, int DOFs, double spacing, double height,
		double depth, std::vector<double> &inputAngles, double inputDensity, double inputE) {


	// Set values and geometry
	ID = i;
	length0 = spacing;
	length = length0;
	area = height * depth;

	// Set angles
	angles.resize(2);
	angles = inputAngles;

	// Get the second moment areas
	I.resize(L_DIMS-1);
	I[0] = depth * height * height * height / 12.0;
#if (L_DIMS == 3)
	I[1] = height * depth * depth * depth / 12.0;
#endif

	// Material properties
	E = inputE;
	density = inputDensity;

	// Initialise internal forces to zero
	F.resize(DOFs, 0.0);
}
