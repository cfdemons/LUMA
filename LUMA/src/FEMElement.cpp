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
	angles = 0.0;
	I = 0.0;
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
	angles = inputAngles[0];

	// Get the second moment areas
	I = depth * height * height * height / 12.0;

	// Material properties
	E = inputE;
	density = inputDensity;

	// Initialise internal forces to zero
	F.resize(DOFs, 0.0);
}


// \brief Constructor for parent element class.
FEMElement::FEMChildNodes::FEMChildNodes() {

	// Set default values
	nodeID = 0;
	zeta1 = 0.0;
	zeta2 = 0.0;
}

// \brief Destructor for parent element class.
FEMElement::FEMChildNodes::~FEMChildNodes() {
}

// \brief Constructor for parent element class.
FEMElement::FEMChildNodes::FEMChildNodes(int node, double zetaA, double zetaB) {

	// Set default values
	nodeID = node;
	zeta1 = zetaA;
	zeta2 = zetaB;
}
