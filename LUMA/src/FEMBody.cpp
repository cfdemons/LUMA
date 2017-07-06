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
#include "../inc/FEMBody.h"

// ***************************************************************************************************
/// \brief Default constructor for FEM body.
FEMBody::FEMBody () {

	// Set members to default values
	iBodyPtr = NULL;
	DOFsPerNode = 0;
	DOFsPerElement = 0;
	systemDOFs = 0;
	NRIterations = 0;
	NRResidual = 0.0;
	BC_DOFs = 0;
}

// ***************************************************************************************************
/// \brief Default destructor for FEM body.
FEMBody::~FEMBody () {
}

// ***************************************************************************************************
/// \brief Custom constructor for building FEM body from inputs.
FEMBody::FEMBody (IBBody *iBody, std::vector<double> &start_position, double length,
		double height, double depth, std::vector<double> &angles, int nElements, bool clamped, double density, double E) {

	// Set members to default values
	iBodyPtr = iBody;
	DOFsPerNode = 3;
	DOFsPerElement = 6;
	systemDOFs = (nElements + 1) * DOFsPerNode;
	NRIterations = 0;
	NRResidual = 0.0;

	// Set number of DOFs to remove in BC
	if (clamped == true)
		BC_DOFs = 3;
	else
		BC_DOFs = 2;

	// Get horizontal and vertical angles
	double body_angle_v = angles[0];
#if (L_DIM == 3)
	double body_angle_h = angles[1];
#else
	double body_angle_h = 0.0;
#endif

	// Compute spacing
	double spacing = length / nElements;							// Physical spacing between markers
	double spacing_h = spacing * cos(body_angle_v * L_PI / 180);	// Local spacing projected onto the horizontal plane

	// Add all FEM nodes
	double x, y, z;
	for (int i = 0; i < nElements + 1; i++) {
		x = start_position[eXDirection] + i * spacing_h * cos(body_angle_h * L_PI / 180.0);
		y = start_position[eYDirection] + i * spacing * sin(body_angle_v * L_PI / 180.0);
		z = start_position[eZDirection] + i * spacing_h * sin(body_angle_h * L_PI / 180.0);
		node.emplace_back(i, x, y, z, angles);
	}

// If 2D then set depth to 1 lattice unit
#if (L_DIMS == 2)
	depth = iBodyPtr->_Owner->dh;
#endif

	// Add all FEM elements
	for (int i = 0; i < nElements; i++) {
		element.emplace_back(i, DOFsPerElement, spacing, height, depth, angles, density, E);
	}

	// Get number of IBM nodes
//	int nIBMNodes = floor(length / iBodyPtr->_Owner->dh) + 1;

	// Resize the matrices and set to zero
	M.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	K_L.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	K_NL.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	R.resize(systemDOFs, 0.0);
	F.resize(systemDOFs, 0.0);
	U.resize(systemDOFs, 0.0);
	delU.resize(systemDOFs, 0.0);
	Udot.resize(systemDOFs, 0.0);
	Udotdot.resize(systemDOFs, 0.0);
}
