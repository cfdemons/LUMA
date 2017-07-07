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
		nodes.emplace_back(i, x, y, z, angles);
	}

// If 2D then set depth to 1 lattice unit
#if (L_DIMS == 2)
	depth = iBodyPtr->_Owner->dh;
#endif

	// Add all FEM elements
	for (int i = 0; i < nElements; i++) {
		elements.emplace_back(i, DOFsPerElement, spacing, height, depth, angles, density, E);
	}

	// Get number of IBM nodes
	int nIBMNodes = floor(length / iBodyPtr->_Owner->dh) + 1;
	int nFEMNodes = nodes.size();

	// Compute IBM-FEM conforming parameters
	computeNodeMapping(nIBMNodes, nFEMNodes);

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

// \brief Main outer routine for solving FEM.
//
void FEMBody::dynamicFEM () {

	std::cout << "FEM TEST" << std::endl;
}


// \brief Compute the mapping parameters for FEM-IBM nodes.
//
///	\param	nIBMNodes	number of IBM nodes in body.
///	\param	nFEMNodes	number of FEM nodes in body.
void FEMBody::computeNodeMapping (int nIBMNodes, int nFEMNodes) {

	// Set the IBM parent elements
	IBNodeParents.resize(nIBMNodes);

	// Set first node first
	IBNodeParents[0].elementID = 0;
	IBNodeParents[0].zeta = -1.0;

	// Now loop through and set values
	for (int i = 1; i < nIBMNodes; i++) {

		// Check if remainder is zero
		if (fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) == 0.0) {
			IBNodeParents[i].elementID = floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0))) - 1;
			IBNodeParents[i].zeta = 1.0;
		}
		else {
			IBNodeParents[i].elementID = floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)));
			IBNodeParents[i].zeta = fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) * 2.0 - 1.0;
		}
	}

	// Loop through elements
	double node1, node2;
	for (int el = 0; el < elements.size(); el++) {

		// Loop through all nodes and scale range to local coordinates for element
		for (int node = 0; node < nIBMNodes; node++) {
			node1 = -1 + 2.0 * ((node - 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - el) / 1.0;
			node2 = -1 + 2.0 * ((node + 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - el) / 1.0;

			// Check if any points lie within element coordinate range
			if ((node1 > -1.0 && node1 < 1.0) || (node2 > -1.0 && node2 < 1.0)) {

				// Sort out end nodes where one point lie outside element
				if (node1 < -1.0)
					node1 = -1.0;
				if (node2 > 1.0)
					node2 = 1.0;

				// Call constructor for chile IB point
				elements[el].IBChildNodes.emplace_back(node, node1, node2);
			}
		}
	}
}


// \brief Constructor for parent element class.
FEMBody::IBMParentElements::IBMParentElements() {

	// Set default values
	elementID = 0;
	zeta = 0.0;
}

// \brief Destructor for parent element class.
FEMBody::IBMParentElements::~IBMParentElements() {
}

