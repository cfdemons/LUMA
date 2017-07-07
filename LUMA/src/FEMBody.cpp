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
	U_n.resize(systemDOFs, 0.0);
	delU.resize(systemDOFs, 0.0);
	Udot.resize(systemDOFs, 0.0);
	Udotdot.resize(systemDOFs, 0.0);
}

// \brief Main outer routine for solving FEM.
//
void FEMBody::dynamicFEM () {

	// While loop parameters
	double res;
	int it;
	double TOL = 1e-10;
	double MAXIT = 20;

	// Set while counter to zero
	it = 0;

	// While loop for FEM solver
	do {

		// Solve and iterate over the system
		newtonRaphsonIterator();

		// Check residual
		res = checkNRConvergence();

		// Increment counter
		it++;

	} while (res > TOL && it < MAXIT);

	// Store the number of iterations
	NRIterations = it;
	NRResidual = res;

	// Calculate velocities and accelerations
	updateVelocityAndAcceleration();
}


// \brief Newton-Raphson routine for solve non-linear FEM
//
void FEMBody::newtonRaphsonIterator () {

	// Construct mass matrix
	constructMassMat();

	// Construct stiffness matrix
	constructStiffMat();

	// Construct the internal force vector
	constructFVector();

	// Construct the full non-linear stiffness matrix
	constructNLStiffMat();

	// Declare reduced matrices with BCs applied
	std::vector<std::vector<double>> M_hat(systemDOFs - BC_DOFs, std::vector<double>(systemDOFs - BC_DOFs, 0.0));
	std::vector<std::vector<double>> K_hat(systemDOFs - BC_DOFs, std::vector<double>(systemDOFs - BC_DOFs, 0.0));
	std::vector<double> RmF_hat(systemDOFs - BC_DOFs, 0.0);
	std::vector<double> delU_hat(systemDOFs - BC_DOFs, 0.0);

	// Apply the boundary conditions
	bcFEM(M_hat, K_hat, RmF_hat);

	// Apply Newmark scheme (using Newmark coefficients)
	setNewmark(M_hat, K_hat, RmF_hat);

	// Solve linear system using LAPACK library
	GridUtils::solveLinearSystem(K_hat, RmF_hat, delU_hat);

	// Assign displacement to delU
	for (int i = 0; i < systemDOFs - BC_DOFs; i++) {
		delU[i+BC_DOFs] = delU_hat[i];
	}

	// Add deltaU to U
	for (int i = 0; i < systemDOFs; i++) {
		U[i] = U[i] + delU[i];
	}

	// Update FEM positions
	updateFEMNodes();
}


// \brief Check convergence of the Newton-Raphson scheme
//
double FEMBody::checkNRConvergence () {

	// Get original length of body
	double L = elements[0].length0 * elements.size();

	// Normalise the displacement
	std::vector<double> delUStar(delU.size(), 0.0);
	for (int i = 0; i < delUStar.size(); i++)
		delUStar[i] = delU[i] / L;

	// Get the norm but divide by length so it is resolution independent and return
	return sqrt(GridUtils::dotprod(delUStar, delUStar) / delUStar.size());
}


// \brief Newmark-Beta scheme for getting FEM velocities and accelerations
//
void FEMBody::updateVelocityAndAcceleration () {

	// Get velocity and acceleration using Newmark-Beta scheme
	double Dt = iBodyPtr->_Owner->dt;
	double UdotdotStar;

	// Newmark coefficients
	double a6 = 1.0 / (L_NB_ALPHA * SQ(Dt));
	double a7 = -1.0 / (L_NB_ALPHA * Dt);
	double a8 = -(1.0 / (2.0 * L_NB_ALPHA) - 1.0);
	double a9 = Dt * (1.0 - L_NB_DELTA);
	double a10 = L_NB_DELTA * Dt;

	// Update velocities and accelerations
	for (int i = 0; i < systemDOFs; i++) {
		UdotdotStar = a6 * (U[i] - U_n[i]) + a7 * Udot[i] + a8 * Udotdot[i];
		Udot[i] = Udot[i] + a9 * Udotdot[i] + a10 * UdotdotStar;
		Udotdot[i] = UdotdotStar;
	}
}


// \brief Construct mass matrix
//
void FEMBody::constructMassMat () {

}


// \brief Construct linear stiffness matrix
//
void FEMBody::constructStiffMat () {

}


// \brief Construct internal force vector
//
void FEMBody::constructFVector () {

}


// \brief Construct non-linear stiffness matrix
//
void FEMBody::constructNLStiffMat () {

}


// \brief Apply BCs by removing elements in global matrices
//
///	\param	M_hat	mass matrix with BCs about to be applied
///	\param	K_hat	stiffness matrix with BCs about to be applied
///	\param	RmF_hat	balanced load vector with BCs about to be applied
void FEMBody::bcFEM (std::vector<std::vector<double>> &M_hat, std::vector<std::vector<double>> &K_hat, std::vector<double> &RmF_hat) {

}


// \brief First step in Newmar-Beta time integration
//
///	\param	M_hat	mass matrix with BCs applied
///	\param	K_hat	stiffness matrix with BCs applied
///	\param	RmF_hat	balanced load vector with BCs applied
void FEMBody::setNewmark (std::vector<std::vector<double>> &M_hat, std::vector<std::vector<double>> &K_hat, std::vector<double> &RmF_hat) {

}


// \brief Update the FEM node data using the new displacements
//
void FEMBody::updateFEMNodes () {

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

