/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

#include "../inc/stdafx.h"
#include "../inc/FEMBody.h"
#include "../inc/IBMarker.h"
#include "../inc/IBBody.h"

// *****************************************************************************
/// \brief	Default constructor for FEM body
FEMBody::FEMBody () {

	// Set members to default values
	iBodyPtr = NULL;
	DOFsPerNode = 0;
	DOFsPerElement = 0;
	systemDOFs = 0;
	it = 0;
	res = 0.0;
	timeav_FEMIterations = 0.0;
	timeav_FEMResidual = 0.0;
	BC_DOFs = 0;
}

// *****************************************************************************
/// \brief	Default destructor for FEM body.
FEMBody::~FEMBody () {
}

// *****************************************************************************
///	\brief	Custom constructor for building FEM body from inputs
///
///	\param	iBody			pointer to owning IBBody
///	\param	start_position	start position of FEMBody
///	\param	length			length of FEMBody
///	\param	height			height of FEMBody
///	\param	depth			depth of FEMBody (set to dh in 2D cases)
///	\param	angles			angles of FEMBody
///	\param	nElements		number of elements in FEMBody
///	\param	clamped			boundary condition
///	\param	density			material density
///	\param	E				Young's Modulus
FEMBody::FEMBody (IBBody *iBody, std::vector<double> &start_position, double length,
		double height, double depth, std::vector<double> &angles, int nElements, bool clamped, double density, double E) {

	// Set members to default values
	iBodyPtr = iBody;
	DOFsPerNode = 3;
	DOFsPerElement = 6;
	systemDOFs = (nElements + 1) * DOFsPerNode;
	it = 0;
	res = 0.0;
	timeav_FEMIterations = 0.0;
	timeav_FEMResidual = 0.0;

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
		elements.emplace_back(this, i, DOFsPerElement, spacing, height, depth, angles, density, E);
	}

	// Get number of IBM nodes
	int nIBMNodes = static_cast<int>(std::floor(length / iBodyPtr->_Owner->dh)) + 1;
	int nFEMNodes = static_cast<int>(nodes.size());

	// Compute IBM-FEM conforming parameters
	computeNodeMapping(nIBMNodes, nFEMNodes);

	// Resize the matrices and set to zero
	M.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	K.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	R.resize(systemDOFs, 0.0);
	F.resize(systemDOFs, 0.0);
	U.resize(systemDOFs, 0.0);
	U_n.resize(systemDOFs, 0.0);
	delU.resize(systemDOFs, 0.0);
	Udot.resize(systemDOFs, 0.0);
	Udot_n.resize(systemDOFs, 0.0);
	Udotdot.resize(systemDOFs, 0.0);
	Udotdot_n.resize(systemDOFs, 0.0);
}


// *****************************************************************************
///	\brief	Main outer routine for FEM solver
void FEMBody::dynamicFEM () {

	// Set values to start of timestep
	U = U_n;
	Udot = Udot_n;
	Udotdot = Udotdot_n;

	// Reset the nodal and elemental values
	updateFEMValues();

	// Set R vector to zero
	fill(R.begin(), R.end(), 0.0);

	// Construct load vector as invariant during Newton-Raphson iterations
	for (size_t el = 0; el < elements.size(); el++)
		elements[el].loadVector();

	// While loop parameters
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

	// Calculate velocities and accelerations
	finishNewmark();

	// Update IBM markers
	updateIBMarkers();
}

// *****************************************************************************
///	\brief	Newton-Raphson routine for solving non-linear FEM
void FEMBody::newtonRaphsonIterator () {

	// Set matrices to zero
	fill(F.begin(), F.end(), 0.0);
	for (int i = 0; i < systemDOFs; i++) {
		fill(M[i].begin(), M[i].end(), 0.0);
		fill(K[i].begin(), K[i].end(), 0.0);
	}

	// Loop through and build global matrices
	for (size_t el = 0; el < elements.size(); el++) {

		// Build force vector
		elements[el].forceVector();

		// Build mass matrix
		elements[el].massMatrix();

		// Build stiffness matrix
		elements[el].stiffMatrix();
	}

	// Apply Newmark scheme (using Newmark coefficients)
	setNewmark();

	// Solve linear system using LAPACK library
	delU = GridUtils::solveLinearSystem(K, F, BC_DOFs);

	// Add deltaU to U
	for (int i = 0; i < systemDOFs; i++) {
		U[i] = U[i] + delU[i];
	}

	// Update FEM positions
	updateFEMValues();
}

// *****************************************************************************
///	\brief	Check convergence of the Newton-Raphson scheme
///
/// \return	residual
double FEMBody::checkNRConvergence () {

	// Get original length of body
	double L = elements[0].length0 * elements.size();

	// Get the norm but divide by length so it is resolution independent and return
	return sqrt(GridUtils::dotprod(delU, delU) / static_cast<double> (delU.size())) / L;
}


// *****************************************************************************
///	\brief	Newmark-Beta scheme for getting FEM velocities and accelerations
void FEMBody::finishNewmark () {

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


// *****************************************************************************
///	\brief	First step in Newmark-Beta time integration
///
///	\param	M_hat		mass matrix with BCs applied
///	\param	K_hat		stiffness matrix with BCs applied
///	\param	RmF_hat		balanced load vector with BCs applied
void FEMBody::setNewmark () {

	// Newmark-beta method for time integration
	double Dt = iBodyPtr->_Owner->dt;
	double a0, a2, a3;
	a0 = 1.0 / (L_NB_ALPHA * SQ(Dt));
	a2 = 1.0 / (L_NB_ALPHA * Dt);
	a3 = 1.0 / (2.0 * L_NB_ALPHA) - 1.0;

	// Calculate effective load vector
	std::vector<double> Meff_hat(systemDOFs, 0.0);
	for (int i = 0; i < systemDOFs; i++) {
		Meff_hat[i] = a0 * (U_n[i] - U[i]) + a2 * Udot[i] + a3 * Udotdot[i];
	}

	// Multiply with mass matrix to get inertia forces
	std::vector<double> MF_hat = GridUtils::matrix_multiply(M, Meff_hat);

	// Calculate effective load vector and stiffness matrix
	for (int i = 0; i < systemDOFs; i++) {

		// Effective load
		F[i] = R[i] - F[i] + MF_hat[i];

		// Effective stiffness
		for (int j = 0; j < systemDOFs; j++) {
			K[i][j] += a0 * M[i][j];
		}
	}
}

// *****************************************************************************
///	\brief	Update the new FEM node data using the displacements
void FEMBody::updateFEMValues () {

	// Set the new positions in the particle_struct
	for (size_t n = 0; n < nodes.size(); n++) {

		// Positions
		for (int d = 0; d < L_DIMS; d++)
			nodes[n].position[d] = nodes[n].position0[d] + U[n*DOFsPerNode+d];

		// Set angle
		nodes[n].angles = nodes[n].angles0 + U[n*DOFsPerNode+L_DIMS];
	}

	// Set the new angles and lengths of the elements
	std::vector<double> elVector;
	for (size_t el = 0; el < elements.size(); el++) {

		// Set the new angles and lengths of the elements
		elVector = GridUtils::subtract(nodes[el+1].position, nodes[el].position);
		elements[el].angles = atan2(elVector[eYDirection], elVector[eXDirection]);
		elements[el].length = GridUtils::vecnorm(elVector);

		// Set new transformation matrix
		elements[el].T[0][0] = elements[el].T[1][1] =  elements[el].T[3][3] = elements[el].T[4][4] = cos(elements[el].angles);
		elements[el].T[0][1] = elements[el].T[3][4] = sin(elements[el].angles);
		elements[el].T[1][0] = elements[el].T[4][3] = -sin(elements[el].angles);
		elements[el].T[2][2] = elements[el].T[5][5] =  1.0;
	}
}

// *****************************************************************************
///	\brief	Update the new IBMaker data using the FEM data
void FEMBody::updateIBMarkers() {

	// Parameters
	std::vector<double> dashU;
	std::vector<double> dashUdot;
	std::vector<std::vector<double>> T(L_DIMS, std::vector<double>(L_DIMS, 0.0));

	// Loop through all IBM nodes
	for (size_t node = 0; node < IBNodeParents.size(); node++) {

		// Get pointer to element and zeta value
		FEMElement *el = &(elements[IBNodeParents[node].elementID]);
		double zeta = IBNodeParents[node].zeta;

		// Get the element values
		dashU = el->disassembleGlobalMat(U);
		dashUdot = el->disassembleGlobalMat(Udot);

		// Get element values in local coordinates
		dashU = GridUtils::matrix_multiply(el->T, dashU);
		dashUdot = GridUtils::matrix_multiply(el->T, dashUdot);

		// Multiply by shape functions
		dashU = el->shapeFuns(dashU, zeta);
		dashUdot = el->shapeFuns(dashUdot, zeta);

		// Get subset of transformation matrix
		T = {{el->T[0][0], el->T[0][1]},
			 {el->T[1][0], el->T[1][1]}};

		// Get the IB node displacements in global coordinates
		dashU = GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), dashU);
		dashUdot = GridUtils::vecmultiply(iBodyPtr->_Owner->dt / iBodyPtr->_Owner->dh, GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), dashUdot));

		// Set the IBM node
		for (int d = 0; d < L_DIMS; d++) {
			iBodyPtr->markers[node].position[d] = iBodyPtr->markers[node].position0[d] + dashU[d];
			iBodyPtr->markers[node].markerVel_km1[d] = iBodyPtr->markers[node].markerVel[d];
			iBodyPtr->markers[node].markerVel[d] = L_RELAX * dashUdot[d] + (1.0 - L_RELAX) * iBodyPtr->markers[node].markerVel_km1[d];
		}
	}
}


// *****************************************************************************
///	\brief	Compute the mapping parameters for FEM-IBM nodes
///
///	\param nIBMNodes	number of IBM nodes in body
///	\param nFEMNodes	number of FEM nodes in body
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
			IBNodeParents[i].elementID = static_cast<int>(std::floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)))) - 1;
			IBNodeParents[i].zeta = 1.0;
		}
		else {
			IBNodeParents[i].elementID = static_cast<int>(std::floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0))));
			IBNodeParents[i].zeta = fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) * 2.0 - 1.0;
		}
	}

	// Loop through elements
	double node1, node2;
	for (size_t el = 0; el < elements.size(); el++) {

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

// *****************************************************************************
///	\brief	Default constructor for parent element class
FEMBody::IBMParentElements::IBMParentElements() {

	// Set default values
	elementID = 0;
	zeta = 0.0;
}

// *****************************************************************************
///	\brief	Default destructor for parent element class
FEMBody::IBMParentElements::~IBMParentElements() {
}

