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
#include "../inc/FEMElement.h"
#include "../inc/FEMBody.h"
#include "../inc/IBBody.h"
#include "../inc/IBMarker.h"

// *****************************************************************************
///	\brief	Default constructor for finite element class
FEMElement::FEMElement () {

	// Set members to default values
	fPtr = NULL;
	ID = 0;
	length0 = 0;
	length = 0;
	area = 0;
	E = 0;
	density = 0;
	angles = 0.0;
	I = 0.0;
}


// *****************************************************************************
///	\brief	Default destructor for finite element class
FEMElement::~FEMElement () {
}


// *****************************************************************************
///	\brief	Custom constructor to build FEM element from inputs
///
///	\param	i				element ID
///	\param	DOFs			number of degrees of freedom per element
///	\param	spacing			element length
///	\param	height			element height
///	\param	depth			element depth (set to dh for 2D cases)
///	\param	inputAngles		angle of element
///	\param	inputDensity	material density
///	\param	inputE			Young's modulus
FEMElement::FEMElement (FEMBody *f, int i, int elDOFs, double spacing, double height,
		double depth, std::vector<double> &inputAngles, double inputDensity, double inputE) {

	// Set pointer
	fPtr = f;

	// Set values and geometry
	ID = i;
	length0 = spacing;
	length = length0;
	area = height * depth;

	// Set angles
	angles = inputAngles[0] * L_PI / 180.0;

	// Get the second moment areas
	I = depth * TH(height) / 12.0;

	// Material properties
	E = inputE;
	density = inputDensity;

	// Get DOFs
	for (int i = 0; i < fPtr->DOFsPerElement; i++)
		DOFs.push_back(ID * fPtr->DOFsPerNode + i);

	// Resize transformation matrix
	T.resize(elDOFs, std::vector<double>(elDOFs, 0.0));

	// Set to correct values
	T[0][0] = T[1][1] =  T[3][3] = T[4][4] = cos(angles);
	T[0][1] = T[3][4] = sin(angles);
	T[1][0] = T[4][3] = -sin(angles);
	T[2][2] = T[5][5] =  1.0;

	// Initialise internal forces to zero
	F.resize(elDOFs, 0.0);
}


// *****************************************************************************
///	\brief	Get shape functions given natural value
///
///	\param vec		vector containing the values in natural coordinates
///	\param zeta		the natural coordinate
/// \return value at the coordinate
std::vector<double> FEMElement::shapeFuns(const std::vector<double> &vec, double zeta) {

	// Results vector
	std::vector<double> resVec(L_DIMS);

	// Use shape functions to calculate values
	double N0 = 1.0 - (zeta + 1.0) / 2.0;
	double N1 = 1.0 - 3.0 * SQ((zeta + 1.0) / 2.0) + 2.0 * TH((zeta + 1.0) / 2.0);
	double N2 = ((zeta + 1.0) / 2.0 - 2.0 * SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * length;
	double N3 = (zeta + 1.0) / 2.0;
	double N4 = 3.0 * SQ((zeta + 1.0) / 2.0) - 2.0 * TH((zeta + 1.0) / 2.0);
	double N5 = (-SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * length;

	// Calculate values using shape functions
	resVec[eXDirection] = vec[0] * N0 + vec[3] * N3;
	resVec[eYDirection] = vec[1] * N1 + vec[2] * N2 + vec[4] * N4 + vec[5] * N5;

	// Return
	return resVec;
}


// *****************************************************************************
///	\brief	Construct elemental load vector
void FEMElement::loadVector () {

	// Initialise arrays and matrices for calculating load vector
	std::vector<std::vector<double>> Tsub(L_DIMS, std::vector<double>(L_DIMS, 0.0));
	std::vector<double> Rlocal(fPtr->DOFsPerElement, 0.0);
	std::vector<double> RGlobal;
	std::vector<double> F;

	// Get force scaling parameter
	double forceScale = fPtr->iBodyPtr->_Owner->dm / SQ(fPtr->iBodyPtr->_Owner->dt);

	// Now loop through all child IB nodes this element has
	for (size_t node = 0; node < IBChildNodes.size(); node++) {

		// Get IB node and integration ranges
		int IBnode = IBChildNodes[node].nodeID;
		double a = IBChildNodes[node].zeta1;
		double b = IBChildNodes[node].zeta2;

		// Get subset of transpose matrix
		Tsub = {{T[0][0], T[0][1]},
			    {T[1][0], T[1][1]}};

		// Convert force to local coordinates
		F = GridUtils::matrix_multiply(Tsub, GridUtils::vecmultiply(fPtr->iBodyPtr->markers[IBnode].epsilon * 1.0 * forceScale, fPtr->iBodyPtr->markers[IBnode].force_xyz));

		// Get the nodal values by integrating over range of IB point
		Rlocal[0] = F[0] * 0.5 * length * (0.5 * b - 0.5 * a + 0.25 * SQ(a) - 0.25 * SQ(b));
		Rlocal[1] = F[1] * 0.5 * length * (0.5 * b - 0.5 * a - SQ(a) * SQ(a) / 16.0 + SQ(b) * SQ(b) / 16.0 + 3.0 * SQ(a) / 8.0 - 3.0 * SQ(b) / 8.0);
		Rlocal[2] = F[1] * 0.5 * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 - length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 + length * (b - a) / 8.0);
		Rlocal[3] = F[0] * 0.5 * length * (-0.25 * SQ(a) + 0.25 * SQ(b) + 0.5 * b - 0.5 * a);
		Rlocal[4] = F[1] * 0.5 * length * (0.5 * b - 0.5 * a + SQ(a) * SQ(a) / 16.0 - SQ(b) * SQ(b) / 16.0 - 3.0 * SQ(a) / 8.0 + 3.0 * SQ(b) / 8.0);
		Rlocal[5] = F[1] * 0.5 * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 + length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 - length * (b - a) / 8.0);

		// Get element internal forces
		RGlobal = GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), Rlocal);

		// Assemble into global vector
		assembleGlobalMat(RGlobal, fPtr->R);
	}
}


// *****************************************************************************
///	\brief	Construct elemental mass matrix
void FEMElement::massMatrix () {

	// Initialise arrays and matrices for calculating mass matrix
	std::vector<std::vector<double>> Mlocal(fPtr->DOFsPerElement, std::vector<double>(fPtr->DOFsPerElement, 0.0));

	// Coefficients
	double C1 = density * area * length0 / 420.0;

	// Set local mass matrix
	Mlocal[0][0] = C1 * 140.0;
	Mlocal[0][3] = C1 * 70.0;
	Mlocal[1][1] = C1 * 156.0;
	Mlocal[1][2] = C1 * 22.0 * length0;
	Mlocal[1][4] = C1 * 54;
	Mlocal[1][5] = C1 * (-13.0 * length0);
	Mlocal[2][2] = C1 * 4.0 * SQ(length0);
	Mlocal[2][4] = C1 * 13.0 * length0;
	Mlocal[2][5] = C1 * (-3.0 * SQ(length0));
	Mlocal[3][3] = C1 * 140.0;
	Mlocal[4][4] = C1 * 156.0;
	Mlocal[4][5] = C1 * (-22.0 * length0);
	Mlocal[5][5] = C1 * 4.0 * SQ(length0);

	// Copy to the lower half (symmetrical matrix)
	for (int i = 1; i < fPtr->DOFsPerElement; i++) {
		for (int j = 0; j < i; j++) {
			Mlocal[i][j] = Mlocal[j][i];
		}
	}

	// Multiply by transformation matrices to get global matrix for single element
	std::vector<std::vector<double>> Mglobal = GridUtils::matrix_multiply(GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), Mlocal), T);

	// Assemble into global matrix
	assembleGlobalMat(Mglobal, fPtr->M);
}


// *****************************************************************************
///	\brief	Construct elemental stiffness matrix
void FEMElement::stiffMatrix () {

	// Initialise arrays and matrices for calculating linear stiffness matrix
	std::vector<std::vector<double>> Klocal(fPtr->DOFsPerElement, std::vector<double>(fPtr->DOFsPerElement, 0.0));

	// Construct upper half of local stiffness matrix for single element
	Klocal[0][0] = E * area / length0;
	Klocal[0][3] = -E * area / length0;
	Klocal[1][1] = 12.0 * E * I / TH(length0);
	Klocal[1][2] = 6.0 * E * I / SQ(length0);
	Klocal[1][4] = -12.0 * E * I / TH(length0);
	Klocal[1][5] = 6.0 * E * I / SQ(length0);
	Klocal[2][2] = 4.0 * E * I / length0;
	Klocal[2][4] = -6.0 * E * I / SQ(length0);
	Klocal[2][5] = 2.0 * E * I / length0;
	Klocal[3][3] = E * area / length0;
	Klocal[4][4] = 12.0 * E * I / TH(length0);
	Klocal[4][5] = -6.0 * E * I / SQ(length0);
	Klocal[5][5] = 4.0 * E * I / length0;

	// Copy to the lower half (symmetrical matrix)
	for (int i = 1; i < fPtr->DOFsPerElement; i++) {
		for (int j = 0; j < i; j++) {
			Klocal[i][j] = Klocal[j][i];
		}
	}

	// Now add nonlinear part
	double F0 = -F[0];
	double V0 = F[4];

	// Nonlinear contributions to stiffness matrix
	Klocal[0][1] += -V0 / length0;
	Klocal[0][4] += V0 / length0;
	Klocal[1][0] += -V0 / length0;
	Klocal[1][1] += F0 / length0;
	Klocal[1][3] += V0 / length0;
	Klocal[1][4] += -F0 / length0;
	Klocal[3][1] += V0 / length0;
	Klocal[3][4] += -V0 / length0;
	Klocal[4][0] += V0 / length0;
	Klocal[4][1] += -F0 / length0;
	Klocal[4][3] += -V0 / length0;
	Klocal[4][4] += F0 / length0;

	// Multiply by transformation matrices to get global matrix for single element
	std::vector<std::vector<double>> Kglobal = GridUtils::matrix_multiply(GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), Klocal), T);

	// Assemble into global matrix
	assembleGlobalMat(Kglobal, fPtr->K);
}


// *****************************************************************************
///	\brief	Construct elemental internal force vector
void FEMElement::forceVector () {

	// Get node angles
	double angle1 = fPtr->nodes[ID].angles;
	double angle2 = fPtr->nodes[ID+1].angles;

	// Calculate the local displacements
	double u = (SQ(length) - SQ(length0)) / (length + length0);
	double theta1 = atan2(cos(angles) * sin(angle1) - sin(angles) * cos(angle1), cos(angles) * cos(angle1) + sin(angles) * sin(angle1));
	double theta2 = atan2(cos(angles) * sin(angle2) - sin(angles) * cos(angle2), cos(angles) * cos(angle2) + sin(angles) * sin(angle2));

	// Calculate internal forces for beam
	double F0 = (E * area / length0) * u;
	double M1 = (2 * E * I / length0) * (2.0 * theta1 + theta2);
	double M2 = (2 * E * I / length0) * (theta1 + 2.0 * theta2);

	// Set the internal local nodal forces
	F[0] = -F0;
	F[1] = (1.0 / length0) * (M1 + M2);
	F[2] = M1;
	F[3] = F0;
	F[4] = -(1.0 / length0) * (M1 + M2);
	F[5] = M2;

	// Get element internal forces
	std::vector<double> FGlobal = GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), F);

	// Assemble into global vector
	assembleGlobalMat(FGlobal, fPtr->F);
}


// *****************************************************************************
///	\brief	Assemble global vector from local elemental vector
///
///	\param	localVec			elemental vector
///	\param	globalVec			global vector
void FEMElement::assembleGlobalMat (const std::vector<double> &localVec, std::vector<double> &globalVec) {

	// Loop through and set
	for (size_t i = 0; i < localVec.size(); i++)
		globalVec[DOFs[i]] += localVec[i];
}


// *****************************************************************************
///	\brief	Assemble global matrix from local elemental matrix
///
///	\param	localMat			elemental matrix
///	\param	globalMat			global matrix
void FEMElement::assembleGlobalMat (const std::vector<std::vector<double>> &localMat, std::vector<std::vector<double>> &globalMat) {

	// Get rows and cols
	size_t rows = localMat.size();
	size_t cols = localMat[0].size();

	// Now loop through and set
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			globalMat[DOFs[i]][DOFs[j]] += localMat[i][j];
		}
	}
}


// *****************************************************************************
///	\brief	Disassemble global vector into local elemental vector
///
///	\param	globalVec			global vector
/// \return	local vector
std::vector<double> FEMElement::disassembleGlobalMat (const std::vector<double> &globalVec) {

	// Declare vector
	std::vector<double> localVec(fPtr->DOFsPerElement);

	// Loop through and set
	for (size_t i = 0; i < localVec.size(); i++)
		localVec[i] = globalVec[DOFs[i]];

	// Return
	return localVec;
}

// *****************************************************************************
///	\brief	Default constructor for child IB marker class
FEMElement::FEMChildNodes::FEMChildNodes() {

	// Set default values
	nodeID = 0;
	zeta1 = 0.0;
	zeta2 = 0.0;
}


// *****************************************************************************
///	\brief	Default destructor for child IB marker class
FEMElement::FEMChildNodes::~FEMChildNodes() {
}


// *****************************************************************************
///	\brief	Custom constructor for child IB marker class
///
///	\param	node			node ID
///	\param	zetaA			start of integration range
///	\param	zetaB			end of integration range
FEMElement::FEMChildNodes::FEMChildNodes(int node, double zetaA, double zetaB) {

	// Set default values
	nodeID = node;
	zeta1 = zetaA;
	zeta2 = zetaB;
}
