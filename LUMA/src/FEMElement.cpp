/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2018 The University of Manchester
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


// *****************************************************************************
///	\brief	Default constructor for finite element class
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
	angles_n = 0.0;
	length_n = 0.0;
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
FEMElement::FEMElement (int i, int DOFs, double spacing, double height,
		double depth, std::vector<double> &inputAngles, double inputDensity, double inputE) {


	// Set values and geometry
	ID = i;
	length0 = spacing;
	length = length0;
	length_n = length;
	area = height * depth;

	// Set angles
	angles = inputAngles[0] * L_PI / 180.0;
	angles_n = angles;

	// Get the second moment areas
	I = depth * TH(height) / 12.0;

	// Material properties
	E = inputE;
	density = inputDensity;

	// Resize transformation matrix
	T.resize(DOFs, std::vector<double>(DOFs, 0.0));
	T_n.resize(DOFs, std::vector<double>(DOFs, 0.0));

	// Set to correct values
	T[0][0] = T[1][1] =  T[3][3] = T[4][4] = cos(angles);
	T[0][1] = T[3][4] = sin(angles);
	T[1][0] = T[4][3] = -sin(angles);
	T[2][2] = T[5][5] =  1.0;

	// Set start of timestep value
	T_n = T;

	// Initialise internal forces to zero
	F.resize(DOFs, 0.0);
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
