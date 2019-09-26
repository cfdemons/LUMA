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
