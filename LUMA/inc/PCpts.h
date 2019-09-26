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

#ifndef PCPTS_H
#define PCPTS_H

#include "stdafx.h"
/// \brief	Class to hold point cloud data.
///
///			A container class for hold the X, Y and Z positions of 
///			points in a point cloud.
class PCpts {

public:

	/// Default constructor
	PCpts(void) {};

	/// Default destructor
	~PCpts(void) {};

	std::vector<double> x;		///< Vector of X positions
	std::vector<double> y;		///< Vector of Y positions
	std::vector<double> z;		///< Vector of Z positions
	std::vector<int> id;		///< Vector of point IDs
};

#endif