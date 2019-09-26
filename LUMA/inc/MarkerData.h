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

#ifndef MARDAT_H
#define MARDAT_H

#include "stdafx.h"

/// \brief	Container class to hold marker information.
class MarkerData {

public:

	/// \brief Constructor.
	/// \param i i-index of primary support site
	/// \param j j-index of primary support site
	/// \param k k-index of primary support site
	/// \param x x-position of marker
	/// \param y y-position of marker
	/// \param z z-position of marker
	/// \param ID marker number in a given body
	MarkerData(int i, int j, int k, double x, double y, double z, int ID)
		: i(i), j(j), k(k), x(x), y(y), z(z), ID(ID)
	{ };

	/// \brief	Default Constructor.
	///
	///			Initialise with invalid marker indicator which 
	///			is to set the x position to NaN.
	MarkerData(void) {

		// Give invalid ID number at construction time
		this->ID = -1;

	};

	/// Default destructor
	~MarkerData(void) {};

	/// Method to check if marker data structure is valid
	bool isValid()
	{
		if (ID == -1) return false;
		return true;
	}

	// Voxel indices
	int i;	///< i-index of primary support site
	int j;	///< j-index of primary support site
	int k;	///< k-index of primary support site

	/// Marker ID (position in array of markers)
	int ID;

	// Marker position
	double x;	///< x-position of marker
	double y;	///< y-position of marker
	double z;	///< z-position of marker

};

#endif // MARDAT_H