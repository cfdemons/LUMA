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

#ifndef MARKER_H
#define MARKER_H

#include "stdafx.h"

/// \brief	Generic marker class
///
///			Other markers (e.g. IBMarker) inherit from this.
class Marker
{

public:

	/************** Constructors **************/

	// Default constructor
	Marker(void)
	{
		// Set marker position to zero by default
		position.push_back(0.0);
		position.push_back(0.0);
		position.push_back(0.0);

		// Set ID to zero by default
		id = 0;
	};

	// Default destructor
	~Marker(void)
	{
	};

	/// \brief Custom constructor which locates marker
	///
	///			In order to properly initialise during construction, a
	///			grid should be passed in on which primary support point
	///			can be found.
	///
	/// \param	x			X-position of marker
	/// \param	y			Y-position of marker
	/// \param	z			Z-position of marker
	/// \param	markerID	ID of marker within body
	///	\param	body_owner	Grid on which primary support point is to be found.
	Marker(double x, double y, double z, int markerID, GridObj const * const body_owner)
	{

		// Set marker position
		position.push_back(x);
		position.push_back(y);
		position.push_back(z);

		// Get the i, j, and k for the closest voxel to markers
		std::vector<int> ijk;
		eLocationOnRank loc = eNone;
		bool onRank = GridUtils::isOnThisRank(x, y, z, &loc, body_owner, &ijk);

		// Only add the support point information if the marker is on the current rank core
		if (onRank && loc == eCore)
		{
			supp_i.push_back(ijk[0]);
			supp_j.push_back(ijk[1]);
			supp_k.push_back(ijk[2]);

			// Get the position of the support point
			supp_x.push_back(body_owner->XPos[ijk[eXDirection]]);
			supp_y.push_back(body_owner->YPos[ijk[eYDirection]]);
			supp_z.push_back(body_owner->ZPos[ijk[eZDirection]]);

			// Set rank of this centre support point and marker ID in body
			support_rank.push_back(GridUtils::safeGetRank());
		}
		id = markerID;
	}


public:

	/************** Member Data **************/

	int id;			///< ID of marker within its owning body

	std::vector<double> position;	///< Position vector of marker location in physical units

	/* Vector of indices for lattice sites considered to be in support of the marker:
	* In IBM these are all the lattice support sites for spreading and interpolating;
	* In BFL these are the voxel indices in which the BFL marker resides;
	*/
	std::vector<int> supp_i;	///< X-indices of lattice sites in support of this marker
	std::vector<int> supp_j;	///< Y-indices of lattice sites in support of this marker
	std::vector<int> supp_k;	///< Z-indices of lattice sites in support of this marker

	std::vector<double> supp_x;	///< X-position of lattice sites in support of this marker
	std::vector<double> supp_y;	///< Y-position of lattice sites in support of this marker
	std::vector<double> supp_z;	///< Z-position of lattice sites in support of this marker

	std::vector<int> support_rank; ///< Array of indices indicating on which rank the given support point resides
};

#endif
