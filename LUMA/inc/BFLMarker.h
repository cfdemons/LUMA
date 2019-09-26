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

#ifndef BFLMARKER_H
#define BFLMARKER_H

#include "stdafx.h"
#include "Marker.h"

/// \brief	BFL marker.
///
///			This class declaration is for a BFL Lagrange point.
///			A collection of these points form BFL body.
class BFLMarker :
	public Marker
{

	// Make ObjectManager a friend class so it can access the protected data of BFLMarker objects
	friend class ObjectManager;

	// Allow BFLbody to access the marker positions etc.
	friend class BFLBody;

public:
	// Default constructor and destructor
	BFLMarker(void);
	~BFLMarker(void);

	// Custom constructor when positions are passed
	BFLMarker(double x, double y, double z, int markerID, GridObj const * const body_owner);

protected:
	// Per marker force storage
	double forceX = 0.0;	///< Instantaneous X-direction force on marker
	double forceY = 0.0;	///< Instantaneous Y-direction force on marker
	double forceZ = 0.0;	///< Instantaneous Z-direction force on marker

};

#endif
