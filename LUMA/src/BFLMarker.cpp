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
#include "../inc/BFLMarker.h"

// Implementation for the BFLMarker class //


/// Default constructor
BFLMarker::BFLMarker(void)
{
}

/// Default destructor
BFLMarker::~BFLMarker(void)
{
}

/// \brief Custom constructor with position.
/// \param x			x-position of marker
/// \param y			y-position of marker
/// \param z			z-position of marker
/// \param markerID		ID of marker within body
///	\param body_owner	Grid on which primary support is to be found.
BFLMarker::BFLMarker(double x, double y, double z, int markerID, GridObj const * const body_owner) : Marker(x, y, z, markerID, body_owner)
{
}
