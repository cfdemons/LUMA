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
#include "../inc/IBInfo.h"


// *********************** Marker-Owner Comm Methods ***************************

// *****************************************************************************
///	\brief	Default constructor for owner-side marker-owner comm class
MarkerCommOwnerSideClass::MarkerCommOwnerSideClass() {

	// Set default values
	rankComm = 0;
	bodyID = 0;
	markerID = 0;
	nSupportSites = 0;
}


// *****************************************************************************
///	\brief	Custom constructor for owner-side marker-owner comm class
///
///	\param	rank	rank to communicate with
///	\param	body	global body ID
///	\param	marker	marker ID within body
MarkerCommOwnerSideClass::MarkerCommOwnerSideClass(int rank, int body, int marker) {

	// Set values
	rankComm = rank;
	bodyID = body;
	markerID = marker;
	nSupportSites = 0;
}


// *****************************************************************************
///	\brief	Default constructor for marker-side marker-owner comm class
MarkerCommMarkerSideClass::MarkerCommMarkerSideClass() {

	// Set default values
	rankComm = 0;
	bodyID = 0;
	markerIdx = 0;
}

// *****************************************************************************
///	\brief	Custom constructor for marker-side marker-owner comm class
///
///	\param	rank	rank to communicate with
///	\param	body	global body ID
///	\param	idx		marker local index within body
MarkerCommMarkerSideClass::MarkerCommMarkerSideClass(int rank, int body, int idx) {

	// Set values
	rankComm = rank;
	bodyID = body;
	markerIdx = idx;
}



// ********************** Marker-Support Comm Methods **************************



// *****************************************************************************
///	\brief	Default constructor for support-side marker-support comm class
SupportCommSupportSideClass::SupportCommSupportSideClass() {

	// Default values
	rankComm = 0;
	bodyID = 0;
}

// *****************************************************************************
///	\brief	Custom constructor for support-side marker-support comm class
///
///	\param	rankID		rank to communicate with
///	\param	body		global body ID
///	\param	position	position indices of support point
SupportCommSupportSideClass::SupportCommSupportSideClass(int rankID, int body, std::vector<int> &position) {

	// Default values
	rankComm = rankID;
	bodyID = body;
	supportIdx.push_back(position[eXDirection]);
	supportIdx.push_back(position[eYDirection]);
	supportIdx.push_back(position[eZDirection]);
}


// *****************************************************************************
///	\brief	Default constructor for marker-side marker-support comm class
SupportCommMarkerSideClass::SupportCommMarkerSideClass() {

	// Default values
	bodyID = 0;
	markerIdx = 0;
	supportID = 0;
	rankComm = 0;
}


// *****************************************************************************
///	\brief	Custom constructor for marker-side marker-support comm class
///
///	\param	rankID		rank to communicate with
///	\param	body		global body ID
///	\param	marker		marker ID within body
///	\param	support		support ID within marker
SupportCommMarkerSideClass::SupportCommMarkerSideClass(int rankID, int body, int marker, int support) {

	// Default values
	bodyID = body;
	markerIdx = marker;
	supportID = support;
	rankComm = rankID;
}
