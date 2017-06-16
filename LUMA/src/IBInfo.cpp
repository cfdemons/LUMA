/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
*  Copyright (C) 2015, 2016
*  E-mail contact: info@luma.manchester.ac.uk
*
* This software is for academic use only and not available for
* distribution without written consent.
*
*/

#include "../inc/stdafx.h"
#include "../inc/IBInfo.h"

// Default constructor
epsCommOwnerSideClass::epsCommOwnerSideClass() {

	// Set default values
	rankComm = 0;
	bodyID = 0;
	markerID = 0;
	nSupportSites = 0;
}

/// \brief 	Custom constructor for creating epsCommOwnerSideClass object
///	\param	rank				rank receiving from.
///	\param	body				body ID to receive.
///	\param	marker				marker ID to receive.
epsCommOwnerSideClass::epsCommOwnerSideClass(int rank, int body, int marker) {

	// Set values
	rankComm = rank;
	bodyID = body;
	markerID = marker;
	nSupportSites = 0;
}


// Default constructor
epsCommMarkerSideClass::epsCommMarkerSideClass() {

	// Set default values
	rankComm = 0;
	bodyID = 0;
	markerIdx = 0;
}

/// \brief 	Custom constructor for creating epsCommMarkerSideClass object
///	\param	rank				rank sending to.
///	\param	body				body ID to receive.
///	\param	marker				marker ID to receive.
epsCommMarkerSideClass::epsCommMarkerSideClass(int rank, int body, int idx) {

	// Set values
	rankComm = rank;
	bodyID = body;
	markerIdx = idx;
}




/// ******************** ///



// Default constructor
supportCommSupportSideClass::supportCommSupportSideClass() {

	// Default values
	rankComm = 0;
	bodyID = 0;
}

/// \brief Custom constructor for support side communication class.
///	\param	rank to send to		index of body.
///	\param	marker				index of marker.
///	\param	support				index of support.
supportCommSupportSideClass::supportCommSupportSideClass(int rankID, int body, std::vector<int> &position) {

	// Default values
	rankComm = rankID;
	bodyID = body;
	supportIdx.push_back(position[eXDirection]);
	supportIdx.push_back(position[eYDirection]);
	supportIdx.push_back(position[eZDirection]);
}


// Default constructor
supportCommMarkerSideClass::supportCommMarkerSideClass() {

	// Default values
	bodyID = 0;
	markerIdx = 0;
	supportID = 0;
	rankComm = 0;
}

/// \brief Custom constructor for support receiver communication class.
///	\param	body				index of body.
///	\param	marker				index of marker.
///	\param	support				index of support.
supportCommMarkerSideClass::supportCommMarkerSideClass(int body, int marker, int support, int rankID) {

	// Default values
	bodyID = body;
	markerIdx = marker;
	supportID = support;
	rankComm = rankID;
}
