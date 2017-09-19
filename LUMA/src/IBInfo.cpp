/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
*  Copyright (C) The University of Manchester 2017
*  E-mail contact: info@luma.manchester.ac.uk
*
* This software is for academic use only and not available for
* further distribution commericially or otherwise without written consent.
*
*/

#include "../inc/stdafx.h"
#include "../inc/IBInfo.h"


// *********************** Marker-Owner Comm Methods ***************************

// *****************************************************************************
///	\brief	Default constructor for owner-side marker-owner comm class
markerCommOwnerSideClass::markerCommOwnerSideClass() {

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
markerCommOwnerSideClass::markerCommOwnerSideClass(int rank, int body, int marker) {

	// Set values
	rankComm = rank;
	bodyID = body;
	markerID = marker;
	nSupportSites = 0;
}


// *****************************************************************************
///	\brief	Default constructor for marker-side marker-owner comm class
markerCommMarkerSideClass::markerCommMarkerSideClass() {

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
markerCommMarkerSideClass::markerCommMarkerSideClass(int rank, int body, int idx) {

	// Set values
	rankComm = rank;
	bodyID = body;
	markerIdx = idx;
}



// ********************** Marker-Support Comm Methods **************************



// *****************************************************************************
///	\brief	Default constructor for support-side marker-support comm class
supportCommSupportSideClass::supportCommSupportSideClass() {

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
supportCommSupportSideClass::supportCommSupportSideClass(int rankID, int body, std::vector<int> &position) {

	// Default values
	rankComm = rankID;
	bodyID = body;
	supportIdx.push_back(position[eXDirection]);
	supportIdx.push_back(position[eYDirection]);
	supportIdx.push_back(position[eZDirection]);
}


// *****************************************************************************
///	\brief	Default constructor for marker-side marker-support comm class
supportCommMarkerSideClass::supportCommMarkerSideClass() {

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
supportCommMarkerSideClass::supportCommMarkerSideClass(int rankID, int body, int marker, int support) {

	// Default values
	bodyID = body;
	markerIdx = marker;
	supportID = support;
	rankComm = rankID;
}
