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
#ifndef IBINFO_H
#define IBINFO_H

#include "stdafx.h"


/// \brief Class for arranging data for epsilon communication on owner side.
///
///
class markerCommOwnerSideClass {

	// Make MPIManager a friend so it can access body data
	friend class MpiManager;
	friend class IBBody;
	friend class ObjectManager;

public:

	// Default Constructor
	markerCommOwnerSideClass();

	// Custom constructor for creating epsCommOwnerSideClass object
	markerCommOwnerSideClass(int rank, int body, int marker);

private:

	// Member data
	int rankComm;
	int bodyID;
	int markerID;
	int nSupportSites;
};


/// \brief Class for arranging data for epsilon communication on marker side.
///
///
class markerCommMarkerSideClass {

	// Make MPIManager a friend so it can access body data
	friend class MpiManager;

public:

	// Default Constructor
	markerCommMarkerSideClass();

	// Custom constructor for creating epsCommMarkerSideClass object
	markerCommMarkerSideClass(int rank, int body, int idx);

private:

	// Member data
	int rankComm;
	int bodyID;
	int markerIdx;
};



///	\brief	Class for arranging data for support-marker communication on support side.
///
///
class supportCommSupportSideClass {

public:

	// Default Constructor
	supportCommSupportSideClass();

	// Custom constructor for creating supportCommSupportSide object
	supportCommSupportSideClass(int rankID, int bodyID, std::vector<int> &idx);

public:

	// ID data
	int rankComm;
	int bodyID;

	// Support index
	std::vector<int> supportIdx;
};


///	\brief	Class for arranging data for support-marker communication on marker side.
///
///
class supportCommMarkerSideClass {

public:

	// Default Constructor
	supportCommMarkerSideClass();

	// Custom constructor for creating supportCommMarkerSide object
	supportCommMarkerSideClass(int body, int marker, int support, int rankID);

public:

	// ID data
	int rankComm;
	int bodyID;
	int markerIdx;
	int supportID;
};
#endif	// L_IBINFO_H
