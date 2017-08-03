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


/// \brief Owner-side comm class for marker-owner communications
///
///			Class for arranging marker-owner MPI communication.
class markerCommOwnerSideClass {

	/************** Friends **************/
	friend class MpiManager;
	friend class IBBody;
	friend class ObjectManager;

public:

	/************** Constructors **************/
	markerCommOwnerSideClass();

	// Custom constructor for creating epsCommOwnerSideClass object
	markerCommOwnerSideClass(int rank, int body, int marker);

private:

	/************** Member Data **************/
	int rankComm;				///< Rank to be communicating with
	int bodyID;					///< Global body ID for communicating
	int markerID;				///< Marker ID for communicating
	int nSupportSites;			///< Number of support sites that need to be communicated
};


/// \brief Marker-side comm class for marker-owner communications
///
///			Class for arranging marker-owner MPI communication.
class markerCommMarkerSideClass {

	/************** Friends **************/
	friend class MpiManager;

public:

	/************** Constructors **************/
	markerCommMarkerSideClass();

	// Custom constructor for creating epsCommMarkerSideClass object
	markerCommMarkerSideClass(int rank, int body, int idx);

private:

	/************** Member Data **************/
	int rankComm;				///< Rank to be communicating with
	int bodyID;					///< Global body ID for communicating
	int markerIdx;				///< Local (rank) index of marker for communicating
};



/// \brief Support-side comm class for marker-support communications
///
///			Class for arranging marker-support MPI communication.
class supportCommSupportSideClass {

public:

	/************** Constructors **************/
	supportCommSupportSideClass();

	// Custom constructor for creating supportCommSupportSide object
	supportCommSupportSideClass(int rankID, int bodyID, std::vector<int> &idx);

public:

	/************** Member Data **************/
	int rankComm;						///< Rank to be communicating with
	int bodyID;							///< Global body ID for communicating
	std::vector<int> supportIdx;		///< Local (rank) grid indices of support point
};


/// \brief Marker-side comm class for marker-support communications
///
///			Class for arranging marker-support MPI communication.
class supportCommMarkerSideClass {

public:

	/************** Constructors **************/
	supportCommMarkerSideClass();

	// Custom constructor for creating supportCommMarkerSide object
	supportCommMarkerSideClass(int rankID, int body, int marker, int support);

public:

	/************** Member Data **************/
	int rankComm;					///< Rank to be communicating with
	int bodyID;						///< Global body ID for communicating
	int markerIdx;					///< Local (rank) index of marker for communicating
	int supportID;					///< Support index within the marker
};
#endif	// L_IBINFO_H
