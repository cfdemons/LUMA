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


///	\brief	Class for arranging data for support-marker communication on support side.
///
///
class supportCommSupportSideClass {

public:

	// Default Constructor
	supportCommSupportSideClass();

	// Custom constructor for creating supportCommSupportSide object
	supportCommSupportSideClass(int rankID, int bodyID, std::vector<int> &position);

public:

	// ID data
	int rank;
	int bodyID;

	// Force data
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
	int rank;
	int bodyID;
	int markerID;
	int supportID;
};



// ***** USE IBODY INSTEAD WITH A CUSTOM CONSTRUCTOR FOR CREATING THE BODIES **** //
///	\brief	Class for arranging data before epsilon calculation.
///
///
class epsCalcMarkerClass
{

public:
	// Default Constructor
	epsCalcMarkerClass();

	// Custom constructor for creating eps calc marker
	epsCalcMarkerClass(int bodyID, std::vector<double> position, double area, double dilation, std::vector<std::vector<double>> supp_position, std::vector<double> deltaval);

public:

	// Marker data
	int bodyID;
	std::vector<double> position;
	double local_area;
	double dilation;

	// Support data
	std::vector<double> deltaval;
	std::vector<double> supp_x;
	std::vector<double> supp_y;
	std::vector<double> supp_z;
};

#endif	// L_IBINFO_H
