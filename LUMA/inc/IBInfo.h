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

// Forward declare
class IBBody;

///	\brief	Structure for passing IB information between MPI processes.
///
///			This structure has a series of different constructors depending on 
///			what information should be passed.
///
class IBInfo
{

public:
	// Default Constructor
	IBInfo();
	IBInfo(eIBInfoType type, IBBody *iBody, std::vector<int> &buffer);
	IBInfo(eIBInfoType type, IBBody *iBody, std::vector<double> &buffer);

	// Methods


public:

	// Needed for all cases
	MPI_Datatype mpi_type;

	// Data needed for sending to body owner
	int sendingRank;
	int nMarkers;

	// Data needed for epsilon calculation
	double markerX;
	double markerY;
	double markerZ;
	double deltaVals[L_SUPPORT_SIZE];
};

#endif	// L_IBINFO_H
