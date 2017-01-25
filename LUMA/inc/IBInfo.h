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

///	\enum eIBInfoType
///	\brief	Type of container required.
enum eIBInfoType {
	eIBDeltaSum,
	eIBEpsilon,
	eIBVelocityInterpolation,
	eIBVelocitySpreading,
	eIBMarkerPositions
};


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
	IBInfo(IBBody *iBody, int m, eIBInfoType type);

	// Methods
	int mapToMpiStruct(eIBInfoType type, MPI_Datatype *mpi_struct_type);

public:

	// Possible member data
	double markerX;
	double markerY;
	double markerZ;
	//std::vector<double> deltaVals;	// ** Could solve this using arrays with a buffer size that will always have enough space for information (plus a value telling you how many support points there are)

};

#endif	// L_IBINFO_H
