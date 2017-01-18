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
	IBInfo(IBBody *iBody, eIBInfoType type);

	// Methods
	int mapToMpiStruct(eIBInfoType type);

private:
	// Possible member data
	std::vector<double> voxel_centre_X;
	std::vector<double> voxel_centre_Y;
	std::vector<double> voxel_centre_Z;
	std::vector<double> delta_sum;

};

#endif	// L_IBINFO_H