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
#include "../inc/IBBody.h"
#include "../inc/IBMarker.h"

IBInfo::IBInfo()
{
	// Default constructor
}

/// \brief Custom constructor for different types of message containers.
///	\param	iBody	pointer to the iBody being packed.
///	\param	type	the type fo container to be created.
IBInfo::IBInfo(IBBody *iBody, eIBInfoType type)
{

	switch (type)
	{

	case eIBDeltaSum:
		// Pack voxel_centres for each marker and sum of delta functions
		for (size_t m = 0; m < iBody->markers.size(); ++m)
		{
			// First support point should be the closest
			voxel_centre_X.push_back(iBody->_Owner->XPos[iBody->markers[m].supp_i[0]]);
			voxel_centre_Y.push_back(iBody->_Owner->YPos[iBody->markers[m].supp_j[0]]);
			voxel_centre_Z.push_back(iBody->_Owner->ZPos[iBody->markers[m].supp_k[0]]);

			// Sum of deltas
			delta_sum.push_back(
				std::accumulate(iBody->markers[m].deltaval.begin(), iBody->markers[m].deltaval.end(), 0.0)
				);
		}
		break;

		// Add other cases here

	}
}

///	\brief Maps a version of the IBInfo structure to an MPI_Struct datatype.
///	\param	type	type of container you want to map.
///	\return	handle to the MPI struct data.
int IBInfo::mapToMpiStruct(eIBInfoType type)
{

	switch (type)
	{

	case eIBDeltaSum:

		// TODO: Map to a sturcture
		//return MPI_Type_struct();

		break;
	}

	return 0;
}