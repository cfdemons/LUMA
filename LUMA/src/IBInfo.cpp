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

///// \brief Custom constructor for different types of message containers.
/////	\param	iBody	pointer to the iBody being packed.
/////	\param	m		marker.
/////	\param	type	the type of container to be created.
//IBInfo::IBInfo(eIBInfoType type, IBBody *iBody, std::vector<int> &buffer) {
//
//	// Get MPI Manager Instance
//	MpiManager *mpim = MpiManager::getInstance();
//
//	switch (type)
//	{
//
//	case eBodyNumMarkers:
//
//		// Set the sending rank and the number of markers it has to send to the owning rank
//		this->sendingRank = mpim->my_rank;
//		this->nMarkers = iBody->markers.size();
//		this->mpi_type = MPI_INT;
//
//		// Fill the vector with required data
//		buffer.push_back(this->sendingRank);
//		buffer.push_back(this->nMarkers);
//
//		break;
//	}
//}
//
///// \brief Custom constructor for different types of message containers.
/////	\param	iBody	pointer to the iBody being packed.
/////	\param	m		marker.
/////	\param	type	the type of container to be created.
//IBInfo::IBInfo(eIBInfoType type, IBBody *iBody, std::vector<double> &buffer)
//{
//
//	switch (type)
//	{
//
//	case eIBDeltaSum:
//
//		// Pack voxel_centres for each marker and sum of delta functions
//		// First support point should be the closest
//		markerX = iBody->markers[m].position[eXDirection];
//		markerY = iBody->markers[m].position[eYDirection];
//		markerZ = iBody->markers[m].position[eZDirection];
//
//		// Delta values
//		for (int i = 0; i < iBody->markers[m].deltaval.size(); i++)
//			deltaVals[i] = iBody->markers[m].deltaval[i];
//
//		break;
//	}
//}
