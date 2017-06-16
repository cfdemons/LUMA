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

// *****************************************************************************

#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"
#include "../inc/IBInfo.h"


/// \brief	Builds MPI comm class for MPI communication in IBM.
void ObjectManager::ibm_buildMPIComms() {

	// Get MPI manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Build MPI comm classes for markers
	mpim->mpi_getIBMarkers();

	// Build MPI comm classes for epsilon calculation
	mpim->mpi_buildEpsComms();

	// Build MPI comm class for support communication
	mpim->mpi_buildSupportComms();
}


/// \brief	Pass velocity values from support site which exist off-rank.
void ObjectManager::ibm_interpolate_comm() {

//	// Get the mpi manager instance
//	MpiManager *mpim = MpiManager::getInstance();
//
//	// Pack and send data first
//	int sendToRank, ib;
//	std::vector<int> idx(3, 0);
//	std::vector<std::vector<double>> sendVel;
//	sendVel.resize(mpim->num_ranks);
//	for (int i = 0; i < supportCommSupportSide.size(); i++) {
//
//		// Get rank to send to
//		sendToRank = supportCommSupportSide[i].rank;
//		ib = supportCommSupportSide[i].bodyID;
//
//		// Get grid sizes
//		size_t M_lim = iBody[ib]._Owner->M_lim;
//#if (L_DIMS == 3)
//		size_t K_lim = iBody[ib]._Owner->K_lim;
//#endif
//
//		// Get indices
//		idx = supportCommSupportSide[i].supportIdx;
//
//		// Pack velocity into buffer
//#if (L_DIMS == 2)
//		sendVel[sendToRank].push_back(iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], eXDirection, M_lim, L_DIMS));
//		sendVel[sendToRank].push_back(iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], eYDirection, M_lim, L_DIMS));
//#elif (L_DIMS == 3)
//		sendVel[sendToRank].push_back(iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eXDirection, M_lim, K_lim, L_DIMS));
//		sendVel[sendToRank].push_back(iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eYDirection, M_lim, K_lim, L_DIMS));
//		sendVel[sendToRank].push_back(iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eZDirection, M_lim, K_lim, L_DIMS));
//#endif
//	}
//
//	// Set array of request values
//	MPI_Request send_requests[mpim->num_ranks] = {MPI_REQUEST_NULL};
//
//	// Now send the values
//	for (int toRank = 0; toRank < mpim->num_ranks; toRank++) {
//
//		// Check if current rank needs to send stuff to rank toRank
//		if (sendVel[toRank].size() > 0)
//			MPI_Isend(&sendVel[toRank].front(), sendVel[toRank].size(), MPI_DOUBLE, toRank, mpim->my_rank, mpim->world_comm, &send_requests[toRank]);
//	}
//
//	// Get receive sizes from each rank
//	std::vector<int> recvVelSize(mpim->num_ranks, 0);
//	for (int i = 0; i < supportCommMarkerSide.size(); i++) {
//		recvVelSize[supportCommMarkerSide[i].rank] += L_DIMS;
//	}
//
//	// Now receive the velocity values
//	std::vector<std::vector<double>> recvVel;
//	recvVel.resize(mpim->num_ranks);
//	for (int fromRank = 0; fromRank < mpim->num_ranks; fromRank++) {
//
//		// If there is data to receive from rank fromRank
//		if (recvVelSize[fromRank] > 0) {
//			recvVel[fromRank].resize(recvVelSize[fromRank]);
//			MPI_Recv(&recvVel[fromRank].front(), recvVel[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, mpim->world_comm, MPI_STATUS_IGNORE);
//		}
//	}
//
//	// Now interpolate these remaining values onto the marker
//	int m, s, rank;
//	for (int i = 0; i < supportCommMarkerSide.size(); i++) {
//
//		// Get IDs of support site
//		ib = supportCommMarkerSide[i].bodyID;
//		m = supportCommMarkerSide[i].markerID;
//		s = supportCommMarkerSide[i].supportID;
//		rank = supportCommMarkerSide[i].rank;
//
//		// Interpolate these values
//		for (int dir = 0; dir < L_DIMS; dir++)
//			iBody[ib].markers[m].interpVel[dir] += recvVel[rank][dir] * iBody[ib].markers[m].deltaval[s] * iBody[ib].markers[m].local_area;
//
//		// Erase those values
//		recvVel[rank].erase(recvVel[rank].begin(), recvVel[rank].begin()+L_DIMS);
//	}
}

// ****************************************************************************
