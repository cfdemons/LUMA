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

#include "../inc/ObjectManager.h"

// ************************************************************************* //
/// \brief	Do communication required for getting off-rank markers
///
///
/// \param current grid level
void MpiManager::mpi_forceCommGather(int level) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare send buffer
	std::vector<std::vector<double>> sendBuffer(num_ranks, std::vector<double>(0));

	// Pack the data to send
	int toRank, ib, m;
	for (int i = 0; i < markerCommMarkerSide[level].size(); i++) {

		// Get body ID
		ib = objman->bodyIDToIdx[markerCommMarkerSide[level][i].bodyID];

		// Only pack if body belongs to current grid level and is flexible
		if (objman->iBody[ib]._Owner->level == level && objman->iBody[ib].isFlexible) {

			// Get ID info
			toRank = markerCommMarkerSide[level][i].rankComm;
			m = markerCommMarkerSide[level][i].markerIdx;

			// Pack marker data
			sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].force_xyz[eXDirection]);
			sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].force_xyz[eYDirection]);
#if (L_DIMS == 3)
			sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].force_xyz[eZDirection]);
#endif
		}
	}

	// Now loop through send buffer and send to appropriate ranks
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (sendBuffer[toRank].size() > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&sendBuffer[toRank].front(), sendBuffer[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Get buffer sizes
	std::vector<int> bufferSize(num_ranks, 0);
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++) {
		if (objman->iBody[objman->bodyIDToIdx[markerCommOwnerSide[level][i].bodyID]].isFlexible)
			bufferSize[markerCommOwnerSide[level][i].rankComm] += L_DIMS;
	}

	// Now create receive buffer
	std::vector<std::vector<double>> recvBuffer(num_ranks, std::vector<double>(0));

	// Now loop through and receive the buffer
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Only do if there is information to receive
		if (bufferSize[fromRank] > 0) {
			recvBuffer[fromRank].resize(bufferSize[fromRank]);
			MPI_Recv(&recvBuffer[fromRank].front(), recvBuffer[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Index vector for looping through recvBuffer
	std::vector<int> idx(num_ranks, 0);

	// Now unpack
	int fromRank;
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++) {

		// Get body idx
		ib = objman->bodyIDToIdx[markerCommOwnerSide[level][i].bodyID];

		// Only unpack if body is flexible
		if (objman->iBody[ib].isFlexible) {

			// Get ID info
			fromRank = markerCommOwnerSide[level][i].rankComm;
			m = markerCommOwnerSide[level][i].markerID;

			// Loop through and set force
			for (int d = 0; d < L_DIMS; d++)
				objman->iBody[ib].markers[m].force_xyz[d] = recvBuffer[fromRank][idx[fromRank]+d];

			// Increment
			idx[fromRank] += L_DIMS;
		}
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}

// ************************************************************************** //
