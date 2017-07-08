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

// *****************************************************************************
/// \brief	Construct the R vector for all bodies on current level.
///
///
/// \param level current grid level
void ObjectManager::fem_getOffRankForces(int level, std::vector<std::vector<int>> &markerIdx, std::vector<std::vector<std::vector<double>>> &forceBuffer) {

	// Get the mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Perform interpolation communication
	std::vector<std::vector<double>> recvBuffer;
	mpim->mpi_offRankForcesComm(level, recvBuffer);

	// Resize the forceBuffer
	forceBuffer.resize(iBody.size());
	markerIdx.resize(iBody.size());

	// Vector of indices for unpacking receive buffer
	std::vector<int> idx(mpim->num_ranks, 0);

	// Rank to collect data from
	int ib, fromRank, markerID;

	// Vector of values to push back
	std::vector<double> vals(L_DIMS, 0.0);

	// Now loop through all received data
	for (int i = 0; i < mpim->epsCommOwnerSide.size(); i++) {

		// Get body ID
		ib = bodyIDToIdx[mpim->epsCommOwnerSide[i].bodyID];

		// Only pack if body belongs to current grid level
		if (iBody[ib]._Owner->level == level) {

			// Get ID info
			fromRank = mpim->epsCommOwnerSide[i].rankComm;
			markerID = mpim->epsCommOwnerSide[i].markerID;

			// Loop through dimensions
			for (int d = 0; d < L_DIMS; d++)
				vals[d] = recvBuffer[fromRank][idx[fromRank]+d];

			// Increment idx
			idx[fromRank] += L_DIMS;

			// Push the values to end of vector
			forceBuffer[ib].push_back(vals);

			// Add to markerIdx
			markerIdx[ib].push_back(markerID);
		}
	}
}

// ****************************************************************************
