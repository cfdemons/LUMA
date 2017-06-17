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
void ObjectManager::ibm_gatherOffRankVels(int level) {

	// Get the mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Perform interpolation communication
	std::vector<std::vector<double>> interpVels;
	mpim->mpi_interpolateComm(level, interpVels);

	// Create idx vector
	std::vector<int> idx(mpim->num_ranks, 0);

	// Now interpolate these remaining values onto the marker
	int ib, m, s, fromRank;
	for (int i = 0; i < mpim->supportCommMarkerSide.size(); i++) {

		// Get body idx
		ib = bodyIDToIdx[mpim->supportCommMarkerSide[i].bodyID];

		// Only do if body is on this grid level
		if (iBody[ib]._Owner->level == level) {

			// Get IDs of support site
			m = mpim->supportCommMarkerSide[i].markerIdx;
			s = mpim->supportCommMarkerSide[i].supportID;
			fromRank = mpim->supportCommMarkerSide[i].rankComm;

			// Interpolate these values
			for (int dir = 0; dir < L_DIMS; dir++)
				iBody[ib].markers[m].interpVel[dir] += interpVels[fromRank][dir+idx[fromRank]] * iBody[ib].markers[m].deltaval[s] * iBody[ib].markers[m].local_area;

			// Shift index
			idx[fromRank] += L_DIMS;
		}
	}
}

// ****************************************************************************
