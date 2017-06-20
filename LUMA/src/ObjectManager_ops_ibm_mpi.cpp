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
void ObjectManager::ibm_interpolateOffRankVels(int level) {

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

			// Interpolate density
			iBody[ib].markers[m].interpRho += interpVels[fromRank][idx[fromRank]] * iBody[ib].markers[m].deltaval[s] * iBody[ib].markers[m].local_area;

			// Shift index
			idx[fromRank]++;

			// Interpolate these values
			for (int dir = 0; dir < L_DIMS; dir++)
				iBody[ib].markers[m].interpMom[dir] += interpVels[fromRank][dir+idx[fromRank]] * iBody[ib].markers[m].deltaval[s] * iBody[ib].markers[m].local_area;

			// Shift index
			idx[fromRank] += L_DIMS;
		}
	}
}

/// \brief	Pass velocity values from support site which exist off-rank.
void ObjectManager::ibm_spreadOffRankForces(int level) {

	// Get the mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Perform interpolation communication
	std::vector<std::vector<double>> spreadForces;
	mpim->mpi_spreadComm(level, spreadForces);

	// Create idx vector
	std::vector<int> idx(mpim->num_ranks, 0);
	std::vector<int> suppIdx(3, 0);

	// Now interpolate these remaining values onto the marker
	int ib, fromRank;
	for (int i = 0; i < mpim->supportCommSupportSide.size(); i++) {

		// Get body idx
		ib = bodyIDToIdx[mpim->supportCommSupportSide[i].bodyID];

		// Only do if body is on this grid level
		if (iBody[ib]._Owner->level == level) {

			// Get grid sizes
			size_t M_lim = iBody[ib]._Owner->M_lim;
			size_t K_lim = iBody[ib]._Owner->K_lim;

			// Get IDs of support site
			fromRank = mpim->supportCommSupportSide[i].rankComm;
			suppIdx = mpim->supportCommSupportSide[i].supportIdx;

			// Interpolate these values
			for (int dir = 0; dir < L_DIMS; dir++)
				iBody[ib]._Owner->force_xyz(suppIdx[eXDirection], suppIdx[eYDirection], suppIdx[eZDirection], dir, M_lim, K_lim, L_DIMS) -=
						spreadForces[fromRank][idx[fromRank] + dir];

			// Shift index
			idx[fromRank] += L_DIMS;
		}
	}
}

// ****************************************************************************
