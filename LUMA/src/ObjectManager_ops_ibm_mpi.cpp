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


/// \brief	Gathers data required for epsilon calculation and arranges.
/// \param	rootRank			root rank which holds epsilon values.
/// \param	nMarkersOnThisRank	number of markers on current rank.
/// \param	nMarkersOnAllRanks	vector containing number of markers on each rank for all ranks.
/// \param	markerDisps			displacement vector for MPI communication.
/// \param	markerData			class holding all data required for epsilon calculation.
void ObjectManager::ibm_gatherForEpsCalc(int rootRank, int &nMarkersOnThisRank, std::vector<int> &nMarkersOnAllRanks, std::vector<int> &markerDisps, std::vector<epsCalcMarkerClass> &markerData) {

#ifndef L_BUILD_FOR_MPI

	// Fixed size inputs for markerData constructor
	double area, dilation;
	std::vector<double> position(L_DIMS, 0.0);

	// Variable sized inputs for markerData constructor
	std::vector<std::vector<double>> suppPosition;
	std::vector<double> deltaVal;

	// Loop through all bodies
	for (int ib = 0; ib < iBody.size(); ib++) {

		// Loop through all markers within each body
		for (int node = 0; node < iBody[ib].markers.size(); node++) {

			// Resize support vectors
			suppPosition.resize(iBody[ib].markers[node].deltaval.size(), std::vector<double>(L_DIMS, 0.0));
			deltaVal.resize(iBody[ib].markers[node].deltaval.size(), 0.0);

			// Populate marker data
			area = iBody[ib].markers[node].local_area;
			dilation = iBody[ib].markers[node].dilation;
			position[eXDirection] = iBody[ib].markers[node].position[eXDirection];
			position[eYDirection] = iBody[ib].markers[node].position[eYDirection];
#if (L_DIMS == 3)
			position[eZDirection] = iBody[ib].markers[node].position[eZDirection];
#endif

			// Populate support data
			for (int supp = 0; supp < iBody[ib].markers[node].deltaval.size(); supp++) {
				deltaVal[supp] = iBody[ib].markers[node].deltaval[supp];
				suppPosition[supp][eXDirection] = iBody[ib].markers[node].supp_x[supp];
				suppPosition[supp][eYDirection] = iBody[ib].markers[node].supp_y[supp];
#if (L_DIMS == 3)
				suppPosition[supp][eZDirection] = iBody[ib].markers[node].supp_z[supp];
#endif
			}

			// Insert marker data
			markerData.emplace_back(ib, position, area, dilation, suppPosition, deltaVal);
		}
	}

#else

	// Get mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// GET NUMBER OF MARKERS (AND BODYID) AND SUPPORT SITES TO RECEIVE
	std::vector<int> deltasPerMarkerOnThisRank;
	std::vector<int> bodyIDForMarkerOnThisRank;
	for (int ib = 0; ib < iBody.size(); ib++) {

		// Get number of markers on this rank
		nMarkersOnThisRank += iBody[ib].markers.size();

		// Get number of support sites to send
		for (int node = 0; node < iBody[ib].markers.size(); node++) {
			deltasPerMarkerOnThisRank.push_back(iBody[ib].markers[node].deltaval.size());
			bodyIDForMarkerOnThisRank.push_back(ib);
		}
	}

	// Root rank needs to know how many markers to receive from each rank
	if (mpim->my_rank == rootRank)
		nMarkersOnAllRanks.resize(mpim->num_ranks, 0);

	// Perform gather so that root rank knows how many markers to receive
	MPI_Gather(&nMarkersOnThisRank, 1, MPI_INT, &nMarkersOnAllRanks.front(), 1, MPI_INT, rootRank, mpim->world_comm);

	// Get total number of markers across all ranks
	int nTotalMarkers = std::accumulate(nMarkersOnAllRanks.begin(), nMarkersOnAllRanks.end(), 0);

	// Root rank create and resize result and displacement arrays
	std::vector<int> deltasPerMarkerOnAllRanks;
	std::vector<int> bodyIDForMarkerOnAllRanks;

	// If root rank resize the vectors
	if (mpim->my_rank == rootRank) {

		// Resize the results vector and displacement vector
		deltasPerMarkerOnAllRanks.resize(nTotalMarkers, 0);
		bodyIDForMarkerOnAllRanks.resize(nTotalMarkers, 0);
		markerDisps.resize(mpim->num_ranks, 0);

		// Calculate the block displacements
		for (int i = 1; i < mpim->num_ranks; i++)
			markerDisps[i] = std::accumulate(nMarkersOnAllRanks.begin(), nMarkersOnAllRanks.begin()+i, 0);
	}

	// Perform gather so that root rank knows how many support sites to receive and the body ID for each marker
	MPI_Gatherv(&deltasPerMarkerOnThisRank.front(), nMarkersOnThisRank, MPI_INT, &deltasPerMarkerOnAllRanks.front(), &nMarkersOnAllRanks.front(), &markerDisps.front(), MPI_INT, rootRank, mpim->world_comm);
	MPI_Gatherv(&bodyIDForMarkerOnThisRank.front(), nMarkersOnThisRank, MPI_INT, &bodyIDForMarkerOnAllRanks.front(), &nMarkersOnAllRanks.front(), &markerDisps.front(), MPI_INT, rootRank, mpim->world_comm);


	// NOW CALCULATE EXACT NUMBER OF DATA TO RECEIVE AND FROM WHERE
	std::vector<int> recvBufferSizes;
	std::vector<int> recvBufferDisps;
	std::vector<double> recvBuffer;
	if (mpim->my_rank == rootRank) {

		// Resize for number of ranks
		recvBufferSizes.resize(mpim->num_ranks, 0);
		recvBufferDisps.resize(mpim->num_ranks, 0);

		// Counter
		int count = 0;

		// Loop through all ranks
		for (int rank = 0; rank < mpim->num_ranks; rank++) {

			// Loop through all markers which exist on certain rank
			for (int node = 0; node < nMarkersOnAllRanks[rank]; node++) {

				// Get buffer sizes for gatherv later
				recvBufferSizes[rank] += (L_DIMS + 2) + deltasPerMarkerOnAllRanks[count] * (L_DIMS + 1);

				// Increment counter
				count++;
			}
		}

		// Calculate the block displacements
		for (int i = 1; i < mpim->num_ranks; i++)
			recvBufferDisps[i] = std::accumulate(recvBufferSizes.begin(), recvBufferSizes.begin()+i, 0);

		// Resize receive buffer
		recvBuffer.resize(std::accumulate(recvBufferSizes.begin(), recvBufferSizes.end(), 0), 0.0);
	}


	// NOW PACK THE DATA
	std::vector<double> sendBuffer;
	for (int ib = 0; ib < iBody.size(); ib++) {
		for (int node = 0; node < iBody[ib].markers.size(); node++) {

			// Pack the position first then area and dilation
			sendBuffer.push_back(iBody[ib].markers[node].position[eXDirection]);
			sendBuffer.push_back(iBody[ib].markers[node].position[eYDirection]);
#if (L_DIMS == 3)
			sendBuffer.push_back(iBody[ib].markers[node].position[eZDirection]);
#endif
			sendBuffer.push_back(iBody[ib].markers[node].local_area);
			sendBuffer.push_back(iBody[ib].markers[node].dilation);

			// Now pack the support site data
			for (int supp = 0; supp < iBody[ib].markers[node].deltaval.size(); supp++) {
				sendBuffer.push_back(iBody[ib].markers[node].supp_x[supp]);
				sendBuffer.push_back(iBody[ib].markers[node].supp_y[supp]);
#if (L_DIMS == 3)
				sendBuffer.push_back(iBody[ib].markers[node].supp_z[supp]);
#endif
				sendBuffer.push_back(iBody[ib].markers[node].deltaval[supp]);
			}
		}
	}


	// NOW SEND ALL DATA TO ROOT
	MPI_Gatherv(&sendBuffer.front(), sendBuffer.size(), MPI_DOUBLE, &recvBuffer.front(), &recvBufferSizes.front(), &recvBufferDisps.front(), MPI_DOUBLE, rootRank, mpim->world_comm);


	// NOW UNPACK AND CONSTRUCT MARKERDATA VECTOR
	if (mpim->my_rank == rootRank) {

		// Fixed size inputs for markerData constructor
		double area, dilation;
		std::vector<double> position(L_DIMS, 0.0);

		// Variable sized inputs for markerData constructor
		std::vector<std::vector<double>> suppPosition;
		std::vector<double> deltaVal;

		// Counter
		int count = 0;

		// Loop through all markers
		for (int node = 0; node < nTotalMarkers; node++) {

			// Resize support vectors
			suppPosition.resize(deltasPerMarkerOnAllRanks[node], std::vector<double>(L_DIMS, 0.0));
			deltaVal.resize(deltasPerMarkerOnAllRanks[node], 0.0);

			// Populate marker data
			position[eXDirection] = recvBuffer[count]; count++;
			position[eYDirection] = recvBuffer[count]; count++;
#if (L_DIMS == 3)
			position[eZDirection] = recvBuffer[count]; count++;
#endif
			area = recvBuffer[count]; count++;
			dilation = recvBuffer[count]; count++;

			// Populate support data
			for (int supp = 0; supp < deltasPerMarkerOnAllRanks[node]; supp++) {
				suppPosition[supp][eXDirection] = recvBuffer[count]; count++;
				suppPosition[supp][eYDirection] = recvBuffer[count]; count++;
#if (L_DIMS == 3)
				suppPosition[supp][eZDirection] = recvBuffer[count]; count++;
#endif
				deltaVal[supp] = recvBuffer[count]; count++;
			}

			// Insert marker data
			markerData.emplace_back(bodyIDForMarkerOnAllRanks[node], position, area, dilation, suppPosition, deltaVal);
		}
	}
#endif
}

/// \brief	Scatters epsilon back to correct ranks and unpacks it.
/// \param	rootRank			root rank which holds epsilon values.
/// \param	nMarkersOnThisRank	number of markers on current rank.
/// \param	nMarkersOnAllRanks	vector containing number of markers on each rank for all ranks.
/// \param	markerDisps			displacement vector for MPI communication.
/// \param	epsilon				epsilon values which have been calculated.
void ObjectManager::ibm_scatterAfterEpsCalc(int rootRank, int &nMarkersOnThisRank, std::vector<int> &nMarkersOnAllRanks, std::vector<int> &markerDisps, std::vector<double> &epsilon) {

#ifndef L_BUILD_FOR_MPI

	// Loop through all bodies and markers
	int count = 0;
	for (int ib = 0; ib < iBody.size(); ib++) {
		for (int node = 0; node < iBody[ib].markers.size(); node++) {
			iBody[ib].markers[node].epsilon = epsilon[count];
			count++;
		}
	}

#else

	// Get mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Create receive buffer
	std::vector<double> recvBuffer(nMarkersOnThisRank, 0.0);

	// Scatter the epsilon values to the correct rank
	MPI_Scatterv(&epsilon.front(), &nMarkersOnAllRanks.front(), &markerDisps.front(), MPI_DOUBLE, &recvBuffer.front(), nMarkersOnThisRank, MPI_DOUBLE, rootRank, mpim->world_comm);

	// Unpack into correct body and marker
	int count = 0;
	for (int ib = 0; ib < iBody.size(); ib++) {
		for (int node = 0; node < iBody[ib].markers.size(); node++) {
			iBody[ib].markers[node].epsilon = recvBuffer[count];
			count++;
		}
	}

#endif
}


// ****************************************************************************
