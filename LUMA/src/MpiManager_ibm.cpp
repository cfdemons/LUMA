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

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/ObjectManager.h"



// ************************************************************************* //
/// \brief	Do communication required for force spreading
///
///
void MpiManager::mpi_spreadComm(int level, std::vector<std::vector<double>> &spreadForces) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare values
	int toRank, ib , m, s;
	std::vector<std::vector<double>> sendForce(num_ranks, std::vector<double>(0));

	// Pack and send data first
	for (int i = 0; i < supportCommMarkerSide[level].size(); i++) {

		// Get body index
		ib = objman->bodyIDToIdx[supportCommMarkerSide[level][i].bodyID];

		// Only pack if body belongs to current grid level
		if (objman->iBody[ib]._Owner->level == level) {

			// Get rank to send to and support ID info
			toRank = supportCommMarkerSide[level][i].rankComm;
			m = supportCommMarkerSide[level][i].markerIdx;
			s = supportCommMarkerSide[level][i].supportID;

			// Get volume scaling
			double volWidth = objman->iBody[ib].markers[m].epsilon;
			double volDepth = 1.0;
#if (L_DIMS == 3)
			volDepth = objman->iBody[ib].markers[m].epsilon;
#endif

			// Pack into buffer
			for (int dir = 0; dir < L_DIMS; dir++) {
				sendForce[toRank].push_back(objman->iBody[ib].markers[m].deltaval[s] * objman->iBody[ib].markers[m].force_xyz[dir] *
						volWidth * volDepth * objman->iBody[ib].markers[m].ds);
			}
		}
	}

	// Now send the values
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {

		// Check if current rank needs to send stuff to rank toRank
		if (sendForce[toRank].size() > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&sendForce[toRank].front(), sendForce[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Get receive sizes from each rank
	std::vector<int> bufferSize(num_ranks, 0);
	for (int i = 0; i < supportCommSupportSide[level].size(); i++) {
		if (objman->iBody[objman->bodyIDToIdx[supportCommSupportSide[level][i].bodyID]]._Owner->level == level)
			bufferSize[supportCommSupportSide[level][i].rankComm] += L_DIMS;
	}

	// Now receive the velocity values
	spreadForces.resize(num_ranks, std::vector<double>(0));
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// If there is data to receive from rank fromRank
		if (bufferSize[fromRank] > 0) {
			spreadForces[fromRank].resize(bufferSize[fromRank]);
			MPI_Recv(&spreadForces[fromRank].front(), spreadForces[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}


// ************************************************************************* //
/// \brief	Do communication required for velocity interpolation
///
///
void MpiManager::mpi_interpolateComm(int level, std::vector<std::vector<double>> &recvVel) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare values
	int toRank, ib;
	std::vector<int> idx(3, 0);
	std::vector<std::vector<double>> sendVel(num_ranks, std::vector<double>(0));

	// Pack and send data first
	for (int i = 0; i < supportCommSupportSide[level].size(); i++) {

		// Get body ID
		ib = objman->bodyIDToIdx[supportCommSupportSide[level][i].bodyID];

		// Only pack if body belongs to current grid level
		if (objman->iBody[ib]._Owner->level == level) {

			// Get rank to send to
			toRank = supportCommSupportSide[level][i].rankComm;

			// Get grid sizes
			size_t M_lim = objman->iBody[ib]._Owner->M_lim;
#if (L_DIMS == 3)
			size_t K_lim = objman->iBody[ib]._Owner->K_lim;
#endif

			// Get indices
			idx = supportCommSupportSide[level][i].supportIdx;

			// Pack density and momentum into buffer
#if (L_DIMS == 2)
			sendVel[toRank].push_back(objman->iBody[ib]._Owner->rho(idx[eXDirection], idx[eYDirection], M_lim));

			sendVel[toRank].push_back(objman->iBody[ib]._Owner->rho(idx[eXDirection], idx[eYDirection], M_lim) *
									  objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], eXDirection, M_lim, L_DIMS));

			sendVel[toRank].push_back(objman->iBody[ib]._Owner->rho(idx[eXDirection], idx[eYDirection], M_lim) *
									  objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], eYDirection, M_lim, L_DIMS));
#elif (L_DIMS == 3)
			sendVel[toRank].push_back(objman->iBody[ib]._Owner->rho(idx[eXDirection], idx[eYDirection], idx[eZDirection], M_lim, K_lim));

			sendVel[toRank].push_back(objman->iBody[ib]._Owner->rho(idx[eXDirection], idx[eYDirection], idx[eZDirection], M_lim, K_lim) *
									  objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eXDirection, M_lim, K_lim, L_DIMS));

			sendVel[toRank].push_back(objman->iBody[ib]._Owner->rho(idx[eXDirection], idx[eYDirection], idx[eZDirection], M_lim, K_lim) *
									  objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eYDirection, M_lim, K_lim, L_DIMS));

			sendVel[toRank].push_back(objman->iBody[ib]._Owner->rho(idx[eXDirection], idx[eYDirection], idx[eZDirection], M_lim, K_lim) *
									  objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eZDirection, M_lim, K_lim, L_DIMS));
#endif
		}
	}

	// Now send the values
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {

		// Check if current rank needs to send stuff to rank toRank
		if (sendVel[toRank].size() > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&sendVel[toRank].front(), sendVel[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Get receive sizes from each rank
	std::vector<int> bufferSize(num_ranks, 0);
	for (int i = 0; i < supportCommMarkerSide[level].size(); i++) {
		if (objman->iBody[objman->bodyIDToIdx[supportCommMarkerSide[level][i].bodyID]]._Owner->level == level)
			bufferSize[supportCommMarkerSide[level][i].rankComm] += L_DIMS + 1;
	}

	// Now receive the velocity values
	recvVel.resize(num_ranks, std::vector<double>(0));
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// If there is data to receive from rank fromRank
		if (bufferSize[fromRank] > 0) {
			recvVel[fromRank].resize(bufferSize[fromRank]);
			MPI_Recv(&recvVel[fromRank].front(), recvVel[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}

// ************************************************************************* //
/// \brief	Do communication required for epsilon calculation
///
///
void MpiManager::mpi_epsilonCommScatter(int level) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare send buffer
	std::vector<std::vector<double>> sendBuffer(num_ranks, std::vector<double>(0));

	// Pack and send the other epsilon values
	int toRank, ib, markerID;
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++) {

		// Get ID info
		toRank = markerCommOwnerSide[level][i].rankComm;
		ib = objman->bodyIDToIdx[markerCommOwnerSide[level][i].bodyID];
		markerID = markerCommOwnerSide[level][i].markerID;

		// Insert into buffer
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[markerID].epsilon);
	}

	// Loop through and post send if there is data
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (sendBuffer[toRank].size() > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&sendBuffer[toRank].front(), sendBuffer[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Set up receive buffer sizes
	std::vector<int> bufferSize(num_ranks, 0);

	// Size the receive buffer
	for (int i = 0; i < markerCommMarkerSide[level].size(); i++)
		bufferSize[markerCommMarkerSide[level][i].rankComm]++;

	// Declare receive buffer
	std::vector<std::vector<double>> recvBuffer(num_ranks, std::vector<double>(0));

	// Loop through and receive if there is data
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Only do if there is information to receive
		if (bufferSize[fromRank] > 0) {
			recvBuffer[fromRank].resize(bufferSize[fromRank]);
			MPI_Recv(&recvBuffer[fromRank].front(), recvBuffer[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Create idx vector
	std::vector<int> idx(num_ranks, 0);

	// Now unpack into epsilon values
	int fromRank, markerIdx;
	for (int i = 0; i < markerCommMarkerSide[level].size(); i++) {

		// Get ID info
		fromRank = markerCommMarkerSide[level][i].rankComm;
		ib = objman->bodyIDToIdx[markerCommMarkerSide[level][i].bodyID];
		markerIdx = markerCommMarkerSide[level][i].markerIdx;

		// Put into epsilon
		objman->iBody[ib].markers[markerIdx].epsilon = recvBuffer[fromRank][idx[fromRank]];
		idx[fromRank]++;
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}


// ************************************************************************* //
/// \brief	Do communication required for epsilon calculation
///
///
void MpiManager::mpi_epsilonCommGather(int level) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare send buffer
	std::vector<std::vector<double>> sendBuffer(num_ranks, std::vector<double>(0));

	// Pack the data to send
	int toRank, ib, m;
	for (int i = 0; i < markerCommMarkerSide[level].size(); i++) {

		// Get ID info
		toRank = markerCommMarkerSide[level][i].rankComm;
		ib = objman->bodyIDToIdx[markerCommMarkerSide[level][i].bodyID];
		m = markerCommMarkerSide[level][i].markerIdx;

		// Pack support data
		for (int s = 0; s < objman->iBody[ib].markers[m].deltaval.size(); s++) {
			sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].supp_x[s]);
			sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].supp_y[s]);
#if (L_DIMS == 3)
			sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].supp_z[s]);
#endif
			sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].deltaval[s]);
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
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++)
		bufferSize[markerCommOwnerSide[level][i].rankComm] += markerCommOwnerSide[level][i].nSupportSites * (L_DIMS + 1);

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
	int fromRank, nSupports;
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++) {

		// Get ID info
		fromRank = markerCommOwnerSide[level][i].rankComm;
		ib = objman->bodyIDToIdx[markerCommOwnerSide[level][i].bodyID];
		m = markerCommOwnerSide[level][i].markerID;
		nSupports = markerCommOwnerSide[level][i].nSupportSites;

		// Clear the vectors
		objman->iBody[ib].markers[m].supp_x.resize(nSupports, 0.0);
		objman->iBody[ib].markers[m].supp_y.resize(nSupports, 0.0);
		objman->iBody[ib].markers[m].supp_z.resize(nSupports, 0.0);
		objman->iBody[ib].markers[m].deltaval.resize(nSupports, 0.0);

		// Unpack into iBody
		for (int s = 0; s < nSupports; s++) {
			objman->iBody[ib].markers[m].supp_x[s] = recvBuffer[fromRank][idx[fromRank]];
			objman->iBody[ib].markers[m].supp_y[s] = recvBuffer[fromRank][idx[fromRank]+1];
#if (L_DIMS == 3)
			objman->iBody[ib].markers[m].supp_z[s] = recvBuffer[fromRank][idx[fromRank]+2];
#endif
			idx[fromRank] += L_DIMS;

			// Get delta values
			objman->iBody[ib].markers[m].deltaval[s] = recvBuffer[fromRank][idx[fromRank]];
			idx[fromRank]++;
		}
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}


// ************************************************************************* //
/// \brief	Do communication required for universal epsilon calculation
///
///
void MpiManager::mpi_uniEpsilonCommGather(int level, int rootRank, IBBody &iBodyTmp) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Initialise values
	int nMarkersOnThisRank = 0;
	std::vector<int> deltasPerMarkerOnThisRank;

	// Get how many markers and supports exist on each rank
	for (int ib = 0; ib < objman->iBody.size(); ib++) {

		// Number of markers on this rank
		nMarkersOnThisRank += objman->iBody[ib].validMarkers.size();

		// Loop through markers and get number of deltas
		for (auto m : objman->iBody[ib].validMarkers)
			deltasPerMarkerOnThisRank.push_back(objman->iBody[ib].markers[m].deltaval.size());
	}

	// Set receive buffer for how many markers on each rank
	std::vector<int> nMarkersOnEachRank;
	if (my_rank == rootRank)
		nMarkersOnEachRank.resize(num_ranks);

	// Perform gather so root rank knows how many markers on each rank
	MPI_Gather(&nMarkersOnThisRank, 1, MPI_INT, &nMarkersOnEachRank.front(), 1, MPI_INT, rootRank, world_comm);

	// Get number of markers in whole system
	int nSystemMarkers = std::accumulate(nMarkersOnEachRank.begin(), nMarkersOnEachRank.end(), 0);

	// Create result and displacement arrays
	std::vector<int> deltasPerMarkerOnEachRank;
	std::vector<int> markerDisps;

	// Root rank resize result and displacement arrays
	if (my_rank == rootRank) {

		// Resize result vector
		deltasPerMarkerOnEachRank.resize(nSystemMarkers);
		markerDisps.resize(num_ranks);

		// Calculate block displacements
		for (int i = 1; i < num_ranks; i++)
			markerDisps[i] = std::accumulate(nMarkersOnEachRank.begin(), nMarkersOnEachRank.begin()+i, 0);
	}

	// Gather so root rank knows how many support sites each markers has
	MPI_Gatherv(&deltasPerMarkerOnThisRank.front(), nMarkersOnThisRank, MPI_INT, &deltasPerMarkerOnEachRank.front(), &nMarkersOnEachRank.front(), &markerDisps.front(), MPI_INT, rootRank, world_comm);

	// Declare receive buffer vectors
	std::vector<int> recvBufferSizes;
	std::vector<int> recvBufferDisps;
	std::vector<double> recvBuffer;

	// Now create receive buffer, buffer size and buffer displacements
	if (my_rank == rootRank) {

		// Resize for number of ranks
		recvBufferSizes.resize(num_ranks, 0);
		recvBufferDisps.resize(num_ranks, 0);

		// Counter
		int count = 0;

		// Loop through all ranks and sort out receive buffer
		for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

			// Loop through all markers which exist on this rank
			for (int m = 0; m < nMarkersOnEachRank[fromRank]; m++) {

				// Add contribution to buffer size
				recvBufferSizes[fromRank] += (L_DIMS + 3) + deltasPerMarkerOnEachRank[count] * (L_DIMS + 1);

				// Increment counter
				count++;
			}
		}

		// Calculate block displacements
		for (int i = 1; i < num_ranks; i++)
			recvBufferDisps[i] = std::accumulate(recvBufferSizes.begin(), recvBufferSizes.begin()+i, 0);

		// Resize receive buffer
		recvBuffer.resize(std::accumulate(recvBufferSizes.begin(), recvBufferSizes.end(), 0), 0.0);
	}

	// Declare send buffer vector
	std::vector<double> sendBuffer;

	// Now pack data into send buffer
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		for (auto m : objman->iBody[ib].validMarkers) {

			// Pack position first
			sendBuffer.push_back(objman->iBody[ib].markers[m].position[eXDirection]);
			sendBuffer.push_back(objman->iBody[ib].markers[m].position[eYDirection]);
#if (L_DIMS == 3)
			sendBuffer.push_back(objman->iBody[ib].markers[m].position[eZDirection]);
#endif
			sendBuffer.push_back(objman->iBody[ib].markers[m].local_area);
			sendBuffer.push_back(objman->iBody[ib].markers[m].dilation);
			sendBuffer.push_back(objman->iBody[ib].markers[m].ds);

			// Now pack the support site data
			for (int s = 0; s < objman->iBody[ib].markers[m].deltaval.size(); s++) {
				sendBuffer.push_back(objman->iBody[ib].markers[m].supp_x[s]);
				sendBuffer.push_back(objman->iBody[ib].markers[m].supp_y[s]);
#if (L_DIMS == 3)
				sendBuffer.push_back(objman->iBody[ib].markers[m].supp_z[s]);
#endif
				sendBuffer.push_back(objman->iBody[ib].markers[m].deltaval[s]);
			}
		}
	}

	// Now send all data to root
	MPI_Gatherv(&sendBuffer.front(), sendBuffer.size(), MPI_DOUBLE, &recvBuffer.front(), &recvBufferSizes.front(), &recvBufferDisps.front(), MPI_DOUBLE, rootRank, world_comm);

	// Set global values
	iBodyTmp.owningRank = rootRank;
	iBodyTmp.level = level;
	iBodyTmp._Owner = objman->iBody[0]._Owner;

	// Unpack receive buffer
	if (my_rank == rootRank) {

		// Resize number of markers
		iBodyTmp.markers.resize(nSystemMarkers);

		// Counter
		int count = 0;

		// Marker ID
		int m = 0;

		// Now loop through all markers
		for (int rank = 0; rank < num_ranks; rank++) {
			for (int marker = 0; marker < nMarkersOnEachRank[rank]; marker++) {

				// Get owning rank of marker (need this for the scatter later
				iBodyTmp.markers[m].owningRank = rank;

				// Pack position
				for (int d = 0; d < L_DIMS; d++)
					iBodyTmp.markers[m].position[d] = recvBuffer[count+d];

				// Increment
				count += L_DIMS;

				// Pack area and dilation
				iBodyTmp.markers[m].local_area = recvBuffer[count]; count++;
				iBodyTmp.markers[m].dilation = recvBuffer[count]; count++;
				iBodyTmp.markers[m].ds = recvBuffer[count]; count++;

				// Now pack support data
				for (int s = 0; s < deltasPerMarkerOnEachRank[m]; s++) {
					iBodyTmp.markers[m].supp_x.push_back(recvBuffer[count]); count++;
					iBodyTmp.markers[m].supp_y.push_back(recvBuffer[count]); count++;
#if (L_DIMS == 3)
					iBodyTmp.markers[m].supp_z.push_back(recvBuffer[count]); count++;
#endif
					iBodyTmp.markers[m].deltaval.push_back(recvBuffer[count]); count++;
				}

				// Increment m
				m++;
			}
		}
	}
}


// ************************************************************************* //
/// \brief	Do communication required for universal epsilon calculation
///
///
void MpiManager::mpi_uniEpsilonCommScatter(int rootRank, IBBody &iBodyTmp) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare send buffer vectors
	std::vector<int> sendBufferSizes;
	std::vector<int> sendBufferDisps;
	std::vector<double> sendBuffer;

	// Root rank pack data into send buffer and sort sizes
	if (my_rank == rootRank) {

		// Set sizes
		sendBufferSizes.resize(num_ranks, 0);
		sendBufferDisps.resize(num_ranks, 0);

		// Loop through pack epsilon and insert buffer sizes for each rank
		for (int m = 0; m < iBodyTmp.markers.size(); m++) {
			sendBuffer.push_back(iBodyTmp.markers[m].epsilon);
			sendBufferSizes[iBodyTmp.markers[m].owningRank]++;
		}

		// Calculate block displacements
		for (int i = 1; i < num_ranks; i++)
			sendBufferDisps[i] = std::accumulate(sendBufferSizes.begin(), sendBufferSizes.begin()+i, 0);
	}

	// Get how many markers and supports exist on each rank
	int nMarkersOnThisRank = 0;
	for (int ib = 0; ib < objman->iBody.size(); ib++)
		nMarkersOnThisRank += objman->iBody[ib].validMarkers.size();

	// Create receive buffer
	std::vector<double> recvBuffer(nMarkersOnThisRank, 0.0);

	// Scatter the epsilon values to the correct rank
	MPI_Scatterv(&sendBuffer.front(), &sendBufferSizes.front(), &sendBufferDisps.front(), MPI_DOUBLE, &recvBuffer.front(), nMarkersOnThisRank, MPI_DOUBLE, rootRank, world_comm);

	// Now unpack
	int count = 0;
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		for (auto m : objman->iBody[ib].validMarkers) {
			objman->iBody[ib].markers[m].epsilon = recvBuffer[count];
			count++;
		}
	}

	// Update owning rank's epsilon values
	std::vector<std::vector<double>> sendEpsBuffer(num_ranks, std::vector<double>(0));

	// Pack and send the other epsilon values
	int toRank, ib, markerIdx;
	for (int lev = 0; lev < markerCommMarkerSide.size(); lev++) {
		for (int i = 0; i < markerCommMarkerSide[lev].size(); i++) {

			// Get ID info
			toRank = markerCommMarkerSide[lev][i].rankComm;
			ib = objman->bodyIDToIdx[markerCommMarkerSide[lev][i].bodyID];
			markerIdx = markerCommMarkerSide[lev][i].markerIdx;

			// Insert into buffer
			sendEpsBuffer[toRank].push_back(objman->iBody[ib].markers[markerIdx].epsilon);
		}
	}

	// Loop through and post send if there is data
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (sendEpsBuffer[toRank].size() > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&sendEpsBuffer[toRank].front(), sendEpsBuffer[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Set up receive buffer sizes
	std::vector<int> bufferSize(num_ranks, 0);

	// Size the receive buffer
	for (int lev = 0; lev < markerCommOwnerSide.size(); lev++) {
		for (int i = 0; i < markerCommOwnerSide[lev].size(); i++) {
			bufferSize[markerCommOwnerSide[lev][i].rankComm]++;
		}
	}

	// Declare receive buffer
	std::vector<std::vector<double>> recvEpsBuffer(num_ranks, std::vector<double>(0));

	// Loop through and receive if there is data
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Only do if there is information to receive
		if (bufferSize[fromRank] > 0) {
			recvEpsBuffer[fromRank].resize(bufferSize[fromRank]);
			MPI_Recv(&recvEpsBuffer[fromRank].front(), recvEpsBuffer[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Create idx vector
	std::vector<int> idx(num_ranks, 0);

	// Now unpack into epsilon values
	int fromRank, markerID;
	for (int lev = 0; lev < markerCommOwnerSide.size(); lev++) {
		for (int i = 0; i < markerCommOwnerSide[lev].size(); i++) {

			// Get ID info
			fromRank = markerCommOwnerSide[lev][i].rankComm;
			ib = objman->bodyIDToIdx[markerCommOwnerSide[lev][i].bodyID];
			markerID = markerCommOwnerSide[lev][i].markerID;

			// Put into epsilon
			objman->iBody[ib].markers[markerID].epsilon = recvEpsBuffer[fromRank][idx[fromRank]];
			idx[fromRank]++;
		}
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}

// ************************************************************************* //
/// \brief	Build helper classes for bodyOwner-marker communication
///
///
void MpiManager::mpi_buildMarkerComms(int level) {

	// Clear the vector
	markerCommOwnerSide[level].clear();
	markerCommMarkerSide[level].clear();

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Build owner side for marker comms
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		if (my_rank == objman->iBody[ib].owningRank && level == objman->iBody[ib]._Owner->level) {
			for (int m = 0; m < objman->iBody[ib].markers.size(); m++) {
				if (my_rank != objman->iBody[ib].markers[m].owningRank)
					markerCommOwnerSide[level].emplace_back(objman->iBody[ib].markers[m].owningRank, objman->iBody[ib].id, objman->iBody[ib].markers[m].id);
			}
		}
	}

	// Build marker side for marker comms
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		if (my_rank != objman->iBody[ib].owningRank && level == objman->iBody[ib]._Owner->level) {
			for (auto m : objman->iBody[ib].validMarkers)
				markerCommMarkerSide[level].emplace_back(objman->iBody[ib].owningRank, objman->iBody[ib].id, m);
		}
	}

	// Send number of support sites
	std::vector<std::vector<int>> sendBuffer(num_ranks, std::vector<int>(0));

	// Pack data
	int toRank, ib, m;
	for (int i = 0; i < markerCommMarkerSide[level].size(); i++) {

		// Get marker info
		toRank = markerCommMarkerSide[level][i].rankComm;
		ib = objman->bodyIDToIdx[markerCommMarkerSide[level][i].bodyID];
		m = markerCommMarkerSide[level][i].markerIdx;

		// Pack into send buffer
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].deltaval.size());
	}

	// Now loop through and send to required ranks
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (sendBuffer[toRank].size() > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&sendBuffer[toRank].front(), sendBuffer[toRank].size(), MPI_INT, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Now create and size the receive buffer
	std::vector<std::vector<int>> recvBuffer(num_ranks, std::vector<int>(0));
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++)
		recvBuffer[markerCommOwnerSide[level][i].rankComm].push_back(0);

	// Now loop through and receive
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
		if (recvBuffer[fromRank].size() > 0) {
			MPI_Recv(&recvBuffer[fromRank].front(), recvBuffer[fromRank].size(), MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Now unpack
	std::vector<int> idx(num_ranks, 0);
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++) {
		markerCommOwnerSide[level][i].nSupportSites = recvBuffer[markerCommOwnerSide[level][i].rankComm][idx[markerCommOwnerSide[level][i].rankComm]];
		idx[markerCommOwnerSide[level][i].rankComm]++;
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}

// ************************************************************************* //
/// \brief	Build communication information required for support communication.
///
///			The owning rank of each marker needs to know where to get support
///			information from.
///
void MpiManager::mpi_buildSupportComms(int level) {

	// Clear the vector
	supportCommMarkerSide[level].clear();
	supportCommSupportSide[level].clear();

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Loop through all bodies, markers and support sites and get places where communication is necessary
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		if (objman->iBody[ib]._Owner->level == level) {
			for (auto m : objman->iBody[ib].validMarkers) {
				for (int s = 0; s < objman->iBody[ib].markers[m].deltaval.size(); s++) {

					// If this rank does not own support site then add new element in comm vector
					if (my_rank != objman->iBody[ib].markers[m].support_rank[s])
						supportCommMarkerSide[level].emplace_back(objman->iBody[ib].id, m, s, objman->iBody[ib].markers[m].support_rank[s]);
				}
			}
		}
	}

	// Get how many supports site that need to be received from each rank
	std::vector<int> nSupportToRecv(num_ranks, 0);
	for (int i = 0; i < supportCommMarkerSide[level].size(); i++) {
		nSupportToRecv[supportCommMarkerSide[level][i].rankComm]++;
	}

	// Do global comm so that each rank knows how many support sites to send and where
	std::vector<int> nSupportToSend(num_ranks, 0);
	MPI_Alltoall(&nSupportToRecv.front(), 1, MPI_INT, &nSupportToSend.front(), 1, MPI_INT, world_comm);

	// Declare variables for packing positions and body IDs
	int ib, m, s, toRank;
	std::vector<std::vector<double>> supportPositions(num_ranks, std::vector<double>(0));
	std::vector<std::vector<int>> bodyIDs(num_ranks, std::vector<int>(0));

	// Pack the positions of the support sites that need to be received into buffer
	for (int i = 0; i < supportCommMarkerSide[level].size(); i++) {

		// Get IDs of support site
		ib = objman->bodyIDToIdx[supportCommMarkerSide[level][i].bodyID];
		m = supportCommMarkerSide[level][i].markerIdx;
		s = supportCommMarkerSide[level][i].supportID;
		toRank = supportCommMarkerSide[level][i].rankComm;

		// Body ID
		bodyIDs[toRank].push_back(supportCommMarkerSide[level][i].bodyID);

		// Add to send buffer
		supportPositions[toRank].push_back(objman->iBody[ib].markers[m].supp_x[s]);
		supportPositions[toRank].push_back(objman->iBody[ib].markers[m].supp_y[s]);
#if (L_DIMS == 3)
		supportPositions[toRank].push_back(objman->iBody[ib].markers[m].supp_z[s]);
#endif
	}

	// Loop through all sends that are required
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {

		// Check if it has stuff to send to rank toRank
		if (nSupportToRecv[toRank] > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&supportPositions[toRank].front(), supportPositions[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &sendRequests.back());
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&bodyIDs[toRank].front(), bodyIDs[toRank].size(), MPI_INT, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Declare receive buffers
	std::vector<std::vector<double>> supportPositionsRecv(num_ranks, std::vector<double>(0));
	std::vector<std::vector<int>> bodyIDsRecv(num_ranks, std::vector<int>(0));

	// Loop through all receives that are required
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Resize the buffer
		supportPositionsRecv[fromRank].resize(nSupportToSend[fromRank] * L_DIMS);
		bodyIDsRecv[fromRank].resize(nSupportToSend[fromRank]);

		// Check if it has stuff to receive from rank i
		if (nSupportToSend[fromRank] > 0) {
			MPI_Recv(&supportPositionsRecv[fromRank].front(), supportPositionsRecv[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
			MPI_Recv(&bodyIDsRecv[fromRank].front(), bodyIDsRecv[fromRank].size(), MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Get enclosing voxel indices
	double x, y, z;
	std::vector<int> nearijk;
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
		for (int i = 0; i < nSupportToSend[fromRank]; i++) {

			// Unpack positions and body ID
			ib = objman->bodyIDToIdx[bodyIDsRecv[fromRank][i]];
			x = supportPositionsRecv[fromRank][i * L_DIMS];
			y = supportPositionsRecv[fromRank][i * L_DIMS + 1];
			z = 0.0;
#if (L_DIMS == 3)
			z = supportPositionsRecv[fromRank][i * L_DIMS + 2];
#endif

			// Get indices of enclosing voxel
			GridUtils::getEnclosingVoxel(x, y, z, objman->iBody[ib]._Owner, &nearijk);

			// Place in sending vector
			supportCommSupportSide[level].emplace_back(fromRank, bodyIDsRecv[fromRank][i], nearijk);
		}
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}


// ************************************************************************* //
/// \brief	Spread ds values from owning rank to other ranks
///
///
void MpiManager::mpi_dsCommScatter(int level) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare send buffer
	std::vector<std::vector<double>> sendBuffer(num_ranks, std::vector<double>(0));

	// Pack and send the other epsilon values
	int toRank, ib, markerID;
	for (int i = 0; i < markerCommOwnerSide[level].size(); i++) {

		// Get ID info
		toRank = markerCommOwnerSide[level][i].rankComm;
		ib = objman->bodyIDToIdx[markerCommOwnerSide[level][i].bodyID];
		markerID = markerCommOwnerSide[level][i].markerID;

		// Insert into buffer
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[markerID].ds);
	}

	// Loop through and post send if there is data
	std::vector<MPI_Request> sendRequests;
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (sendBuffer[toRank].size() > 0) {
			sendRequests.push_back(MPI_REQUEST_NULL);
			MPI_Isend(&sendBuffer[toRank].front(), sendBuffer[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &sendRequests.back());
		}
	}

	// Set up receive buffer sizes
	std::vector<int> bufferSize(num_ranks, 0);

	// Size the receive buffer
	for (int i = 0; i < markerCommMarkerSide[level].size(); i++)
		bufferSize[markerCommMarkerSide[level][i].rankComm]++;

	// Declare receive buffer
	std::vector<std::vector<double>> recvBuffer(num_ranks, std::vector<double>(0));

	// Loop through and receive if there is data
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Only do if there is information to receive
		if (bufferSize[fromRank] > 0) {
			recvBuffer[fromRank].resize(bufferSize[fromRank]);
			MPI_Recv(&recvBuffer[fromRank].front(), recvBuffer[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Create idx vector
	std::vector<int> idx(num_ranks, 0);

	// Now unpack into epsilon values
	int fromRank, markerIdx;
	for (int i = 0; i < markerCommMarkerSide[level].size(); i++) {

		// Get ID info
		fromRank = markerCommMarkerSide[level][i].rankComm;
		ib = objman->bodyIDToIdx[markerCommMarkerSide[level][i].bodyID];
		markerIdx = markerCommMarkerSide[level][i].markerIdx;

		// Put into epsilon
		objman->iBody[ib].markers[markerIdx].ds = recvBuffer[fromRank][idx[fromRank]];
		idx[fromRank]++;
	}

	// If sending any messages then wait for request status
	MPI_Waitall(sendRequests.size(), &sendRequests.front(), MPI_STATUS_IGNORE);
}

// ************************************************************************** //
