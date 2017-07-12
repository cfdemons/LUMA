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
						volWidth * volDepth * objman->iBody[ib].spacing / objman->iBody[ib]._Owner->dh);
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

// ************************************************************************** //
