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
	for (int i = 0; i < supportCommSupportSide.size(); i++) {

		// Get body ID
		ib = objman->bodyIDToIdx[supportCommSupportSide[i].bodyID];

		// Only pack if body belongs to current grid level
		if (objman->iBody[ib]._Owner->level == level) {

			// Get rank to send to
			toRank = supportCommSupportSide[i].rankComm;

			// Get grid sizes
			size_t M_lim = objman->iBody[ib]._Owner->M_lim;
#if (L_DIMS == 3)
			size_t K_lim = objman->iBody[ib]._Owner->K_lim;
#endif

			// Get indices
			idx = supportCommSupportSide[i].supportIdx;

			// Pack velocity into buffer
#if (L_DIMS == 2)
			sendVel[toRank].push_back(objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], eXDirection, M_lim, L_DIMS));
			sendVel[toRank].push_back(objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], eYDirection, M_lim, L_DIMS));
#elif (L_DIMS == 3)
			sendVel[toRank].push_back(objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eXDirection, M_lim, K_lim, L_DIMS));
			sendVel[toRank].push_back(objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eYDirection, M_lim, K_lim, L_DIMS));
			sendVel[toRank].push_back(objman->iBody[ib]._Owner->u(idx[eXDirection], idx[eYDirection], idx[eZDirection], eZDirection, M_lim, K_lim, L_DIMS));
#endif
		}
	}

	// Now send the values
	for (int toRank = 0; toRank < num_ranks; toRank++) {

		// Check if current rank needs to send stuff to rank toRank
		if (sendVel[toRank].size() > 0) {
			MPI_Request send_requests;
			MPI_Isend(&sendVel[toRank].front(), sendVel[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &send_requests);
		}
	}

	// Get receive sizes from each rank
	std::vector<int> bufferSize(num_ranks, 0);
	for (int i = 0; i < supportCommMarkerSide.size(); i++) {
		if (objman->iBody[objman->bodyIDToIdx[supportCommMarkerSide[i].bodyID]]._Owner->level == level)
			bufferSize[supportCommMarkerSide[i].rankComm] += L_DIMS;
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
}

// ************************************************************************* //
/// \brief	Do communication required for epsilon calculation
///
///
void MpiManager::mpi_epsilonCommScatter(std::vector<std::vector<double>> &epsilon) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// First fill the markers which already exist on this rank
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		if (my_rank == objman->iBody[ib].owningRank) {
			for (int m = 0; m < objman->iBody[ib].markers.size(); m++) {
				objman->iBody[ib].markers[m].epsilon = epsilon[objman->iBody[ib].id][objman->iBody[ib].markers[m].id];
			}
		}
	}

	// Declare send buffer
	std::vector<std::vector<double>> sendBuffer(num_ranks, std::vector<double>(0));

	// Pack and send the other epsilon values
	int toRank, bodyID, markerID;
	for (int i = 0; i < epsCommOwnerSide.size(); i++) {

		// Get ID info
		toRank = epsCommOwnerSide[i].rankComm;
		bodyID = epsCommOwnerSide[i].bodyID;
		markerID = epsCommOwnerSide[i].markerID;

		// Insert into buffer
		sendBuffer[toRank].push_back(epsilon[bodyID][markerID]);
	}

	// Loop through and post send if there is data
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (sendBuffer[toRank].size() > 0) {
			MPI_Request requestStatus;
			MPI_Isend(&sendBuffer[toRank].front(), sendBuffer[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &requestStatus);
		}
	}

	// Set up receive buffer sizes
	std::vector<int> bufferSize(num_ranks, 0);

	// Size the receive buffer
	for (int i = 0; i < epsCommMarkerSide.size(); i++)
		bufferSize[epsCommMarkerSide[i].rankComm]++;

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
	int fromRank, ib, markerIdx;
	for (int i = 0; i < epsCommMarkerSide.size(); i++) {

		// Get ID info
		fromRank = epsCommMarkerSide[i].rankComm;
		ib = objman->bodyIDToIdx[epsCommMarkerSide[i].bodyID];
		markerIdx = epsCommMarkerSide[i].markerIdx;

		// Put into epsilon
		objman->iBody[ib].markers[markerIdx].epsilon = recvBuffer[fromRank][idx[fromRank]];
		idx[fromRank]++;
	}
}


// ************************************************************************* //
/// \brief	Do communication required for epsilon calculation
///
///
void MpiManager::mpi_epsilonCommGather(std::vector<std::vector<double>> &recvBuffer) {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Declare send buffer
	std::vector<std::vector<double>> sendBuffer(num_ranks, std::vector<double>(0));

	// Pack the data to send
	int toRank, ib, m;
	for (int i = 0; i < epsCommMarkerSide.size(); i++) {

		// Get ID info
		toRank = epsCommMarkerSide[i].rankComm;
		ib = objman->bodyIDToIdx[epsCommMarkerSide[i].bodyID];
		m = epsCommMarkerSide[i].markerIdx;

		// Pack marker data
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].position[eXDirection]);
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].position[eYDirection]);
#if (L_DIMS == 3)
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].position[eZDirection]);
#endif
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].local_area);
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].dilation);

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
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (sendBuffer[toRank].size() > 0) {
			MPI_Request requestStatus;
			MPI_Isend(&sendBuffer[toRank].front(), sendBuffer[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &requestStatus);
		}
	}

	// Get buffer sizes
	std::vector<int> bufferSize(num_ranks, 0);
	for (int i = 0; i < epsCommOwnerSide.size(); i++)
		bufferSize[epsCommOwnerSide[i].rankComm] += L_DIMS + 2 + epsCommOwnerSide[i].nSupportSites * (L_DIMS + 1);

	// Now create receive buffer
	recvBuffer.resize(num_ranks);

	// Now loop through and receive the buffer
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Only do if there is information to receive
		if (bufferSize[fromRank] > 0) {
			recvBuffer[fromRank].resize(bufferSize[fromRank]);
			MPI_Recv(&recvBuffer[fromRank].front(), recvBuffer[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}
}

// ************************************************************************* //
/// \brief	Owning rank gets how many IBM markers are where for each body
///
///
void MpiManager::mpi_getIBMarkers() {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Build marker side class (only for markers which belong to bodies this rank does not own)
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		if (my_rank != objman->iBody[ib].owningRank) {
			for (int m = 0; m < objman->iBody[ib].markers.size(); m++)
				epsCommMarkerSide.emplace_back(objman->iBody[ib].owningRank, objman->iBody[ib].id, m);
		}
	}

	// Get how many supports site that need to be received from each rank
	std::vector<int> nMarkersToSend(num_ranks, 0);
	for (int i = 0; i < epsCommMarkerSide.size(); i++) {
		nMarkersToSend[epsCommMarkerSide[i].rankComm]++;
	}

	// Do global comm so that each rank knows how many support sites to send and where
	std::vector<int> nMarkersToRecv(num_ranks, 0);
	MPI_Alltoall(&nMarkersToSend.front(), 1, MPI_INT, &nMarkersToRecv.front(), 1, MPI_INT, world_comm);

	// Create send buffer
	std::vector<std::vector<int>> sendBuffer(num_ranks, std::vector<int>(0));

	// Pack body and marker IDs for sending
	int toRank, ib, m;
	for (int i = 0; i < epsCommMarkerSide.size(); i++) {

		// Body and marker ID
		toRank = epsCommMarkerSide[i].rankComm;
		ib = objman->bodyIDToIdx[epsCommMarkerSide[i].bodyID];
		m = epsCommMarkerSide[i].markerIdx;

		// Fill send buffer
		sendBuffer[toRank].push_back(epsCommMarkerSide[i].bodyID);
		sendBuffer[toRank].push_back(objman->iBody[ib].markers[m].id);
	}

	// Loop through and send
	for (toRank = 0; toRank < num_ranks; toRank++) {
		if (sendBuffer[toRank].size() > 0) {
			MPI_Request requestStatus;
			MPI_Isend(&sendBuffer[toRank].front(), sendBuffer[toRank].size(), MPI_INT, toRank, my_rank, world_comm, &requestStatus);
		}
	}

	// Declare receive buffer
	std::vector<std::vector<int>> recvBuffer(num_ranks, std::vector<int>(0));

	// Loop through all receives that are required
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Check if it has stuff to receive from rank i
		if (nMarkersToRecv[fromRank] > 0) {
			recvBuffer[fromRank].resize(2 * nMarkersToRecv[fromRank]);
			MPI_Recv(&recvBuffer[fromRank].front(), recvBuffer[fromRank].size(), MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Now unpack and build comm class
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
		for (int i = 0; i < nMarkersToRecv[fromRank]; i++) {

			// Get ID
			ib = recvBuffer[fromRank][2 * i];
			m = recvBuffer[fromRank][2 * i + 1];

			// Build class
			epsCommOwnerSide.emplace_back(fromRank, ib, m);
		}
	}
}


// ************************************************************************* //
/// \brief	Build communication information required for epsilon calculation.
///
///			The owning rank of each body needs to know how many markers it is
///			receiving and from where. It also needs to know how many support
///			points exist for that marker too. Also, the ranks sending this info
///			need to know this too. This data will be stored in the epsComm class
///
void MpiManager::mpi_buildEpsComms() {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Send buffer which will hold nSupportSites
	std::vector<std::vector<int>> nSupportSitesSend(num_ranks, std::vector<int>(0));

	// Loop through eps-comm class on marker side and pack nSupportSites for sending to owner
	int toRank, ib, m;
	for (int i = 0; i < epsCommMarkerSide.size(); i++) {

		// Get marker info
		toRank = epsCommMarkerSide[i].rankComm;
		ib = objman->bodyIDToIdx[epsCommMarkerSide[i].bodyID];
		m = epsCommMarkerSide[i].markerIdx;

		// Pack into send buffer
		nSupportSitesSend[toRank].push_back(objman->iBody[ib].markers[m].deltaval.size());
	}

	// Now loop through and send to required ranks
	for (int toRank = 0; toRank < num_ranks; toRank++) {
		if (nSupportSitesSend[toRank].size() > 0) {
			MPI_Request requestStatus;
			MPI_Isend(&nSupportSitesSend[toRank].front(), nSupportSitesSend[toRank].size(), MPI_INT, toRank, my_rank, world_comm, &requestStatus);
		}
	}

	// Now create and size the receive buffer
	std::vector<std::vector<int>> nSupportSitesRecv(num_ranks, std::vector<int>(0));
	for (int i = 0; i < epsCommOwnerSide.size(); i++)
		nSupportSitesRecv[epsCommOwnerSide[i].rankComm].push_back(0);

	// Now loop through and receive
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
		if (nSupportSitesRecv[fromRank].size() > 0) {
			MPI_Recv(&nSupportSitesRecv[fromRank].front(), nSupportSitesRecv[fromRank].size(), MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Now unpack
	std::vector<int> idx(num_ranks, 0);
	for (int i = 0; i < epsCommOwnerSide.size(); i++) {
		epsCommOwnerSide[i].nSupportSites = nSupportSitesRecv[epsCommOwnerSide[i].rankComm][idx[epsCommOwnerSide[i].rankComm]];
		idx[epsCommOwnerSide[i].rankComm]++;
	}
}

// ************************************************************************* //
/// \brief	Build communication information required for support communication.
///
///			The owning rank of each marker needs to know where to get support
///			information from.
///
void MpiManager::mpi_buildSupportComms() {

	// Get object manager instance
	ObjectManager *objman = ObjectManager::getInstance();

	// Loop through all bodies, markers and support sites and get places where communication is necessary
	for (int ib = 0; ib < objman->iBody.size(); ib++) {
		for (int m = 0; m < objman->iBody[ib].markers.size(); m++) {
			for (int s = 0; s < objman->iBody[ib].markers[m].deltaval.size(); s++) {

				// If this rank does not own support site then add new element in comm vector
				if (my_rank != objman->iBody[ib].markers[m].support_rank[s])
					supportCommMarkerSide.emplace_back(objman->iBody[ib].id, m, s, objman->iBody[ib].markers[m].support_rank[s]);
			}
		}
	}

	// Get how many supports site that need to be received from each rank
	std::vector<int> nSupportsToRecv(num_ranks, 0);
	for (int i = 0; i < supportCommMarkerSide.size(); i++) {
		nSupportsToRecv[supportCommMarkerSide[i].rankComm]++;
	}

	// Do global comm so that each rank knows how many support sites to send and where
	std::vector<int> nSupportsToSend(num_ranks, 0);
	MPI_Alltoall(&nSupportsToRecv.front(), 1, MPI_INT, &nSupportsToSend.front(), 1, MPI_INT, world_comm);

	// Declare variables for packing positions and body IDs
	int ib, m, s, toRank;
	std::vector<std::vector<double>> supportPositions(num_ranks, std::vector<double>(0));
	std::vector<std::vector<int>> bodyIDs(num_ranks, std::vector<int>(0));

	// Pack the positions of the support sites that need to be received into buffer
	for (int i = 0; i < supportCommMarkerSide.size(); i++) {

		// Get IDs of support site
		ib = objman->bodyIDToIdx[supportCommMarkerSide[i].bodyID];
		m = supportCommMarkerSide[i].markerIdx;
		s = supportCommMarkerSide[i].supportID;
		toRank = supportCommMarkerSide[i].rankComm;

		// Body ID
		bodyIDs[toRank].push_back(supportCommMarkerSide[i].bodyID);

		// Add to send buffer
		supportPositions[toRank].push_back(objman->iBody[ib].markers[m].supp_x[s]);
		supportPositions[toRank].push_back(objman->iBody[ib].markers[m].supp_y[s]);
#if (L_DIMS == 3)
		supportPositions[toRank].push_back(objman->iBody[ib].markers[m].supp_z[s]);
#endif
	}

	// Loop through all sends that are required
	for (int toRank = 0; toRank < num_ranks; toRank++) {

		// Check if it has stuff to send to rank toRank
		if (nSupportsToRecv[toRank] > 0) {
			MPI_Request requestStatusPos, requestStatusIDs;
			MPI_Isend(&supportPositions[toRank].front(), supportPositions[toRank].size(), MPI_DOUBLE, toRank, my_rank, world_comm, &requestStatusPos);
			MPI_Isend(&bodyIDs[toRank].front(), bodyIDs[toRank].size(), MPI_INT, toRank, my_rank, world_comm, &requestStatusIDs);
		}
	}

	// Declare receive buffers
	std::vector<std::vector<double>> supportPositionsRecv(num_ranks, std::vector<double>(0));
	std::vector<std::vector<int>> bodyIDsRecv(num_ranks, std::vector<int>(0));

	// Loop through all receives that are required
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {

		// Resize the buffer
		supportPositionsRecv[fromRank].resize(nSupportsToSend[fromRank] * L_DIMS);
		bodyIDsRecv[fromRank].resize(nSupportsToSend[fromRank]);

		// Check if it has stuff to receive from rank i
		if (nSupportsToSend[fromRank] > 0) {
			MPI_Recv(&supportPositionsRecv[fromRank].front(), supportPositionsRecv[fromRank].size(), MPI_DOUBLE, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
			MPI_Recv(&bodyIDsRecv[fromRank].front(), bodyIDsRecv[fromRank].size(), MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Get enclosing voxel indices
	double x, y, z;
	std::vector<int> nearijk;
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
		for (int i = 0; i < nSupportsToSend[fromRank]; i++) {

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
			supportCommSupportSide.emplace_back(fromRank, bodyIDsRecv[fromRank][i], nearijk);
		}
	}
}

// ************************************************************************** //
