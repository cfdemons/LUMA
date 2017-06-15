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
				objman->iBody[ib].markers[m].epsilon = epsilon[ib][objman->iBody[ib].markers[m].id];
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
		if (sendBuffer.size() > 0) {
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
	int fromRank, markerIdx;
	for (int i = 0; i < epsCommMarkerSide.size(); i++) {

		// Get ID info
		fromRank = epsCommMarkerSide[i].rankComm;
		bodyID = epsCommMarkerSide[i].bodyID;
		markerIdx = epsCommMarkerSide[i].markerIdx;

		// Put into epsilon
		objman->iBody[bodyID].markers[markerIdx].epsilon = recvBuffer[fromRank][idx[fromRank]];
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
		ib = epsCommMarkerSide[i].bodyID;
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

	// Loop through all bodies
	for (int ib = 0; ib < objman->iBody.size(); ib++) {

		// Get number of markers for this body which exist on this rank
		int nMarkers = objman->iBody[ib].markers.size();

		// Create and size buffer for receiving (also nMarkers to zero as it will not be sending)
		std::vector<int> nMarkersAll;
		if (my_rank == objman->iBody[ib].owningRank) {
			nMarkers = 0;
			nMarkersAll.resize(num_ranks);
		}

		// Gather the number of markers which exist on each rank
		MPI_Gather(&nMarkers, 1, MPI_INT, &nMarkersAll.front(), 1, MPI_INT, objman->iBody[ib].owningRank, world_comm);

		// Send the marker IDs
		if (nMarkers > 0) {

			// Pack into buffer
			std::vector<int> markerIDs;
			for (int m = 0; m < objman->iBody[ib].markers.size(); m++)
				markerIDs.push_back(objman->iBody[ib].markers[m].id);

			// Asynchronous send
			MPI_Request requestStatus;
			MPI_Isend(&markerIDs.front(), markerIDs.size(), MPI_INT, objman->iBody[ib].owningRank, my_rank, world_comm, &requestStatus);
		}


		// Receive marker IDs
		if (my_rank == objman->iBody[ib].owningRank) {

			// Marker ID receive buffer
			std::vector<std::vector<int>> markerIDsAll(num_ranks, std::vector<int>(0));

			// Loop through all ranks to receive from
			for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
				if (nMarkersAll[fromRank] > 0) {

					// Resize buffer and receive
					markerIDsAll[fromRank].resize(nMarkersAll[fromRank]);
					MPI_Recv(&markerIDsAll[fromRank].front(), nMarkersAll[fromRank], MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
				}
			}

			// Now build the eps-comm class on owner side
			for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
				for (int m = 0; m < nMarkersAll[fromRank]; m++) {
					epsCommOwnerSide.emplace_back(fromRank, ib, markerIDsAll[fromRank][m]);
				}
			}
		}

		// Now build eps-comm class on marker side
		for (int m = 0; m < nMarkers; m++)
			epsCommMarkerSide.emplace_back(objman->iBody[ib].owningRank, ib, m);
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
		ib = epsCommMarkerSide[i].bodyID;
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

// ************************************************************************** //
