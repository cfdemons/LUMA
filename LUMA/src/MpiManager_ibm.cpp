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
		if (my_rank != objman->iBody[ib].owningRank && nMarkers > 0) {

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
		toRank = epsCommMarkerSide[i].toRank;
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
		nSupportSitesRecv[epsCommOwnerSide[i].fromRank].push_back(0);

	// Now loop through and receive
	for (int fromRank = 0; fromRank < num_ranks; fromRank++) {
		if (nSupportSitesRecv[fromRank].size() > 0) {
			MPI_Recv(&nSupportSitesRecv[fromRank].front(), nSupportSitesRecv[fromRank].size(), MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
		}
	}

	// Now unpack
	std::vector<int> idx(num_ranks, 0);
	for (int i = 0; i < epsCommOwnerSide.size(); i++) {
		epsCommOwnerSide[i].nSupportSites = nSupportSitesRecv[epsCommOwnerSide[i].fromRank][idx[epsCommOwnerSide[i].fromRank]];
		idx[epsCommOwnerSide[i].fromRank]++;
	}
}

// ************************************************************************** //
