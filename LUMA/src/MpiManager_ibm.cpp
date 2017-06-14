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

		// Create and size buffer for receiving
		std::vector<int> nMarkersAll;
		if (my_rank == objman->iBody[ib].owningRank)
			nMarkersAll.resize(num_ranks);

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
				if (my_rank != fromRank && nMarkersAll[fromRank] > 0) {

					// Resize buffer and receive
					markerIDsAll[fromRank].resize(nMarkersAll[fromRank]);
					MPI_Recv(&markerIDsAll[fromRank].front(), nMarkersAll[fromRank], MPI_INT, fromRank, fromRank, world_comm, MPI_STATUS_IGNORE);
				}
			}
		}
	}

	MPI_Barrier(world_comm);
	exit(0);

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

	//

}

// ************************************************************************** //
