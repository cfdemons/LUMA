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

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/ObjectManager.h"
#include "../inc/MpiManager.h"
#include "../inc/GridUtils.h"


// *****************************************************************************
///	\brief Performs MPI communication depending on enumeration type case.
///	\param	type		type of communication.
///	\param	iBody		the current IB body.
///	\param	numMarkers	where the buffer is being unpacked to
void ObjectManager::ibm_getEpsInfo(std::vector<IBBody> &iBody, std::vector<IBBody> &iBodyEps) {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	// Values needed for send
	MPI_Request send_requests[iBody.size()] = {MPI_REQUEST_NULL};
	MPI_Request barrier_requests[iBody.size()] = {MPI_REQUEST_NULL};

	// Each rank needs a vector with the body IDs it is responsible for
	std::vector<int> ownedBodies;

	// Vector of where each rank has sent messages
	std::vector<int> sentToRank(mpim->num_ranks, 0);

	// Create a vector sized for number of total iBodies as can't be reused
	std::vector<std::vector<double>> bufferSend;
	bufferSend.resize(iBody.size());

	// Loop through each iBody
	for (int ib = 0; ib < iBody.size(); ib++) {

		// Check if it is responsible for this body
		if (mpim->my_rank == iBody[ib].owningRank)
			ownedBodies.push_back(ib);

		// Only send if there is marker data for this body
		if (iBody[ib].markers.size() > 0) {

			// Pack into send buffer (see implementation for details on data storage)
			ibm_mpi_pack(&iBody[ib], bufferSend[ib]);

			// Send the data with a non-blocking send
			MPI_Isend(&bufferSend[ib].front(), bufferSend[ib].size(), MPI_DOUBLE, iBody[ib].owningRank, ib, mpim->world_comm, &send_requests[ib]);
			std::cout << "Rank " << mpim->my_rank << " has sent a message to rank " << iBody[ib].owningRank << std::endl;

			// Add to the sentToRank vector
			sentToRank[iBody[ib].owningRank] += 1;
		}
	}


	// Values needed for send
	int bufferSize;
	std::vector<double> bufferRecvTemp;
	std::vector<std::vector<double>> bufferRecv;
	bufferRecv.resize(ownedBodies.size());
	MPI_Status probeStatus;
	int probeFlag;
	int receiveCount = 0;
	std::vector<int>::iterator elIdx;

	// Resize the iBodyEps vector
	iBodyEps.resize(ownedBodies.size());

	// Do an AllReduce so that each process knows how many sends it is expecting to receive
	std::vector<int> expectingReceives(mpim->num_ranks, 0);
	MPI_Allreduce(&sentToRank.front(), &expectingReceives.front(), mpim->num_ranks, MPI_INT, MPI_SUM, mpim->world_comm);


	// If there is a message waiting to be received
	while (receiveCount < expectingReceives[mpim->my_rank]) {

		// Now check for sends this rank needs to receive
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpim->world_comm, &probeFlag, &probeStatus);

		// If a message is there
		if (probeFlag == 1 && std::find(ownedBodies.begin(), ownedBodies.end(), probeStatus.MPI_TAG) != ownedBodies.end()) {

			// Get the element index
			elIdx = std::find(ownedBodies.begin(), ownedBodies.end(), probeStatus.MPI_TAG);

			// Get the size of the message
			MPI_Get_count(&probeStatus, MPI_DOUBLE, &bufferSize);

			// Resize the receving buffer
			bufferRecvTemp.resize(bufferSize);

			// Receive the message
			MPI_Recv(&bufferRecvTemp.front(), bufferSize, MPI_DOUBLE, probeStatus.MPI_SOURCE, probeStatus.MPI_TAG, mpim->world_comm, MPI_STATUS_IGNORE);
			std::cout << "Rank " << mpim->my_rank << " has received a message from rank " << probeStatus.MPI_SOURCE << " addressed to rank " << probeStatus.MPI_TAG << std::endl;

			// Put it into the main buffer
			bufferRecv[*elIdx].insert(bufferRecv[*elIdx].end(), bufferRecvTemp.begin(), bufferRecvTemp.end());

			// Increment the receive counter
			receiveCount++;
		}
	}

	// Unpack into the proper iBody format
	ibm_mpi_unpack(bufferRecv, iBodyEps);


	// Write out the message that has been received
	if (mpim->my_rank == 0) {
		std::cout << "Rank " << mpim->my_rank << std::endl;
		for (int i = 0; i < bufferRecv.size(); i++) {
			for (int j = 0; j < bufferRecv[i].size(); j++) {
				std::cout << bufferRecv[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}


	MPI_Barrier(mpim->world_comm);
	exit(0);
}



///	\brief Packs the send buffer into bufferSend.
///	\param	type			type of communication.
///	\param	iBody			pointer to the current iBody
///	\param	buffer			the sender buffer that is being packed.
void ObjectManager::ibm_mpi_pack(IBBody *iBody, std::vector<double> &buffer) {


	// The data storage for this buffer is x, y, (z), area, dilation for each marker
	// Loop through each marker and insert values
	for (int m = 0; m < iBody->markers.size(); m++) {
		buffer.push_back(iBody->markers[m].position[eXDirection]);
		buffer.push_back(iBody->markers[m].position[eYDirection]);
#if (L_DIMS == 3)
		buffer.push_back(iBody->markers[m].position[eZDirection]);
#endif
		buffer.push_back(iBody->markers[m].local_area);
		buffer.push_back(iBody->markers[m].dilation);

		// Now add the support point positions and delta values
		for (int s = 0; s < iBody->markers[m].supp_x.size(); s++) {
			buffer.push_back(iBody->markers[m].supp_x[s]);
			buffer.push_back(iBody->markers[m].supp_y[s]);
#if (L_DIMS == 3)
			buffer.push_back(iBody->markers[m].supp_x[s]);
#endif
			buffer.push_back(iBody->markers[m].deltaval[s]);
		}
	}
}


///	\brief Unpacks the receive buffer into iBodyEps.
///	\param	bufferRecv			receive buffer which contains all the info for a particular iBody
///	\param	iBodyEps			vector of the final IBBody class which will contain all info for epsilon calculation.
void ObjectManager::ibm_mpi_unpack(std::vector<std::vector<double>> &bufferRecv, std::vector<IBBody> &iBodyEps) {



}
// ****************************************************************************
