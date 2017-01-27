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
	int sendCount = 0;

	// Each rank needs a vector with the body IDs it is responsible for
	std::vector<int> ownedBodies;

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
			MPI_Isend(&bufferSend[ib].front(), bufferSend[ib].size(), MPI_DOUBLE, iBody[ib].owningRank, iBody[ib].owningRank, mpim->world_comm, &send_requests[ib]);
			sendCount++;
			std::cout << "Rank " << mpim->my_rank << " has sent a message to rank " << iBody[ib].owningRank << std::endl;
		}
	}


	// Values needed for send
	int bufferSize;
	std::vector<double> bufferRecvTemp;
	std::vector<std::vector<double>> bufferRecv;
	bufferRecv.resize(iBody.size());
	MPI_Status probeStatus;
	int probeFlag;
	int testFlag;

	// Loop through each iBody
	for (int ib = 0; ib < ownedBodies.size(); ib++) {

		// Wait until all messages for a particular body have been sent
		MPI_Test(&send_requests[ownedBodies[ib]], &testFlag, MPI_STATUS_IGNORE);
		std::cout << "Testflag = " << testFlag << std::endl;

		// While there is a send for this rank waiting to be received
		while (testFlag == 0) {

			// Now check for sends this rank needs to receive
			MPI_Iprobe(MPI_ANY_SOURCE, mpim->my_rank, mpim->world_comm, &probeFlag, &probeStatus);

			std::cout << "Rank " << mpim->my_rank << " has found a message from " << probeStatus.MPI_SOURCE << " addressed to " << probeStatus.MPI_TAG << std::endl;

			// Get buffer size for receiving and resize the vector
			MPI_Get_count(&probeStatus, MPI_DOUBLE, &bufferSize);
			bufferRecvTemp.resize(bufferSize);

			// Receive the message
			MPI_Recv(&bufferRecvTemp.front(), bufferSize, MPI_DOUBLE, probeStatus.MPI_SOURCE, probeStatus.MPI_TAG, mpim->world_comm, MPI_STATUS_IGNORE);

			// Put it into the main buffer
			bufferRecv[ib].insert(bufferRecv[ib].end(), bufferRecvTemp.begin(), bufferRecvTemp.end());

			// Now check for sends this rank needs to receive
			MPI_Test(&send_requests[ownedBodies[ib]], &testFlag, MPI_STATUS_IGNORE);
			std::cout << "Testflag = " << testFlag << std::endl;

		}
	}


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
	}
}


// *****************************************************************************
///	\brief Performs MPI communication depending on enumeration type case.
///	\param	type		type of communication.
///	\param	iBody		the current IB body.
///	\param	numMarkers	where the buffer is being unpacked to
void ObjectManager::ibm_getMarkerPositions(eIBInfoType type, IBBody *iBody, std::vector<std::vector<double>> &markerPos, std::vector<int> &rankID, std::vector<int> &numMarkers) {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	switch (type)
	{
	case eBodyMarkerPositions:
		{
			// Data type to be sent
			MPI_Datatype mpi_type = MPI_DOUBLE;

			// Get how many markers from each rank the owning rank should ask for and pack into buffer
			std::vector<double> bufferSend;
			if (iBody->markers.size() > 0)
				ibm_mpi_pack(type, iBody, bufferSend);

			// Instantiate the receive buffer
			std::vector<std::vector<double>> bufferRecv;

			// Only resize buffer if this is the owning rank
			if (mpim->my_rank == iBody->owningRank) {
				bufferRecv.resize(rankID.size());
				for (int i = 0; i < bufferRecv.size(); i++) {
					bufferRecv[i].resize(numMarkers[i] * L_DIMS);
				}
			}

			// Send the message
			if (iBody->markers.size() > 0 && mpim->my_rank != iBody->owningRank)
				MPI_Send(&bufferSend.front(), bufferSend.size(), mpi_type, iBody->owningRank, mpim->my_rank, mpim->world_comm);

			// Receive the message
			if (mpim->my_rank == iBody->owningRank) {
				for (int i = 0; i < rankID.size(); i++) {
					MPI_Recv(&bufferRecv[i].front(), bufferRecv[i].size(), mpi_type, rankID[i], rankID[i], mpim->world_comm, MPI_STATUS_IGNORE);
				}
			}

			// Unpack back into required vectors class
			if (mpim->my_rank == iBody->owningRank)
				ibm_mpi_unpack(type, bufferRecv, rankID, numMarkers, markerPos);

			break;
		}
	}
}

///	\brief Packs the send buffer into bufferSend.
///	\param	type			type of communication.
///	\param	iBody			pointer to the current iBody
///	\param	buffer			the sender buffer that is being packed.
void ObjectManager::ibm_mpi_pack(eIBInfoType type, IBBody *iBody, std::vector<int> &buffer) {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	switch (type)
	{

	case eBodyNumMarkers:
		{

			// Fill the vector with required data
			buffer.push_back(mpim->my_rank);
			buffer.push_back(iBody->markers.size());

			break;
		}
	}
}


///	\brief Packs the send buffer into bufferSend.
///	\param	type			type of communication.
///	\param	iBody			pointer to the current iBody
///	\param	buffer			the sender buffer that is being packed.
void ObjectManager::ibm_mpi_pack(eIBInfoType type, IBBody *iBody, std::vector<double> &buffer) {

	switch (type)
	{

	case eBodyMarkerPositions:
		{

			// Resize the buffer first
			buffer.resize(iBody->markers.size() * L_DIMS);

			// Loop through all markers to be sent and add them to the buffer
			for (int m = 0; m < iBody->markers.size(); m++) {
				buffer[m*L_DIMS] = iBody->markers[m].position[eXDirection];
				buffer[m*L_DIMS+1] = iBody->markers[m].position[eYDirection];
#if (L_DIMS == 3)
				buffer[m*L_DIMS+2] = iBody->markers[m].position[eZDirection];
#endif

			}

			break;
		}
	}
}


///	\brief Unpacks the receive buffer back into IBInfo class.
///	\param	type			type of communication.
///	\param	bufferRecv		the receiver buffer that is being unpacked.
///	\param	numMarkers		where the buffer is being unpacked to
void ObjectManager::ibm_mpi_unpack(eIBInfoType type, std::vector<std::vector<double>> &bufferRecv, std::vector<int> &rankID, std::vector<int> &numMarkers, std::vector<std::vector<double>> &markerPos) {



	switch (type)
	{
	case eBodyMarkerPositions:
		{

			// Get the total number of markers
			int nMarkers = std::accumulate(numMarkers.begin(), numMarkers.end(), 0);

			// Resize the marker position vector
			markerPos.resize(nMarkers);

			// Loop through all markers and feed in x, y (and z) positions
			int count = 0;
			for (int i = 0; i < bufferRecv.size(); i++) {
				for (int j = 0; j < numMarkers[i]; j++) {
					markerPos[count].push_back(bufferRecv[i][j*L_DIMS]);
					markerPos[count].push_back(bufferRecv[i][j*L_DIMS+1]);
#if (L_DIMS == 3)
					markerPos[count].push_back(bufferRecv[i][j*L_DIMS+2]);
#endif
					count++;
				}
			}

			break;
		}
	}
}

///	\brief Unpacks the receive buffer back into IBInfo class.
///	\param	type			type of communication.
///	\param	bufferRecv		the receiver buffer that is being unpacked.
///	\param	numMarkers		where the buffer is being unpacked to
void ObjectManager::ibm_mpi_unpack(eIBInfoType type, std::vector<int> &bufferRecv, int bufferSendSize, std::vector<int> &rankID, std::vector<int> &numMarkers) {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	switch (type)
	{
	case eBodyNumMarkers:
		{

			// Loop through each rank
			for (int i = 0; i < mpim->num_ranks; i++) {

				// Check which ranks it should be receiving non-zero number of markers from
				if (bufferRecv[i*bufferSendSize] != mpim->my_rank && bufferRecv[i*bufferSendSize+1] > 0) {
					rankID.push_back(bufferRecv[i*bufferSendSize]);
					numMarkers.push_back(bufferRecv[i*bufferSendSize+1]);
				}
			}

			break;
		}
	}
}

// ****************************************************************************
