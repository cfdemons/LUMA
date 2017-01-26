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
void ObjectManager::ibm_mpi_communicate(eIBInfoType type, IBBody *iBody, std::vector<IBInfo> &numMarkers) {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	switch (type)
	{
	case eBodyNumMarkers:
		{
			// Get how many markers from each rank the owning rank should ask for
			std::vector<int> bufferSend;
			IBInfo msg_send(type, iBody, bufferSend);

			// Instantiate the receive buffer
			std::vector<int> bufferRecv;

			// Only resize buffer if this is the owning rank
			if (mpim->my_rank == iBody->owningRank)
				bufferRecv.resize(bufferSend.size() * mpim->num_ranks);

			// Send to owning rank (root process)
			MPI_Gather(&bufferSend.front(), bufferSend.size(), msg_send.mpi_type, &bufferRecv.front(), bufferSend.size(), msg_send.mpi_type, iBody->owningRank, mpim->world_comm);

			// Unpack back into IBInfo class
			if (mpim->my_rank == iBody->owningRank)
				ibm_mpi_unpack(type, bufferRecv, bufferSend.size(), numMarkers);

			break;
		}
	}
}
///	\brief Unpacks the receive buffer back into IBInfo class.
///	\param	type			type of communication.
///	\param	bufferRecv		the receiver buffer that is being unpacked.
///	\param	numMarkers		where the buffer is being unpacked to
void ObjectManager::ibm_mpi_unpack(eIBInfoType type, std::vector<int> &bufferRecv, int bufferSendSize, std::vector<IBInfo> &numMarkers) {

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();

	switch (type)
	{
	case eBodyNumMarkers:
		{

			// Create temporary IBInfo class
			IBInfo tempIBInfo;

			// Loop through each rank
			for (int i = 0; i < mpim->num_ranks; i++) {

				// Check which ranks it should be receiving non-zero number of markers from
				if (bufferRecv[i*bufferSendSize] != mpim->my_rank && bufferRecv[i*bufferSendSize+1] > 0) {
					tempIBInfo.sendingRank = bufferRecv[i*bufferSendSize];
					tempIBInfo.nMarkers = bufferRecv[i*bufferSendSize+1];
					numMarkers.push_back(tempIBInfo);
				}
			}

			break;
		}
	}
}

// ****************************************************************************
