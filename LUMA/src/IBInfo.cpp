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
#include "../inc/IBInfo.h"
#include "../inc/IBBody.h"
#include "../inc/IBMarker.h"

IBInfo::IBInfo()
{
	// Default constructor
}

/// \brief Custom constructor for different types of message containers.
///	\param	iBody	pointer to the iBody being packed.
///	\param	type	the type of container to be created.
IBInfo::IBInfo(IBBody *iBody, int m, eIBInfoType type)
{

	switch (type)
	{

	case eIBDeltaSum:

		// Pack voxel_centres for each marker and sum of delta functions
		// First support point should be the closest
		markerX = iBody->markers[m].position[eXDirection];
		markerY = iBody->markers[m].position[eYDirection];
		markerZ = iBody->markers[m].position[eZDirection];

		// Delta values
		//deltaVals = iBody->markers[m].deltaval;

		break;

		// Add other cases here

	}
}

///	\brief Maps a version of the IBInfo structure to an MPI_Struct datatype.
///	\param	type	type of container you want to map.
///	\return	handle to the MPI struct data.
int IBInfo::mapToMpiStruct(eIBInfoType type, MPI_Datatype *mpi_struct_type)
{

	switch (type)
	{

	case eIBDeltaSum:

		// Create the necessary info for creating the MPI_Type_struct
		// Number of items to send
		int count = 3;

		// Size of each array that is being sent
		int array_of_block_lengths[count];
		array_of_block_lengths[0] = 1;
		array_of_block_lengths[1] = 1;
		array_of_block_lengths[2] = 1;
		//array_of_block_lengths[3] = deltaVals.size();

		// Get the address of each element that is being sent
		MPI_Aint markX_addr, markY_addr, markZ_addr, markDelta_addr;
		MPI_Get_address(&(this->markerX), &markX_addr);
		MPI_Get_address(&(this->markerY), &markY_addr);
		MPI_Get_address(&(this->markerZ), &markZ_addr);
		//MPI_Get_address(&(this->deltaVals.front()), &markDelta_addr);

		// Offsets of the addresses
		MPI_Aint array_of_displacements[count];
		array_of_displacements[0] = 0;
		array_of_displacements[1] = markY_addr - markX_addr;
		array_of_displacements[2] = markZ_addr - markX_addr;
		//array_of_displacements[3] = markDelta_addr - markX_addr;

//		std::cout << markX_addr << " " << markY_addr << " " << markZ_addr << " " << markDelta_addr << std::endl;
//		std::cout << array_of_displacements[0] << " " << array_of_displacements[1] << " " << array_of_displacements[2] << " " << array_of_displacements[3] << std::endl;
//		std::cout << array_of_block_lengths[0] << " " << array_of_block_lengths[1] << " " << array_of_block_lengths[2] << " " << array_of_block_lengths[3] << std::endl;


		// Specify the datatypes
		MPI_Datatype array_of_types[count];
		array_of_types[0] = MPI_DOUBLE;
		array_of_types[1] = MPI_DOUBLE;
		array_of_types[2] = MPI_DOUBLE;
		//array_of_types[3] = MPI_DOUBLE;

		// Create and commit the new datatype
		MPI_Type_create_struct(count, array_of_block_lengths, array_of_displacements, array_of_types, mpi_struct_type);
		MPI_Type_commit(mpi_struct_type);

		break;
	}
}
