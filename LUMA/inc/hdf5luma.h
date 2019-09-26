/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

#ifndef HDFLUMA_H
#define HDFLUMA_H

#if (defined L_BUILD_FOR_MPI && !defined H5_HAVE_PARALLEL)
	#define H5_HAVE_PARALLEL	///< Enable parallel HDF5
#endif

#include "stdafx.h"
#include "hdf5.h"	// Load C API

#define H5_BUILT_AS_DYNAMIC_LIB
#define HDF5_EXT_ZLIB
#define HDF5_EXT_SZIP


//***************************************************************************//
/// \brief	Helper method to write out using HDF5.
///
///			Automatically selects the correct slab arrangement and buffers the 
///			data accordingly before writing to structured file.
/// \param	memspace		memory dataspace id.
/// \param	filespace		file dataspace id.
/// \param	dataset_id		dataset id.
/// \param	slab_type		slab type enum.
/// \param	g				pointer to grid which we are writing out.
/// \param	data			pointer to the start of the array to be written.
/// \param	hdf_datatype	HDF5 datatype being written.
///	\param	TL_present		pointer to array of flags indicating whether a lower TL is 
///							present on this grid in given direction so offset in 
///							file can be computed.
/// \param	TL_thickness	the thickness of the TL on this grid level in local lattice units.
///	\param	minEdges		pointer to double array containing position of grid edge (from GM).
/// \param	hdf_data		the data structure containing information about local halos.
template <typename T>
void hdf5_writeDataSet(hid_t& memspace, hid_t& filespace, hid_t& dataset_id,
	eHdf5SlabType slab_type, GridObj *g, T *data, hid_t hdf_datatype,
	bool *TL_present, int TL_thickness, double *minEdges,
	HDFstruct hdf_data) {


	// Writable region indicies from the MPIM
	int i_start = hdf_data.i_start;
	int i_end = hdf_data.i_end;
	int j_start = hdf_data.j_start;
	int j_end = hdf_data.j_end;
	int k_start = hdf_data.k_start;
	int k_end = hdf_data.k_end;

	// Index for grid manager quantities
	int idx;
	if (g->level == 0) idx = 0;
	else idx = g->level + g->region_number * L_NUM_LEVELS;

	// Create buffer (block of data we aim to write from this process)
	T *buffer = (T*)malloc(hdf_data.writable_data_count * sizeof(T));

	// DEBUG //
#ifdef L_HDF_DEBUG
	*GridUtils::logfile << "Writing...Writable data size = " 
		<< (i_end - i_start + 1) << "," 
		<< (j_end - j_start + 1) << 
#if (L_DIMS == 3)
		"," << (k_end - k_start + 1) << 
#endif

#ifdef L_BUILD_FOR_MPI
		" = " << hdf_data.writable_data_count << std::endl;
#else
		" = " << (i_end - i_start + 1) * (j_end - j_start + 1) 
#if (L_DIMS == 3)
		* (k_end - k_start + 1)
#endif		
		<< std::endl;
#endif
#endif

	// Create status
	herr_t status;

	// Set slice counter
	int i = i_start;

	// Create property list
	hid_t plist_id = static_cast<hid_t>(NULL);
#ifdef L_BUILD_FOR_MPI
	// Create property template for parallel dataset
	plist_id = H5Pcreate(H5P_DATASET_XFER);

	/* Set data access mode (collective or independent I/O)
	 * Collective IO requires the same number of calls to be made by each MPI
	 * process or MPI I/O will hang. */
	status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Set file access mode failed: " << status << std::endl;
#else
	// Serial dataset
	plist_id = H5P_DEFAULT;
#endif

	/* Hyperslab variables:
	 * offset	= where to start reading/writing within a dataspace
	 * block	= the size of a block in the pattern
	 * count	= how many times the block is repeated in the pattern
	 * stride	= number of elements between start of one block and next */

	// File hyperslab variables
	hsize_t f_offset[L_DIMS], f_block[L_DIMS], f_count[L_DIMS], f_stride[L_DIMS];

	// Memory hyperslab variables (for strided copy)
	size_t m_count, m_stride, m_offset, m_block;

	// Get starting positions in file space
#ifdef L_BUILD_FOR_MPI

	/* Get global offsets for start of file space from the number of cells 
	 * between the origin and the first writable cell.
	 * Correct the offset due to TL presence as TL is not written out. */
	f_offset[0] = static_cast<int>(std::round((g->XPos[i_start] - minEdges[eXDirection] - (g->dh / 2.0)) / g->dh)) 
		- TL_present[eXDirection] * TL_thickness;
	f_offset[1] = static_cast<int>(std::round((g->YPos[j_start] - minEdges[eYDirection] - (g->dh / 2.0)) / g->dh))
		- TL_present[eYDirection] * TL_thickness;
#if (L_DIMS == 3)
	f_offset[2] = static_cast<int>(std::round((g->ZPos[k_start] - minEdges[eZDirection] - (g->dh / 2.0)) / g->dh))
		- TL_present[eZDirection] * TL_thickness;
#endif

#else
	// In serial, only a single process so start writing at the beginning of the file
	for (int d = 0; d < L_DIMS; d++) f_offset[d] = 0;

#endif // L_BUILD_FOR_MPI

	// Block size based on local writable data
	f_block[0] = i_end - i_start + 1;
	f_block[1] = j_end - j_start + 1;

	f_count[0] = 1;
	f_count[1] = 1;

	f_stride[0] = f_block[0];
	f_stride[1] = f_block[1];

#if (L_DIMS == 3)
	f_block[2] = k_end - k_start + 1;
	f_count[2] = 1;
	f_stride[2] = f_block[2];
#endif

	// DEBUG //
#ifdef L_HDF_DEBUG
#if (L_DIMS == 3)
	*GridUtils::logfile << "f_offset = (" << f_offset[0] << " " << f_offset[1] << " " << f_offset[2] << ")" << std::endl;
	*GridUtils::logfile << "f_stride = (" << f_stride[0] << " " << f_stride[1] << " " << f_stride[2] << ")" << std::endl;
	*GridUtils::logfile << "f_count = (" << f_count[0] << " " << f_count[1] << " " << f_count[2] << ")" << std::endl;
	*GridUtils::logfile << "f_block = (" << f_block[0] << " " << f_block[1] << " " << f_block[2] << ")" << std::endl;
#else
	*GridUtils::logfile << "f_offset = (" << f_offset[0] << " " << f_offset[1] << ")" << std::endl;
	*GridUtils::logfile << "f_stride = (" << f_stride[0] << " " << f_stride[1] << ")" << std::endl;
	*GridUtils::logfile << "f_count = (" << f_count[0] << " " << f_count[1] << ")" << std::endl;
	*GridUtils::logfile << "f_block = (" << f_block[0] << " " << f_block[1] << ")" << std::endl;
#endif
#endif
	

	// Select filespace slab
	status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, f_offset, f_stride, f_count, f_block);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Selection of file space hyperslab failed: " << status << std::endl;


	switch (slab_type)
	{

	case eScalar:

		/* One variable per grid site */

#if (L_DIMS == 3)
		// Copy data a 2D slice at a time for 3D cases
		for (i = i_start; i <= i_end; i++)
#endif
		{

			// Get memory space slab parameters
#if (L_DIMS == 3)
			m_offset = k_start + j_start * g->K_lim + i * g->M_lim * g->K_lim;
			m_block = k_end - k_start + 1;
			m_count = j_end - j_start + 1;
			m_stride = g->K_lim;
#else
			m_offset = j_start + i_start * g->M_lim;	// 1D starting index
			m_block = j_end - j_start + 1;			// Size of block in 1D pattern
			m_count = i_end - i_start + 1;			// Number of blocks in pattern
			m_stride = g->M_lim;					// Number of elements between start of two blocks

#endif
			// Copy slab of memory to buffer
#if (L_DIMS == 3)
			size_t buffer_offset = (i - i_start) * m_block * m_count;
#else
			size_t buffer_offset = 0;
#endif
			GridUtils::stridedCopy(
				buffer, data, 
				m_block, m_offset, m_stride, m_count,
				buffer_offset
				);

		}	// End 2D slice loop
		break;



	case eVector:

		/* L_DIMS variables per grid site */

#if (L_DIMS == 3)
		// Copy data a 2D slice at a time for 3D cases
		for (int i = i_start; i <= i_end; i++)
#endif
		{


			// Loop through "columns" of data and pattern fastest dimension
#if (L_DIMS == 3)
			for (int j = j_start; j <= j_end; j++)
#else
			for (int i = i_start; i <= i_end; i++)
#endif
			{

				// Get memory space slab parameters
#if (L_DIMS == 3)
				m_offset = 0 + k_start * L_DIMS + 
					j * L_DIMS * g->K_lim + i * L_DIMS * g->M_lim * g->K_lim;
				m_block = 1;
				m_count = k_end - k_start + 1;
				m_stride = L_DIMS;
#else
				m_offset = 0 + j_start * L_DIMS + i * L_DIMS * g->M_lim;	// Memory offset handled by incoming pointer
				m_block = 1;
				m_count = j_end - j_start + 1;
				m_stride = L_DIMS;
#endif				
				// Strided copy
#if (L_DIMS == 3)
				size_t buffer_offset = (j - j_start) * m_count + 
					(i - i_start) * (j_end - j_start + 1) * m_count;
#else
				size_t buffer_offset = (i - i_start) * m_count;
#endif
				GridUtils::stridedCopy(
					buffer, data,
					m_block, m_offset, m_stride, m_count,
					buffer_offset
					);
			}
		
		} // End of 2D slice loop
		break;



	case eProductVector:

		/* 3*L_DIMS-3 variables per grid site */


#if (L_DIMS == 3)
		// Copy data a 2D slice at a time for 3D cases
		for (int i = i_start; i <= i_end; i++)
#endif
		{

			// Loop through "columns" of data and pattern fastest dimension
#if (L_DIMS == 3)
			for (int j = j_start; j <= j_end; j++)
#else
			for (int i = i_start; i <= i_end; i++)
#endif
			{

				// Get memory space slab parameters
#if (L_DIMS == 3)
				m_offset = 0 + k_start * (3 * L_DIMS - 3) + 
					j * (3 * L_DIMS - 3) * g->K_lim + i * (3 * L_DIMS - 3) * g->M_lim * g->K_lim;
				m_block = 1;
				m_count = k_end - k_start + 1;
				m_stride = 3 * L_DIMS - 3;
#else
				m_offset = 0 + j_start * (3 * L_DIMS - 3) + i * (3 * L_DIMS - 3) * g->M_lim;
				m_block = 1;
				m_count = j_end - j_start + 1;
				m_stride = 3 * L_DIMS - 3;
#endif
				// Strided copy
#if (L_DIMS == 3)
				size_t buffer_offset = (j - j_start) * m_count +
					(i - i_start) * (j_end - j_start + 1) * m_count;
#else
				size_t buffer_offset = (i - i_start) * m_count;
#endif
				GridUtils::stridedCopy(
					buffer, data,
					m_block, m_offset, m_stride, m_count,
					buffer_offset
					);
			}

		} // End of 2D slice loop
		break;


	case ePosX:

		for (int i = i_start; i <= i_end; i++) {
			for (int j = j_start; j <= j_end; j++) {
#if (L_DIMS == 3)
				for (int k = k_start; k <= k_end; k++) {

					memcpy(
						buffer +
						(k - k_start) +
						(j - j_start) * (k_end - k_start + 1) +
						(i - i_start) * (k_end - k_start + 1) * (j_end - j_start + 1),
						data + i, sizeof(T)
						);
				}
#else

				memcpy(
					buffer + 
					(j - j_start) + 
					(i - i_start) * (j_end - j_start + 1),
					data + i, sizeof(T)
					);
#endif
			}
		}
		break;

	case ePosY:

		for (int i = i_start; i <= i_end; i++) {
#if (L_DIMS == 3)
			for (int j = j_start; j <= j_end; j++) {
				for (int k = k_start; k <= k_end; k++) {

					memcpy(
						buffer +
						(k - k_start) +
						(j - j_start) * (k_end - k_start + 1) +
						(i - i_start) * (k_end - k_start + 1) * (j_end - j_start + 1),
						data + j, sizeof(T)
						);
				}
			}
#else
			memcpy(
				buffer + 
				(i - i_start) * (j_end - j_start + 1),
				data + j_start, sizeof(T) * (j_end - j_start + 1)
				);
#endif
		}
		break;

	case ePosZ:

		for (int i = i_start; i <= i_end; i++) {
			for (int j = j_start; j <= j_end; j++) {

				memcpy(
					buffer + 
					(j - j_start) * (k_end - k_start + 1) + 
					(i - i_start) * (k_end - k_start + 1) * (j_end - j_start + 1),
					data + k_start, sizeof(T) * (k_end - k_start + 1)
					);
			}
		}
		break;

	}	// End switch on slab_type

	// DEBUG //
#ifdef L_HDF_DEBUG
	*GridUtils::logfile << "Lev " << g->level << " writing slab type = " << std::to_string(slab_type) << "...";
#endif

	// Write data
	status = H5Dwrite(dataset_id, hdf_datatype, memspace, filespace, plist_id, buffer);
	if (status != 0) {
		*GridUtils::logfile << "HDF5 ERROR: Write data failed: " << status << std::endl;
		H5Eprint(H5E_DEFAULT, stderr);
	}
	else
#ifdef L_HDF_DEBUG
		*GridUtils::logfile << "Write Successful." << std::endl;
#endif
	H5Sselect_none(filespace);



	// Close property list
	status = H5Pclose(plist_id);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close file access mode list failed: " << status << std::endl;

	// Free the buffer
	free(buffer);

};

#endif

#ifdef L_BUILD_FOR_MPI

// ************************************************************************** //
/// \brief	Checks for consistent file space calculations before creating HDF5 
///			datasets.
///
///			This was created as the H5Fclose() method can hang if a dataset is
///			created in a file where all ranks contributing to the file do not
///			agree on the file size. It was hard to debug.
///
///	\param	dimsf	pointer to the file space dimensions on this rank.
///	\param	comm	communicator.
void hdf_checkFileSpace(hsize_t * dimsf, MPI_Comm& comm)
{
	// Get size of the communicator
	int numProc;
	MPI_Comm_size(comm, &numProc);

	// Create send and recv buffers
	std::vector<int> sendBuf(numProc * L_DIMS);
	for (size_t i = 0; i < sendBuf.size(); i += L_DIMS)
	{
		for (int d = 0; d < L_DIMS; d++)
		{
			sendBuf[i + d] = static_cast<int>(dimsf[d]);
		}
	}
	std::vector<int> recvBuf(sendBuf);

	// Send all to all
	MPI_Alltoall(&sendBuf[0], L_DIMS, MPI_INT, &recvBuf[0], L_DIMS, MPI_INT, comm);

	// Check the recv buffer to make sure all the value sets are equal
	for (size_t i = 0; i < recvBuf.size(); i += L_DIMS)
	{
		for (int d = 0; d < L_DIMS; d++)
		{
			if (recvBuf[i + d] != sendBuf[i + d])
			{
				L_ERROR("Filespace received = " + std::to_string(recvBuf[i + d]) +
					", Filespace sent = " + std::to_string(sendBuf[i + d]), GridUtils::logfile);
			}
		}
	}

}

#endif