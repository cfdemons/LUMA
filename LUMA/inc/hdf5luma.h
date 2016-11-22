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

#ifndef HDFLUMA_H
#define HDFLUMA_H

#if (defined L_BUILD_FOR_MPI && !defined H5_HAVE_PARALLEL)
	#define H5_HAVE_PARALLEL	///< Enable parallel HDF5
#endif

#include "stdafx.h"
#include "hdf5.h"	// Load C API
#include "MpiManager.h"

#define H5_BUILT_AS_DYNAMIC_LIB
#define HDF5_EXT_ZLIB
#define HDF5_EXT_SZIP

/// \enum	eHdf5SlabType
///	\brief	Defines the type of storage arrangement of the variable in memory.
///
///			The write wrapper can then extract the data from memeory and write it
///			to an HDF5 file using a particular hyperslab selection.
enum eHdf5SlabType {
	eScalar,		///< 2/3D data	-- One variable per grid site
	eVector,		///< 2/3D data	-- L_DIMS variables per grid site
	eProductVector,	///< 1D data	-- 3*L_DIMS-3 variables per grid site
	ePosX,			///< 1D data	-- Single L_dim vector per dimension
	ePosY,			///< 1D data	-- Single L_dim vector per dimension
	ePosZ			///< 1D data	-- Single L_dim vector per dimension
};

//***************************************************************************//
/// \brief	Helper method to write out using HDF5.
///
///			Automatically selects the correct slab arrangement and buffers the 
///			data accordingly before writing to structured file.
/// \param memspace		memory dataspace id.
/// \param filespace	file dataspace id.
/// \param dataset_id	dataset id.
/// \param slab_type	slab type enum.
/// \param N_lim		number of X-direction sites on the local grid.
/// \param M_lim		number of Y-direction sites on the local grid.
/// \param K_lim		number of Z-direction sites on the local grid.
/// \param N_mod		number of X-direction sites excluding TL sites.
/// \param M_mod		number of Y-direction sites excluding TL sites.
/// \param K_mod		number of Z-direction sites excluding TL sites.
/// \param g			pointer to grid which we are writing out.
/// \param data			pointer to the start of the array to be written.
/// \param hdf_datatype HDF5 datatype being written.
/// \param TL_thickness the thickness of the TL on this grid level in local lattice units.
/// \param hdf_data		the data structure containing information about local halos.
template <typename T>
void hdf5_writeDataSet(hid_t& memspace, hid_t& filespace, hid_t& dataset_id,
	eHdf5SlabType slab_type,
	int N_lim, int M_lim, int K_lim,
	int N_mod, int M_mod, int K_mod,
	GridObj *g, T *data, hid_t hdf_datatype,
	int TL_thickness, MpiManager::phdf5_struct hdf_data) {


	
#ifdef L_BUILD_FOR_MPI

	// Writable region indicies
	int i_start = hdf_data.i_start;
	int i_end = hdf_data.i_end;
	int j_start = hdf_data.j_start;
	int j_end = hdf_data.j_end;
	int k_start = hdf_data.k_start;
	int k_end = hdf_data.k_end;

	/* HALO DESCRIPTORS WILL BE REMOVED IN A FUTURE RELEASE */

	//int halo_min = hdf_data.halo_min;
	//int halo_max = hdf_data.halo_max;

	// Create buffer (block of data we aim to write from this process)
	T *buffer = (T*)malloc(hdf_data.writable_data_count * sizeof(T));

#else

	// Writable region indicies
	int i_start = TL_thickness;
	int i_end = N_lim - TL_thickness - 1;
	int j_start = TL_thickness;
	int j_end = M_lim - TL_thickness - 1;
	int k_start = TL_thickness;
	int k_end = K_lim - TL_thickness - 1;
	
	/* HALO DESCRIPTORS WILL BE REMOVED IN A FUTURE RELEASE */

	//int halo_min = 2;
	//int halo_max = 2;

	// L0 grids do not have TL
	/*if (g->level == 0) {
		halo_min = 0;
		halo_max = 0;
	}*/

	// Create buffer excluding TL
	T *buffer = (T*)malloc(N_mod * M_mod * K_mod * sizeof(T));
#endif

	// DEBUG //
#ifdef L_HDF_DEBUG
	*GridUtils::logfile << "Writable data size = " 
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

	/* Get global offsets for start of file space from the indexing vector.
	 * Correct the offset due to TL presence. i.e. if the start index is element
	 * number 2 we need to write to element 0 in the file since the real domain
	 * will be 2 * TL wider than the file space so elements need shifting back by
	 * a single TL_thickness. */
	f_offset[0] = g->XInd[i_start] - TL_thickness;
	f_offset[1] = g->YInd[j_start] - TL_thickness;
#if (L_DIMS == 3)
	f_offset[2] = g->ZInd[k_start] - TL_thickness;
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
			m_offset = k_start + j_start * K_lim + i * M_lim * K_lim;
			m_block = k_end - k_start + 1;
			m_count = j_end - j_start + 1;
			m_stride = K_lim;
#else
			m_offset = j_start + i_start * M_lim;	// 1D starting index
			m_block = j_end - j_start + 1;			// Size of block in 1D pattern
			m_count = i_end - i_start + 1;			// Number of blocks in pattern
			m_stride = M_lim;						// Number of elements between start of two blocks

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
					j * L_DIMS * K_lim + i * L_DIMS * M_lim * K_lim;
				m_block = 1;
				m_count = k_end - k_start + 1;
				m_stride = L_DIMS;
#else
				m_offset = 0 + j_start * L_DIMS + i * L_DIMS * M_lim;	// Memory offset handled by incoming pointer
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
					j * (3 * L_DIMS - 3) * K_lim + i * (3 * L_DIMS - 3) * M_lim * K_lim;
				m_block = 1;
				m_count = k_end - k_start + 1;
				m_stride = 3 * L_DIMS - 3;
#else
				m_offset = 0 + j_start * (3 * L_DIMS - 3) + i * (3 * L_DIMS - 3) * M_lim;
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
	*GridUtils::logfile << "Lev " << g->level << " WRITING slab type " << std::to_string(slab_type) << std::endl;
#endif

	// Write data
	status = H5Dwrite(dataset_id, hdf_datatype, memspace, filespace, plist_id, buffer);
	if (status != 0) {
		*GridUtils::logfile << "HDF5 ERROR: Write data failed: " << status << std::endl;
		H5Eprint(H5E_DEFAULT, stderr);
	}
	H5Sselect_none(filespace);

	// Close property list
	status = H5Pclose(plist_id);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close file access mode list failed: " << status << std::endl;

	// Free the buffer
	free(buffer);

};

#endif