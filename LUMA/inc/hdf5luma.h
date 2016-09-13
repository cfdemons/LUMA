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

#ifdef L_BUILD_FOR_MPI
	#define H5_HAVE_PARALLEL
	#include "hdf5.h"	// Load C API
#else
	#include "H5Cpp.h"	// Load C++ API
	#ifndef H5_NO_NAMESPACE
		using namespace H5;
	#endif
#endif

#define H5_BUILT_AS_DYNAMIC_LIB
#define HDF5_EXT_ZLIB
#define HDF5_EXT_SZIP

enum ePhdf5SlabType {
	eScalar,		/* 2/3D data	-- One variable per grid site			*/
	eVector,		/* 2/3D data	-- L_dims variables per grid site		*/
	eProductVector,	/* 1D data		-- 3*L_dims-3 variables per grid site	*/
	ePosX,			/* 1D data		-- Single L_dim vector per dimension	*/
	ePosY,			/* 1D data		-- Single L_dim vector per dimension	*/
	ePosZ			/* 1D data		-- Single L_dim vector per dimension	*/
};


#pragma once	// Only define the structure once to prevent redefinition
struct phdf5_struct {
	/* Structure storing indices and boolean flags 
	 * related to the presence of halo regions in 
	 * MPI PHDF5 writing */
	int i_start;
	int i_end;
	int j_start;
	int j_end;
	int k_start;
	int k_end;
	int halo_min;
	int halo_max;
};

// ***************************************************************************//
// Helper method to write out using PHDF5 which automatically selects the 
// correct slab arrangement. Templated so in header.
template <typename T>
void phdf5_writeDataSet(hid_t& memspace, hid_t& filespace, hid_t& fileset,
	ePhdf5SlabType slab_type, phdf5_struct hdf_data,
	int N_lim, int M_lim, int K_lim, 
	GridObj *g, T *data, hid_t hdf_datatype) {

	// Create buffer (block of data we aim to write from this process)
	const unsigned int buffer_size =
		(hdf_data.i_end - hdf_data.i_start + 1) *
		(hdf_data.j_end - hdf_data.j_start + 1) *
		(hdf_data.k_end - hdf_data.k_start + 1);
	//std::cout << (hdf_data.i_end - hdf_data.i_start + 1) << "," << (hdf_data.j_end - hdf_data.j_start + 1) << "," << (hdf_data.k_end - hdf_data.k_start + 1) << std::endl;
	//std::cout << "BS " << buffer_size << std::endl;
	T *buffer = (T*)malloc(buffer_size * sizeof(T));

	// Create status
	herr_t status;

	// Set slice counter
	int i = hdf_data.i_start;

	// Create property template for parallel dataset
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

	/* Set data access mode (collective or independent I/O)
	 * Collective IO requires the same number of calls to be made by each MPI
	 * process or MPI I/O will hang. */
	status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Set file access mode failed: " << status << std::endl;

	/* Hyperslab variables:
	 * offset	= where to start reading/writing within a dataspace
	 * block	= the size of a block in the pattern
	 * count	= how many times the block is repeated in the pattern
	 * stride	= number of elements between start of one block and next */
	hsize_t f_count[L_dims], f_stride[L_dims], f_offset[L_dims], f_block[L_dims];

	// Strided copy variables
	size_t m_count, m_stride, m_offset, m_block;

	// Get global offsets for start of file space from the indexing vector
	f_offset[0] = g->XInd[hdf_data.i_start];
	f_offset[1] = g->YInd[hdf_data.j_start];

	f_block[0] = hdf_data.i_end - hdf_data.i_start + 1;
	f_block[1] = hdf_data.j_end - hdf_data.j_start + 1;

	f_count[0] = 1; f_count[1] = 1;
	f_stride[0] = f_block[0]; f_stride[1] = f_block[1];
#if (L_dims == 3)
	f_offset[2] = g->ZInd[hdf_data.k_start];
	f_block[2] = hdf_data.k_end - hdf_data.k_start + 1;
	f_count[2] = 1;
	f_stride[2] = f_block[2];
#endif

	// Select filespace slab
	status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, f_offset, f_stride, f_count, f_block);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Selection of file space hyperslab failed: " << status << std::endl;


	switch (slab_type)
	{

	case eScalar:

		/* One variable per grid site */

#if (L_dims == 3)
		// Copy data a 2D slice at a time for 3D cases
		for (i = hdf_data.i_start; i <= hdf_data.i_end; i++)
#endif
		{

			// Get memory space slab parameters from halo descriptors
#if (L_dims == 3)
			m_offset = hdf_data.k_start + hdf_data.j_start * K_lim + i * M_lim * K_lim;
			m_block = hdf_data.k_end - hdf_data.k_start + 1;
			m_count = hdf_data.j_end - hdf_data.j_start + 1;
			m_stride = m_block + (hdf_data.halo_min + hdf_data.halo_max) * static_cast<hsize_t>(pow(2, g->level));
#else
			m_offset = hdf_data.j_start + hdf_data.i_start * M_lim;
			m_block = hdf_data.j_end - hdf_data.j_start + 1;
			m_count = hdf_data.i_end - hdf_data.i_start + 1;
			m_stride = m_block + (hdf_data.halo_min + hdf_data.halo_max) * static_cast<hsize_t>(pow(2, g->level));
#endif
			// Copy slab of memory to buffer
#if (L_dims == 3)
			size_t buffer_offset = (i - hdf_data.i_start) * m_block * m_count;
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

		/* L_dims variables per grid site */

#if (L_dims == 3)
		// Copy data a 2D slice at a time for 3D cases
		for (int i = hdf_data.i_start; i <= hdf_data.i_end; i++)
#endif
		{


			// Loop through "columns" of data and pattern fastest dimension
#if (L_dims == 3)
			for (int j = hdf_data.j_start; j <= hdf_data.j_end; j++)
#else
			for (int i = hdf_data.i_start; i <= hdf_data.i_end; i++)
#endif
			{

				// Get memory space slab parameters from halo descriptors
#if (L_dims == 3)
				m_offset = 0 + hdf_data.k_start * L_dims + 
					j * L_dims * K_lim + i * L_dims * M_lim * K_lim;
				m_block = 1;
				m_count = hdf_data.k_end - hdf_data.k_start + 1;
				m_stride = L_dims;
#else
				m_offset = 0 + hdf_data.j_start * L_dims + i * L_dims * M_lim;	// Memory offset handled by incoming pointer
				m_block = 1;
				m_count = hdf_data.j_end - hdf_data.j_start + 1;
				m_stride = L_dims;
#endif				
				// Strided copy
#if (L_dims == 3)
				size_t buffer_offset = (j - hdf_data.j_start) * m_count + 
					(i - hdf_data.i_start) * (hdf_data.j_end - hdf_data.j_start + 1) * m_count;
#else
				size_t buffer_offset = (i - hdf_data.i_start) * m_count;
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

		/* 3*L_dims-3 variables per grid site */


#if (L_dims == 3)
		// Copy data a 2D slice at a time for 3D cases
		for (int i = hdf_data.i_start; i <= hdf_data.i_end; i++)
#endif
		{

			// Loop through "columns" of data and pattern fastest dimension
#if (L_dims == 3)
			for (int j = hdf_data.j_start; j <= hdf_data.j_end; j++)
#else
			for (int i = hdf_data.i_start; i <= hdf_data.i_end; i++)
#endif
			{

				// Get memory space slab parameters from halo descriptors
#if (L_dims == 3)
				m_offset = 0 + hdf_data.k_start * (3 * L_dims - 3) + 
					j * (3 * L_dims - 3) * K_lim + i * (3 * L_dims - 3) * M_lim * K_lim;
				m_block = 1;
				m_count = hdf_data.k_end - hdf_data.k_start + 1;
				m_stride = 3 * L_dims - 3;
#else
				m_offset = 0 + hdf_data.j_start * (3 * L_dims - 3) + i * (3 * L_dims - 3) * M_lim;
				m_block = 1;
				m_count = hdf_data.j_end - hdf_data.j_start + 1;
				m_stride = 3 * L_dims - 3;
#endif
				// Strided copy
#if (L_dims == 3)
				size_t buffer_offset = (j - hdf_data.j_start) * m_count +
					(i - hdf_data.i_start) * (hdf_data.j_end - hdf_data.j_start + 1) * m_count;
#else
				size_t buffer_offset = (i - hdf_data.i_start) * m_count;
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

		for (int i = hdf_data.i_start; i <= hdf_data.i_end; i++) {
			for (int j = hdf_data.j_start; j <= hdf_data.j_end; j++) {
#if (L_dims == 3)
				for (int k = hdf_data.k_start; k <= hdf_data.k_end; k++) {

					memcpy(
						buffer +
						(k - hdf_data.k_start) +
						(j - hdf_data.j_start) * (hdf_data.k_end - hdf_data.k_start + 1) +
						(i - hdf_data.i_start) * (hdf_data.k_end - hdf_data.k_start + 1) * (hdf_data.j_end - hdf_data.j_start + 1),
						data + i, sizeof(T)
						);
				}
#else

				memcpy(
					buffer + 
					(j - hdf_data.j_start) + 
					(i - hdf_data.i_start) * (hdf_data.j_end - hdf_data.j_start + 1),
					data + i, sizeof(T)
					);
#endif
			}
		}
		break;

	case ePosY:

		for (int i = hdf_data.i_start; i <= hdf_data.i_end; i++) {
#if (L_dims == 3)
			for (int j = hdf_data.j_start; j <= hdf_data.j_end; j++) {
				for (int k = hdf_data.k_start; k <= hdf_data.k_end; k++) {

					memcpy(
						buffer +
						(k - hdf_data.k_start) +
						(j - hdf_data.j_start) * (hdf_data.k_end - hdf_data.k_start + 1) +
						(i - hdf_data.i_start) * (hdf_data.k_end - hdf_data.k_start + 1) * (hdf_data.j_end - hdf_data.j_start + 1),
						data + j, sizeof(T)
						);
				}
			}
#else
			memcpy(
				buffer + 
				(i - hdf_data.i_start) * (hdf_data.j_end - hdf_data.j_start + 1),
				data + hdf_data.j_start, sizeof(T) * (hdf_data.j_end - hdf_data.j_start + 1)
				);
#endif
		}
		break;

	case ePosZ:

		for (int i = hdf_data.i_start; i <= hdf_data.i_end; i++) {
			for (int j = hdf_data.j_start; j <= hdf_data.j_end; j++) {

				memcpy(
					buffer + 
					(j - hdf_data.j_start) * (hdf_data.k_end - hdf_data.k_start + 1) + 
					(i - hdf_data.i_start) * (hdf_data.k_end - hdf_data.k_start + 1) * (hdf_data.j_end - hdf_data.j_start + 1),
					data + hdf_data.k_start, sizeof(T) * (hdf_data.k_end - hdf_data.k_start + 1)
					);
			}
		}
		break;

	}	// End switch on slab_type


	// Write data
	status = H5Dwrite(fileset, hdf_datatype, memspace, filespace, plist_id, buffer);
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

}