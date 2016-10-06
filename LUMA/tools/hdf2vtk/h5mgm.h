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

/* H5 Multi-Grid Merge Tool for post-processing HDF5 files written by LUMA */

#define H5MGM_VERSION "0.1-alpha"

#include "hdf5.h"

#define H5_BUILT_AS_DYNAMIC_LIB
#define HDF5_EXT_ZLIB
#define HDF5_EXT_SZIP

#include <string>
#include <iostream>
#include <vector>

// Typing enumeration from LUMA
enum eType
{
	eSolid,					// Rigid, solid site
	eFluid,					// Fluid site
	eRefined,				// Fluid site which is represented on a finer grid
	eTransitionToCoarser,	// Fluid site coupled to a coarser grid
	eTransitionToFiner,		// Fluid site coupled to a finer grid
	eBFL,					// Site containing a BFL marker
	eSymmetry,				// Symmetry boundary
	eInlet,					// Inlet boundary
	eOutlet,				// Outlet boundary
	eRefinedSolid,			// Rigid, solid site represented on a finer grid
	eRefinedSymmetry,		// Symmtery boundary represented on a finer grid
	eRefinedInlet			// Inlet site represented on a finer grid
};

/* Templated function to map data from a coarse grid to a fine grid buffer*/
template<typename T>
herr_t bufferData(std::string VAR, std::string TIME_STRING, std::string OUT_FILE_NAME,
	int *gridsize, unsigned int levels, unsigned int regions,
	hid_t& output_fid, hid_t H5type, int*& LatTyp, T* dummy) {

	// Declarations
	herr_t status = 0;
	T *in_buffer_start = nullptr;
	T *in_buffer = nullptr;
	hsize_t dims_output[1];
	hsize_t dims_input[1];
	std::string variable_string = TIME_STRING + VAR;
	size_t in_buffer_size = 0;
	size_t out_count = 0;

	// Loop over the Grids
	for (unsigned int lev = 0; lev < levels; lev++) {
		for (unsigned int reg = 0; reg < regions; reg++) {

			// Declarations
			hid_t input_fid = NULL;
			hid_t input_aid = NULL;
			hid_t input_did = NULL;
			hid_t input_sid = NULL;

			// L0 doesn't have different regions
			if (lev == 0 && reg != 0) continue;

			// Construct input file name
			std::string IN_FILE_NAME("./hdf_R" + std::to_string(reg) + "N" + std::to_string(lev) + ".h5");
			input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (input_fid <= 0) std::cout << "HDF5 ERROR: Cannot open input file!" << std::endl;

			// Get local grid size
			input_aid = H5Aopen(input_fid, "GridSize", H5P_DEFAULT);
			if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
			status = H5Aread(input_aid, H5T_NATIVE_INT, gridsize);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
			status = H5Aclose(input_aid);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

			if (lev == 0 && reg == 0) {
				// Allocate input buffer for the first grid
				in_buffer = (T*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(T));
				if (in_buffer == NULL) {
					std::cout << "Not enough Memory!!!!" << std::endl;
					exit(EXIT_FAILURE);
				}
				in_buffer_start = in_buffer;	// Store pointer to start of array
				in_buffer_size = gridsize[0] * gridsize[1] * gridsize[2];	// Store size
			}
			else {
				// Reallocate input buffer to provide space for sub-grid grid to be appended					
				in_buffer_start = (T*)realloc(in_buffer_start, (in_buffer_size + (gridsize[0] * gridsize[1] * gridsize[2])) * sizeof(T));
				if (in_buffer_start == NULL) {
					std::cout << "Out of Memory!!!!" << std::endl;
					exit(EXIT_FAILURE);
				}
				in_buffer = in_buffer_start + in_buffer_size;	// Set pointer to correct place in new buffer
				in_buffer_size += gridsize[0] * gridsize[1] * gridsize[2];	// Update size
			}

			// Open input dataset
			input_did = H5Dopen(input_fid, variable_string.c_str(), H5P_DEFAULT);
			if (input_did <= 0) {
				std::cout << "HDF5 ERROR: Cannot open input dataset -- no output data for " <<
					TIME_STRING << std::endl;
				free(in_buffer_start);
				return -1;
			}

			// Create input dataspace
			dims_input[0] = gridsize[0] * gridsize[1] * gridsize[2];
			input_sid = H5Screate_simple(1, dims_input, NULL);
			if (input_did <= 0) std::cout << "HDF5 ERROR: Cannot create input dataspace!" << std::endl;

			// Read input datset into new section of 1D buffer
			status = H5Dread(input_did, H5type, input_sid, H5S_ALL, H5P_DEFAULT, in_buffer);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot read input dataset!" << std::endl;

			// Close input dataset
			status = H5Dclose(input_did);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot close input dataset!" << std::endl;

			// Close input dataset
			status = H5Sclose(input_sid);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot close input dataspace!" << std::endl;

			// Close input file
			status = H5Fclose(input_fid);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot close input file!" << std::endl;

		}

	}

	// Store LatTyp matrix on first time through before LatTyp has been assigned
	if (LatTyp == nullptr)
		LatTyp = reinterpret_cast<int*>(in_buffer_start);	// Will never be called if in_buffer_start is not an integer
	
	// Allocate max output space
	T *out_buffer = (T*)malloc(in_buffer_size * sizeof(T));

	// Loop over typing and decide what to store
	for (size_t i = 0; i < in_buffer_size; i++) {

		if (LatTyp[i] != eTransitionToCoarser || LatTyp[i] == eTransitionToFiner || LatTyp[i] == eRefined) {

			// Assign value
			*out_buffer = LatTyp[i];

			// Move pointer and increment total
			out_buffer++;
			out_count++;

		}
	}

	// Rewind pointer
	out_buffer -= out_count;

	// Data space and set ids
	hid_t output_sid = NULL;
	hid_t output_did = NULL;

	// Create dataspace
	dims_output[0] = out_count;
	output_sid = H5Screate_simple(1, dims_output, NULL);

	// Open output file on first grid pass only, already open on other variables
	if (output_fid == NULL) {
		output_fid = H5Fcreate(OUT_FILE_NAME.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (output_fid <= 0) std::cout << "HDF5 ERROR: Cannot create output file!" << std::endl;
	}

	// Create output dataset and write
	output_did = H5Dcreate(output_fid, VAR.c_str(), H5type, output_sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (output_did <= 0) std::cout << "HDF5 ERROR: Cannot create output dataset!" << std::endl;

	// Write variable to output file
	status = H5Dwrite(output_did, H5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_buffer);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot write output dataset!" << std::endl;

	// Close output datasets		
	status = H5Dclose(output_did);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close output dataset!" << std::endl;

	// Close output dataspace
	status = H5Sclose(output_sid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close output dataspace!" << std::endl;

	/* Free buffers.
	 * Only free the input buffer if it wasn't the LatTyp as we need to keep 
	 * this for other variables */
	if (LatTyp != reinterpret_cast<int*>(in_buffer_start))
		free(in_buffer_start);
	free(out_buffer);

	return status;
}

/* End of header */