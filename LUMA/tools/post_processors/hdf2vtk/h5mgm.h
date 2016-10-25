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

#define H5MGM_VERSION "0.2.2"

#include "hdf5.h"
#define H5_BUILT_AS_DYNAMIC_LIB
#define HDF5_EXT_ZLIB
#define HDF5_EXT_SZIP

#include <string>
#include <iostream>
#include <vector>

#ifdef _WIN32
	#include <Windows.h>
#endif

#include "vtkVersion.h"
#include "vtkSmartPointer.h"
#include "vtkVoxel.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"

// Unit vectors for node positions on each cell
const int e[3][27] =
{
	{ -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1 },
	{ -1, -1, -1,  0,  0,  0,  1,  1,  1, -1, -1, -1,  0,  0,  0,  1,  1,  1, -1, -1, -1,  0,  0,  0,  1,  1,  1 },
	{ -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1 }
};

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

// Method to read a dataset with a given name into a buffer
template <typename T>
herr_t readDataset(std::string VAR, std::string TIME_STRING, hid_t input_fid, hid_t input_sid, hid_t H5Type, T *buffer) {

	herr_t status;
	std::string variable_string = TIME_STRING + VAR;
	hid_t input_did = H5Dopen(input_fid, variable_string.c_str(), H5P_DEFAULT);
	if (input_did <= 0) {
		std::cout << "HDF5 ERROR: Cannot open input dataset -- no output data for " <<
			TIME_STRING << std::endl;
		exit(EXIT_FAILURE);
	}
	status = H5Dread(input_did, H5Type, input_sid, H5S_ALL, H5P_DEFAULT, buffer);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read input dataset!" << std::endl;
	status = H5Dclose(input_did);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close input dataset!" << std::endl;

	return status;

}

// Method to compile and add arrays of cell data to the mesh
template<typename T, typename vtkT>
herr_t addDataToGrid(std::string VAR, std::string TIME_STRING, 
	int levels, int regions, int* gridsize, 
	vtkSmartPointer<vtkUnstructuredGrid> grid, T *dummy, hid_t H5Type, 
	vtkSmartPointer<vtkT> vtkArray) {

	// Array ID counter
	int count = 0;

	for (int lev = 0; lev < levels; lev++) {
		for (int reg = 0; reg < regions; reg++) {

			// L0 doesn't have different regions
			if (lev == 0 && reg != 0) continue;

			// Construct input file name
			std::string IN_FILE_NAME("./hdf_R" + std::to_string(reg) + "N" + std::to_string(lev) + ".h5");
			hid_t input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (input_fid <= 0) std::cout << "HDF5 ERROR: Cannot open input file!" << std::endl;

			// Get local grid size
			hid_t input_aid = H5Aopen(input_fid, "GridSize", H5P_DEFAULT);
			if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
			herr_t status = H5Aread(input_aid, H5T_NATIVE_INT, gridsize);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
			status = H5Aclose(input_aid);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

			// Allocate space for typing matrix
			int *Type = (int*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(int));
			if (Type == NULL) {
				std::cout << "Not enough Memory!!!!" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Allocate space for data
			T *data = (T*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(T));
			if (data == NULL) {
				std::cout << "Not enough Memory!!!!" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Create input dataspace
			hsize_t dims_input[1];
			dims_input[0] = gridsize[0] * gridsize[1] * gridsize[2];
			hid_t input_sid = H5Screate_simple(1, dims_input, NULL);
			if (input_sid <= 0) std::cout << "HDF5 ERROR: Cannot create input dataspace!" << std::endl;

			// Read in the typing matrix
			status = readDataset("/LatTyp", TIME_STRING, input_fid, input_sid, H5T_NATIVE_INT, Type);
			if (status != 0) goto readFailed;

			// Open, read and close input dataset
			status = readDataset(VAR, TIME_STRING, input_fid, input_sid, H5Type, data);
			if (status != 0) goto readFailed;

			// Loop over grid sites
			for (int c = 0; c < gridsize[0] * gridsize[1] * gridsize[2]; c++) {

				// If not one of the include sites, ignore
				if (
					Type[c] == eRefined ||
					Type[c] == eRefinedInlet ||
					Type[c] == eRefinedSolid ||
					Type[c] == eRefinedSymmetry ||
					Type[c] == eTransitionToCoarser
					) continue;

					// Insert data into array
					vtkArray->InsertValue(count, data[c]);
					count++;
			}

			// Close input dataspace
			status = H5Sclose(input_sid);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot close input dataspace!" << std::endl;

			// Close input file
			status = H5Fclose(input_fid);
			if (status != 0) std::cout << "HDF5 ERROR: Cannot close input file!" << std::endl;

			// Free memory
			free(Type);
			free(data);

		}
	}

	// Add complete data set to grid
	grid->GetCellData()->AddArray(vtkArray);

	return 0;

readFailed:
	std::cout << "Read failed -- exiting early.";
	return 998;
}

/* End of header */