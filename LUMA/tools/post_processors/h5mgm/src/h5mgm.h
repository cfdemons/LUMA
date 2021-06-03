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

/* H5 Multi-Grid Merge Tool for post-processing HDF5 files written by LUMA */

#define H5MGM_VERSION "0.3.2"

#include "hdf5.h"
#define H5_BUILT_AS_DYNAMIC_LIB
#define HDF5_EXT_ZLIB
#define HDF5_EXT_SZIP

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

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
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkPointData.h"

#define DATASET_READ_FAIL -321
#define EARLY_EXIT 998
#define H5MGM_OUTPUT_PATH "./postprocessedoutput"

// Static variables
static bool bQuiet = false;
static bool bLoud = false;
static std::ofstream logfile;
static bool bCutSolid = false;
static bool bLegacy = false;
static bool bSorter = false;

// Unit vectors for node positions on each cell
const int e[3][8] =
{
	{ -1, -1, -1, -1,  1,  1,  1,  1 },
	{ -1, -1,  1,  1, -1, -1,  1,  1 },
	{ -1,  1, -1,  1, -1,  1, -1,  1 }
};

const int e2[2][4] =
{
	{ -1, -1,  1, 1 },
	{ -1,  1, -1, 1 }
};

// Typing enumeration from LUMA
enum eType
{
	eSolid,					///< Rigid, solid site with no-slip BC
	eFluid,					///< Fluid site
	eRefined,				///< Fluid site which is represented on a finer grid
	eTransitionToCoarser,	///< Fluid site coupled to a coarser grid
	eTransitionToFiner,		///< Fluid site coupled to a finer grid
	eBFL,					///< Site containing a BFL marker
	eVelocity,				///< Velocity boundary
	eExtrapolation,			///< Extrapolation (outlet) site
	eSlip					///< Slip boundary
};

// Error enumeration
enum eError
{
	eHDF,				// HDF5 IO error
	eFatal				// Fatal error (non-HDF5)
};

void writeInfo(std::string info, eError type) {

	// Create header label
	std::string label;
	switch (type)
	{
	case eHDF:
		label = "HDF ERROR: ";
		break;

	case eFatal:
		label = "FATAL ERROR: ";
		break;
	}

	// Perform appropriate action
	if (bQuiet == false)
	{
		// Write to log
		logfile << label << info << std::endl;
	}
	
	if (bLoud == true)
	{
		// Write to screen as well
		logfile << label << info << std::endl;
	}

}

bool isOnIgnoreList(eType cell_type)
{

	if (
		cell_type == eRefined ||
		cell_type == eTransitionToCoarser ||
		(bCutSolid = true && cell_type == eSolid)
		) return true;

	return false;
}

// Method to convert seconds to hours minutes and seconds
template <typename T>
static void secs2hms(T total_seconds, int *hms)
{
	// Round to the nearest second
	total_seconds = std::round(total_seconds);

	// Compute number of hours and minutes represented
	double hours = static_cast<double>(total_seconds) / (60.0 * 60.0);
	double minutes = static_cast<double>(total_seconds) / 60.0;

	// Compute components
	hms[0] = static_cast<int>(std::floor(hours));
	hms[1] = static_cast<int>(std::floor(minutes - (static_cast<double>(hms[0]) * 60.0)));
	hms[2] = static_cast<int>(std::floor(
		static_cast<double>(total_seconds)-
		static_cast<double>(hms[1]) * 60.0 -
		static_cast<double>(hms[0]) * 60.0 * 60.0
		));
}

// Method to read a dataset with a given name into a buffer
template <typename T>
herr_t readDataset(std::string VAR, std::string TIME_STRING, hid_t input_fid, hid_t input_sid, hid_t H5Type, T *buffer) {

	herr_t status;
	std::string variable_string = TIME_STRING + VAR;
	hid_t input_did = H5Dopen(input_fid, variable_string.c_str(), H5P_DEFAULT);
	if (input_did <= 0)
	{
		writeInfo("Cannot open input dataset: " + variable_string, eHDF);
		return DATASET_READ_FAIL;
	}
	status = H5Dread(input_did, H5Type, input_sid, H5S_ALL, H5P_DEFAULT, buffer);
	if (status != 0) writeInfo("Cannot read input dataset!", eHDF);
	status = H5Dclose(input_did);
	if (status != 0) writeInfo("Cannot close input dataset!", eHDF);

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
			if (input_fid <= 0) writeInfo("Cannot open input file!", eHDF);

			// Get local grid size
			hid_t input_aid = H5Aopen(input_fid, "GridSize", H5P_DEFAULT);
			if (input_aid <= 0) writeInfo("Cannot open attribute!", eHDF);
			herr_t status = H5Aread(input_aid, H5T_NATIVE_INT, gridsize);
			if (status != 0) writeInfo("Cannot read attribute!", eHDF);
			status = H5Aclose(input_aid);
			if (status != 0) writeInfo("Cannot close attribute!", eHDF);

			// Allocate space for typing matrix
			int *Type = (int*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(int));
			if (Type == NULL)
			{
				writeInfo("Not enough Memory!!!!", eFatal);
				exit(EXIT_FAILURE);
			}

			// Allocate space for data
			T *data = (T*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(T));
			if (data == NULL)
			{
				writeInfo("Not enough Memory!!!!", eFatal);
				exit(EXIT_FAILURE);
			}

			// Create input dataspace
			hsize_t dims_input[1];
			dims_input[0] = gridsize[0] * gridsize[1] * gridsize[2];
			hid_t input_sid = H5Screate_simple(1, dims_input, NULL);
			if (input_sid <= 0) writeInfo("Cannot create input dataspace!", eHDF);

			// Read in the typing matrix
			status = readDataset("/LatTyp", TIME_STRING, input_fid, input_sid, H5T_NATIVE_INT, Type);
			if (status == DATASET_READ_FAIL) return status;

			// Open, read and close input dataset
			status = readDataset(VAR, TIME_STRING, input_fid, input_sid, H5Type, data);
			if (status == DATASET_READ_FAIL) return status;

			// Loop over grid sites
			for (int c = 0; c < gridsize[0] * gridsize[1] * gridsize[2]; c++) {

				// If not one of the include sites, ignore
				if (isOnIgnoreList(static_cast<eType>(Type[c]))) continue;

				// Insert data into array
				vtkArray->InsertValue(count, data[c]);
				count++;
			}

			// Close input dataspace
			status = H5Sclose(input_sid);
			if (status != 0) writeInfo("Cannot close input dataspace!", eHDF);

			// Close input file
			status = H5Fclose(input_fid);
			if (status != 0) writeInfo("Cannot close input file!", eHDF);

			// Free memory
			free(Type);
			free(data);

		}
	}

	// Add complete data set to grid
	grid->GetCellData()->AddArray(vtkArray);

	return 0;
}

/* End of header */