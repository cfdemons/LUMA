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

#include "h5mgm.h"

/* H5 Multi-Grid Merge Tool for post-processing HDF5 files written by LUMA */
int main(int argc, char* argv[])
{
	// Print out to screen
	std::cout << "H5MultiGridMerge (h5mgm) Version " << H5MGM_VERSION << std::endl;
	std::cout << "Merging files..." << std::endl;

	// Turn auto error printing off
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	// Construct L0 filename
	std::string IN_FILE_NAME("./hdf_R0N0.h5");

	// T = 0 output file name string
	std::string OUT_FILE_NAME("./hdf_merged.0.h5");

	// Set T = 0 time string
	std::string TIME_STRING = "/Time_0";

	// Declarations
	herr_t status = 0;
	hid_t output_fid = NULL;
	hid_t input_fid = NULL;
	hid_t input_aid = NULL;
	int dimensions_p, levels, regions, timesteps, out_every;
	std::string VAR;
	int* LatTyp = nullptr;
	double *dummy_d = nullptr;
	int *dummy_i = nullptr;

	// Open L0 input file
	input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (input_fid == NULL) std::cout << "HDF5 ERROR: Cannot open input file!" << std::endl;

	// Read in key attributes from the file
	input_aid = H5Aopen(input_fid, "Dimensions", H5P_DEFAULT);
	if (input_aid == NULL) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &dimensions_p);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "NumberOfGrids", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &levels);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "NumberOfRegions", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &regions);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "Timesteps", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &timesteps);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "OutputFrequency", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &out_every);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	// Store basic data
	int *gridsize = (int*)malloc(3 * sizeof(int));	// Space for grid size
	gridsize[2] = 1;	// Set 3D dimension to 1, will get overwritten if actually 3D
	int num_grids = ((levels - 1) * regions + 1);	// Total number of grids

	// Close file
	status = H5Fclose(input_fid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close file!" << std::endl;
	
	// Time loop
	int count = 0;
	for (size_t t = 0; t <= (size_t)timesteps; t += out_every) {

		// Do data read, write one variable at a time (always do LatTyp first!)
		status = bufferData("/LatTyp", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_INT, LatTyp, dummy_i);
		if (status != 0) goto readFailed;
		status = bufferData("/Rho", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/Rho_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/Ux", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/Uy", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/Ux_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/Uy_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/UxUx_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/UxUy_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/UyUy_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/XPos", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;
		status = bufferData("/YPos", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
		if (status != 0) goto readFailed;

		if (dimensions_p == 3) {
			status = bufferData("/Uz", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
			if (status != 0) goto readFailed;
			status = bufferData("/Uz_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
			if (status != 0) goto readFailed;
			status = bufferData("/UxUz_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
			if (status != 0) goto readFailed;
			status = bufferData("/UyUz_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
			if (status != 0) goto readFailed;
			status = bufferData("/UzUz_TimeAv", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
			if (status != 0) goto readFailed;
			status = bufferData("/ZPos", TIME_STRING, OUT_FILE_NAME, gridsize, levels, regions, output_fid, H5T_NATIVE_DOUBLE, LatTyp, dummy_d);
			if (status != 0) goto readFailed;
		}

		// Close output file
		status = H5Fclose(output_fid);
		output_fid = NULL;
		if (status != 0) std::cout << "HDF5 ERROR: Cannot close output file!" << std::endl;

		// Update filename for next loop
		OUT_FILE_NAME = "./hdf_merged." + std::to_string(t + out_every) + ".h5";
		TIME_STRING = "/Time_" + std::to_string(t + out_every);

		// Print progress to screen
		std::cout << std::to_string((int)(((float)(t + out_every) / (float)timesteps) * 100.0f)) << "% complete." << std::endl;

	}

	return 0;

readFailed:
	std::cout << "Read failed -- exiting early.";
	return 999;
}

