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

#include "h5mgm.h"

template <typename T>
class VelocitySorter
{

public:

	VelocitySorter() {};
	~VelocitySorter() {};

	// Wrapper
	void readAndSort()
	{
		// Construct L0 filename
		std::string IN_FILE_NAME("./hdf_R0N0.h5");

		// Set T = 0 time string
		std::string TIME_STRING = "/Time_0";

		// Declarations
		hid_t output_fid = NULL;
		hid_t input_fid = NULL;
		hid_t input_aid = NULL;
		int dimensions_p, timesteps, out_every;
		double dx;
		int gridsize[3] = { 1, 1, 1 };

		// Open first file
		input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

		// Read in relevant attributes from the file
		input_aid = H5Aopen(input_fid, "Dimensions", H5P_DEFAULT);
		H5Aread(input_aid, H5T_NATIVE_INT, &dimensions_p);
		H5Aclose(input_aid);

		input_aid = H5Aopen(input_fid, "Timesteps", H5P_DEFAULT);
		H5Aread(input_aid, H5T_NATIVE_INT, &timesteps);
		H5Aclose(input_aid);

		input_aid = H5Aopen(input_fid, "OutputFrequency", H5P_DEFAULT);
		H5Aread(input_aid, H5T_NATIVE_INT, &out_every);
		H5Aclose(input_aid);

		input_aid = H5Aopen(input_fid, "Dx", H5P_DEFAULT);
		H5Aread(input_aid, H5T_NATIVE_DOUBLE, &dx);
		H5Aclose(input_aid);

		input_aid = H5Aopen(input_fid, "GridSize", H5P_DEFAULT);
		H5Aread(input_aid, H5T_NATIVE_INT, &gridsize[0]);
		H5Aclose(input_aid);

		// Resize data matrix
		totalSites = gridsize[0] * gridsize[1] * gridsize[2];
		if (!xyzuvw) delete xyzuvw;
		xyzuvw = new std::vector< std::vector<double> >();
		xyzuvw->resize(totalSites, std::vector<double>(6, 0.0));

		// Create input dataspace
		hsize_t dims_input[1];
		dims_input[0] = gridsize[0] * gridsize[1] * gridsize[2];
		hid_t input_sid = H5Screate_simple(1, dims_input, NULL);

		// Allocate buffers
		std::vector<double> *Xunsorted = new std::vector<double>(totalSites, 0.0);
		std::vector<double> *Yunsorted = new std::vector<double>(*Xunsorted);
		std::vector<double> *Zunsorted = new std::vector<double>(*Xunsorted);
		std::vector<double> *UXunsorted = new std::vector<double>(*Xunsorted);
		std::vector<double> *UYunsorted = new std::vector<double>(*Xunsorted);
		std::vector<double> *UZunsorted = new std::vector<double>(*Xunsorted);

		// Read in X, Y and Z
		readDataset("/XPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, &(*Xunsorted)[0]);
		readDataset("/YPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, &(*Yunsorted)[0]);
		if (dimensions_p == 3)
			readDataset("/ZPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, &(*Zunsorted)[0]);

		// Loop over timesteps
		clock_t startTime;
		double outerLoopTime = 0.0;
		for (size_t t = 0; t <= (size_t)timesteps; t += out_every)
		{
			// Start clock
			startTime = clock();

			// Create output filename
			std::string sortedFilename = H5MGM_OUTPUT_PATH;
			sortedFilename += "/sortedVelocity_" + std::to_string(t) + ".txt";

			// Read in velocity arrays for this timestep
			readDataset("/Ux", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, &(*UXunsorted)[0]);
			readDataset("/Uy", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, &(*UYunsorted)[0]);
			if (dimensions_p == 3)
				readDataset("/Uz", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, &(*UZunsorted)[0]);

			// Pack information into format required (what a pain!)
			for (int cell = 0; cell < totalSites; ++cell)
			{
				
				*((xyzuvw->begin() + cell)->begin() + 0) = *(Xunsorted->begin() + cell);
				*((xyzuvw->begin() + cell)->begin() + 1) = *(Yunsorted->begin() + cell);
				*((xyzuvw->begin() + cell)->begin() + 2) = *(Zunsorted->begin() + cell);
				*((xyzuvw->begin() + cell)->begin() + 3) = *(UXunsorted->begin() + cell);
				*((xyzuvw->begin() + cell)->begin() + 4) = *(UYunsorted->begin() + cell);
				*((xyzuvw->begin() + cell)->begin() + 5) = *(UZunsorted->begin() + cell);
			}

			// Sort data z then y then x
			sortrows(*xyzuvw, 2);
			sortrows(*xyzuvw, 1);
			sortrows(*xyzuvw, 0);

			// Write file
			std::ofstream outfile;
			outfile.open(sortedFilename, std::ios::out);
			for (std::vector<double> &row : *xyzuvw)
			{
				for (double &val : row)
				{
					outfile << val << "\t";
				}
				outfile << std::endl;
			}
			outfile.close();

			// Update time string to next timestep
			TIME_STRING = "/Time_" + std::to_string(t + out_every);

			// Print progress to screen
			outerLoopTime *= static_cast<double>(t);
			outerLoopTime += (static_cast<double>(clock() - startTime)) / CLOCKS_PER_SEC * 1000;
			outerLoopTime /= static_cast<double>(t + out_every);
			int hms[3];
			secs2hms((timesteps - (t + out_every)) * outerLoopTime / 1000, &hms[0]);
			std::cout << "\r" << std::to_string((int)(((float)(t + out_every) /
				(float)(timesteps + out_every)) * 100.0f)) << "% complete. Time to complete approx. " 
				<< hms[0] << " [h] " << hms[1] << " [m] " << hms[2] << " [s]      " << std::flush;

		}

		// Close input dataspace
		H5Sclose(input_sid);

		// Close input file
		H5Fclose(input_fid);

	};

private:
	// Data arranged as X Y Z U V W column vectors
	std::vector< std::vector<double> > *xyzuvw;

	// Others
	size_t totalSites;

	// Method to read a dataset with a given name into a buffer (taken from h5mgm)
	void readDataset(std::string VAR, std::string TIME_STRING, hid_t input_fid, hid_t input_sid, hid_t H5Type, T *buffer)
	{
		std::string variable_string = TIME_STRING + VAR;
		hid_t input_did = H5Dopen(input_fid, variable_string.c_str(), H5P_DEFAULT);
		H5Dread(input_did, H5Type, input_sid, H5S_ALL, H5P_DEFAULT, buffer);
		H5Dclose(input_did);
	};

	// Method to perform a vector sort favouring value in column col
	void sortrows(std::vector<std::vector<T>>& matrix, int col)
	{
		std::stable_sort(
			matrix.begin(),
			matrix.end(),
			[col](const std::vector<T>& lhs, const std::vector<T>& rhs)
		{
			return lhs[col] < rhs[col];
		});
	};

};