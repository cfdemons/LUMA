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

// Routines for reading and writing operations.

#include "../inc/stdafx.h"
#include <sstream>
#include "../inc/GridObj.h"
#include "../inc/MpiManager.h"
#include "../inc/ObjectManager.h"
#include "../inc/definitions.h"
#include "../inc/globalvars.h"
#include "../inc/hdf5luma.h"

// Static declarations
std::vector<H5File*> GridUtils::hdf_file;

using namespace std;

// *****************************************************************************
// Writes all the contents of the class at time t and call recursviely for any 
// subgrids. Writes to text file "Grids.out" by default.
void GridObj::io_textout(std::string output_tag) {

	// Get limits for current level
	size_t N_lim = XPos.size();
	size_t M_lim = YPos.size();
	size_t K_lim = ZPos.size();

	// Create stream and open text file
	ofstream gridoutput;
	gridoutput.precision(L_output_precision);

	// Construct File Name
	string FNameG, N_str, M_str, K_str, ex_str, L_NumLev_str, NumReg_str, mpirank;
	N_str = to_string((int)L_N);
	M_str = to_string((int)L_M);
	K_str = to_string((int)L_K);
	L_NumLev_str = to_string(level);
	if (L_NumLev == 0) ex_str = to_string(0);
	else ex_str = to_string(CoarseLimsX[0]) + string("_") + to_string(CoarseLimsY[0]) + string("_") + to_string(CoarseLimsZ[0]);
	if (L_NumLev == 0) NumReg_str = to_string(0);
	else NumReg_str = to_string(region_number);
	mpirank = to_string(MpiManager::my_rank);
	// Build string
	FNameG = string(GridUtils::path_str + "/Grids")
			+ string("D") +  to_string(L_dims)
			+ string("x") + N_str
			+ string("y") + M_str
			+ string("z") + K_str
			+ string("Lev") + L_NumLev_str
			+ string("Reg") + NumReg_str
			+ string("P") + ex_str
			+ string("Rnk") + mpirank
			+ string(".out");
	// Get character pointer
	const char* FNameG_c = FNameG.c_str();

	// If new simulation then overwrite if old file exists
	if (t == 0 && level == 0 && output_tag == "INITIALISATION") gridoutput.open(FNameG_c, ios::out);
	else gridoutput.open(FNameG_c, ios::out |ios::app);

	if ( gridoutput.is_open() ) {

		// Draw a line to begin
		gridoutput << "\n-------------------------------------------------------------------------------------" << endl;
		gridoutput << "-----------------------------------START OF OUTPUT-----------------------------------" << endl;
		gridoutput << "-------------------------------------------------------------------------------------" << endl;

		// Add tag
		gridoutput << output_tag << std::endl;

		// Print Grid Size header
		gridoutput << "L0 Grid Size = " << L_N << " x " << L_M << " x " << L_K << endl;
		gridoutput << "Local Grid Size = " << XPos.size() << " x " << YPos.size() << " x " << ZPos.size() << " (including any MPI overlap)" << std::endl;

		if (level == 0) {
			// If refined levels exist, print refinement ratio
			if (subGrid.size() != 0) {
				gridoutput << "Grid is refined." << endl;
				// Get size of regions
				for (size_t reg = 0; reg < subGrid.size(); reg++) {
					int finex = subGrid[reg].CoarseLimsX[1] - subGrid[reg].CoarseLimsX[0] + 1;
					int finey = subGrid[reg].CoarseLimsY[1] - subGrid[reg].CoarseLimsY[0] + 1;
					int finez = subGrid[reg].CoarseLimsZ[1] - subGrid[reg].CoarseLimsZ[0] + 1;
					gridoutput << "Local region # " << reg << " refinement = " << (((float)finex)*((float)finey)*((float)finez)*100) / (L_N*L_M*L_K) << "%" << endl;
				}
			}
			// Print time step
			string t_str = to_string(t);
			gridoutput << "Time Step = " << t << endl;
			gridoutput << "-------------------------------------------------------------------------------------" << endl;
		}

		// Print Grid Level
		string r_str = to_string(level);
		gridoutput << "Grid Level = " << r_str << endl;

		// Print region number
		string reg_str = to_string(region_number);
		gridoutput << "Region number = " << reg_str << endl;

		// Now print omega
		gridoutput << "Omega = " << omega << endl;

		// Index Vectors
		gridoutput << "X Index: ";
		for (size_t i = 0; i < N_lim; i++) {
			gridoutput << XInd[i] << "\t";
		}
		gridoutput << "\nY Index: ";
		for (size_t j = 0; j < M_lim; j++) {
			gridoutput << YInd[j] << "\t";
		}
		gridoutput << "\nZ Index: ";
		for (size_t k = 0; k < K_lim; k++) {
			gridoutput << ZInd[k] << "\t";
		}

		// Position Vectors
		gridoutput << "\nX Position: \t";
		for (size_t i = 0; i < N_lim; i++) {
			gridoutput << XPos[i] << "\t";
		}
		gridoutput << "\nY Position: \t";
		for (size_t j = 0; j < M_lim; j++) {
			gridoutput << YPos[j] << "\t";
		}
		gridoutput << "\nZ Position: \t";
		for (size_t k = 0; k < K_lim; k++) {
			gridoutput << ZPos[k] << "\t";
		}

		// Typing Matrix
		gridoutput << "\n\nTyping Matrix";
		for (size_t k = 0; k < K_lim; k++) {
			// New line with z-coordinate
			gridoutput << "\nz = " << ZPos[k] << "\n";

			for (size_t j = 0; j < M_lim; j++) {
				// New line
				gridoutput << "\n";
				for (size_t i = 0; i < N_lim; i++) {

					// Output
					gridoutput << LatTyp(i,j,k,M_lim,K_lim) << "\t";

				}
			}
		}

		// Populations (f, feq)
		gridoutput << "\n\nf Values";
		for (size_t v = 0; v < L_nVels; v++) {
			// Particular velocity
			gridoutput << "\nc = " << v+1;


			for (size_t k = 0; k < K_lim; k++) {
				// New line with z-coordinate
				gridoutput << "\nz = " << ZPos[k] << "\n";

				for (size_t j = 0; j < M_lim; j++) {
					// New line
					gridoutput << "\n";
					for (size_t i = 0; i < N_lim; i++) {

						// Output
						gridoutput << f(i,j,k,v,M_lim,K_lim,L_nVels) << "\t";

					}
				}
			}
		}

		gridoutput << "\n\nfeq Values";
		for (size_t v = 0; v < L_nVels; v++) {
			// Particular velocity
			gridoutput << "\nc = " << v+1;


			for (size_t k = 0; k < K_lim; k++) {
				// New line with z-coordinate
				gridoutput << "\nz = " << ZPos[k] << "\n";

				for (size_t j = 0; j < M_lim; j++) {
					// New line
					gridoutput << "\n";
					for (size_t i = 0; i < N_lim; i++) {

						// Output
						gridoutput << feq(i,j,k,v,M_lim,K_lim,L_nVels) << "\t";

					}
				}
			}
		}

		// Macroscopic (u, rho)
		gridoutput << "\n\nVelocity Values";
		for (size_t n = 0; n < L_dims; n++) {
			// Particular component
			gridoutput << "\nu(" << n+1 << ")";

			for (size_t k = 0; k < K_lim; k++) {
				// New line with z-coordinate
				gridoutput << "\nz = " << ZPos[k] << "\n";

				for (size_t j = 0; j < M_lim; j++) {
					// New line
					gridoutput << "\n";
					for (size_t i = 0; i < N_lim; i++) {

						// Output
						gridoutput << u(i,j,k,n,M_lim,K_lim,L_dims) << "\t";

					}
				}
			}
		}

		gridoutput << "\n\nDensity";
		for (size_t k = 0; k < K_lim; k++) {
			// New line with z-coordinate
			gridoutput << "\nz = " << ZPos[k] << "\n";

			for (size_t j = 0; j < M_lim; j++) {
				// New line
				gridoutput << "\n";
				for (size_t i = 0; i < N_lim; i++) {

					// Output
					gridoutput << rho(i,j,k,M_lim,K_lim) << "\t";

				}
			}
		}

		// Draw a line underneath
		gridoutput << "\n-------------------------------------------------------------------------------------" << endl;
		gridoutput << "------------------------------------END OF OUTPUT------------------------------------" << endl;
		gridoutput << "-------------------------------------------------------------------------------------" << endl;


		// Call recursively for all child subgrids
		size_t regions = subGrid.size();
		if (regions != 0) {
			for (size_t reg = 0; reg < regions; reg++) {

				subGrid[reg].io_textout(output_tag);

			}
		}


		// Close file
		gridoutput.close();

	} else {

		*GridUtils::logfile << "Cannot open file" << endl;

	}

}

// *****************************************************************************
// This routine writes/reads the current rank's data in the custom restart file 
// format to the file whose handle is provided.
void GridObj::io_restart(bool IO_flag) {

	if (IO_flag) {

		// Restart file stream
		std::ofstream file;
		*GridUtils::logfile << "Writing grid level " << level << " region " << region_number << " to restart file..." << endl;


		///////////////////////
		// LBM Data -- WRITE //
		///////////////////////

		if (MpiManager::my_rank == 0 && level == 0) { // Overwrite as first to write
			file.open(GridUtils::path_str + "/restart_LBM.out", std::ios::out);
		} else { // Append
			file.open(GridUtils::path_str + "/restart_LBM.out", std::ios::out | std::ios::app);
		}

		// Get grid sizes
		int N_lim, M_lim, K_lim;
		N_lim = XInd.size();
		M_lim = YInd.size();
		K_lim = ZInd.size();

		// Counters
		int i,j,k,v;

		// Write out global grid indices and then the values of f, u and rho
		for (k = 0; k < K_lim; k++) {
			for (j = 0; j < M_lim; j++) {
				for (i = 0; i < N_lim; i++) {
					
					// Don't write out the receiver layer sites to avoid duplication
					if (GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k])) continue;

					// Grid level and region
					file << level << "\t" << region_number << "\t";

					// Global Index
					file << XInd[i] << "\t" << YInd[j] << "\t" << ZInd[k] << "\t";

					// f values
					for (v = 0; v < L_nVels; v++) {
						file << f(i,j,k,v,YInd.size(),ZInd.size(),L_nVels) << "\t";
					}

					// u values
					for (v = 0; v < L_dims; v++) {
						file << u(i,j,k,v,YInd.size(),ZInd.size(),L_dims) << "\t";
					}

					// rho value
					file << rho(i,j,k,YInd.size(),ZInd.size()) << std::endl;

				}
			}
		}

		// Close file
		file.close();


		///////////////////////
		// IBM Data -- WRITE //
		///////////////////////

#ifdef L_IBM_ON

		ObjectManager::getInstance()->io_restart(IO_flag, level);

#endif



	} else {

		// Input stream
		std::ifstream file;
		*GridUtils::logfile << "Initialising grid level " << level << " region " << region_number << " from restart file..." << endl;


		//////////////////////
		// LBM Data -- READ //
		//////////////////////

		file.open("./restart_LBM.out", std::ios::in);
		if (!file.is_open()) {
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "Error opening LBM restart file. Exiting." << std::endl;
			exit(LUMA_FAILED);
		}
		// Counters, sizes and indices
		int i,j,k,v;
		int N_lim = XInd.size(), M_lim = YInd.size(), K_lim = ZInd.size();
		int gi, gj, gk;
		int in_level, in_regnum;
		std::vector<int> ind;

		// Read in one line of file at a time
		std::string line_in;	// String to store line
		std::istringstream iss;	// Buffer stream

		while( !file.eof() ) {

			// Get line and put in buffer
			std::getline(file,line_in,'\n');
			iss.str(line_in);
			iss.seekg(0); // Reset buffer position to start of buffer

			// Read in level and region
			iss >> in_level >> in_regnum;

			// Check level and region number
			if (in_level != level || in_regnum != region_number) continue;

			// Read in global indices
			iss >> gi >> gj >> gk;

			// Check on this rank before proceding
			if ( !GridUtils::isOnThisRank(gi,gj,gk,*this) ) continue;

			// Get local indices
			ind.clear();
			GridUtils::global_to_local(gi,gj,gk,this,ind);
			i = ind[0];
			j = ind[1];
#if (L_dims == 3)
			k = ind[2];
#else
			k = 0;
#endif


			// Read in f values
			for (v = 0; v < L_nVels; v++) {
				iss >> f(i,j,k,v,M_lim,K_lim,L_nVels);
			}

			// Read in u values
			for (v = 0; v < L_dims; v++) {
				iss >> u(i,j,k,v,M_lim,K_lim,L_dims);
			}

			// Read in rho value
			iss >> rho(i,j,k,M_lim,K_lim);

		}

		// Reached end of file so close file
		file.close();


		//////////////////////
		// IBM Data -- READ //
		//////////////////////

#ifdef L_IBM_ON

		ObjectManager::getInstance()->io_restart(IO_flag, level);

#endif


	}

	// Call recursively for subgrids
	if (level < L_NumLev && subGrid.size()) {
		for (size_t g = 0; g < subGrid.size(); g++) {
			subGrid[g].io_restart(IO_flag);
		}
	}


}
// *****************************************************************************
// Custom routine for writing out point probes or other high frequency, low 
// volume data
void GridObj::io_probeOutput() {

	// Declarations
	std::ofstream probefile;
	int i,j,d,i_local,j_local,k_local,i_global,j_global,k_global;
	int M_lims = YInd.size(), K_lims = ZInd.size();

	if (t == L_out_every_probe && MpiManager::my_rank == 0) {
		// Overwrite existing first time through
		probefile.open("./output/probe.out", std::ios::out);
	} else {
		// Append to existing
		probefile.open("./output/probe.out", std::ios::out | std::ios::app);
	}

	// Start a new line if first rank
	if (MpiManager::my_rank == 0) probefile << std::endl;


	// Probe spacing in each direction
	int pspace[L_dims];
	pspace[0] = abs(xProbeLims[1] - xProbeLims[0]) / (nProbes[0] - 1);
	pspace[1] = abs(yProbeLims[1] - yProbeLims[0]) / (nProbes[1] - 1);
#if (L_dims == 3)
	pspace[2] = abs(zProbeLims[1] - zProbeLims[0]) / (nProbes[2] - 1);
#endif

	// Loop over probe points
	for (i = 0; i < nProbes[0]; i++) {
		i_global = xProbeLims[0] + i*pspace[0];

		for (j = 0; j < nProbes[1]; j++) {
			j_global = yProbeLims[0] + j*pspace[1];

#if (L_dims == 3)
			for (int k = 0; k < nProbes[2]; k++) {
				k_global = zProbeLims[0] + k*pspace[2];
#else
			k_global = 0; {
#endif

				// DEBUG
				*GridUtils::logfile << "Writing probe point: " << 
					i_global << "," << j_global << "," << k_global << 
					std::endl;

				// Is point on rank, if not move on
				if (!GridUtils::isOnThisRank(i_global,j_global,k_global,*this)) {
					continue;
				}

				// Convert global to local if necessary
#ifdef L_BUILD_FOR_MPI
				i_local = i_global - XInd[1] + 1;
				j_local = j_global - YInd[1] + 1;
#if (L_dims == 3)
				k_local = k_global - ZInd[1] + 1;
#else
				k_local = k_global;
#endif

#else
				i_local = i_global; j_local = j_global; k_local = k_global;
#endif
				

				// Write out value and add tab
				for (d = 0; d < L_dims; d++) {
					probefile << u(i_local,j_local,k_local,d,M_lims,K_lims,L_dims) << "\t";
				}

			}

		}

	}

	probefile.close();


}

// ************************************************************** //
// Generic writer for each rank to write out all data row-wise to be 
// processed using a new post-processing application into a suitable 
// output format.
void GridObj::io_lite(double tval, std::string TAG) {

	std::ofstream litefile;

	// Filename
	std::string filename ("./" + GridUtils::path_str + "/io_lite.Lev" + std::to_string(level) + ".Reg" + std::to_string(region_number)
			+ ".Rnk" + std::to_string(MpiManager::my_rank) + "." + std::to_string((int)tval) + ".dat");

	// Create file
	litefile.open(filename, std::ios::out);

	// Set precision and force fixed formatting
	litefile.precision(L_output_precision);
	litefile.setf(std::ios::fixed);
	litefile.setf(std::ios::showpoint);

	// Write simple header
	litefile << "L" << level << " R" << region_number << " P" << std::to_string(MpiManager::my_rank) << " -- " << TAG << std::endl;
	litefile << "L_Timesteps = " << std::to_string(tval) << std::endl;
#ifdef L_MEGA_DEBUG
	litefile << "RANK TYPE X Y Z RHO UX UY UZ F FEQ TA_RHO TA_UX TA_UY TA_UZ TA_UXUX TA_UXUY TA_UXUZ TA_UYUY TA_UYUZ TA_UZUZ" << std::endl;
#else
	litefile << "RANK TYPE X Y Z RHO UX UY UZ TA_RHO TA_UX TA_UY TA_UZ TA_UXUX TA_UXUY TA_UXUZ TA_UYUY TA_UYUZ TA_UZUZ" << std::endl;
#endif
	
	// Indices
	size_t i,j,k,v;
		
	// Write out values
	for (k = 0; k < ZInd.size(); k++) {
		for (j = 0; j < YInd.size(); j++) {
			for (i = 0; i < XInd.size(); i++) {


#if (defined L_MEGA_DEBUG && !defined L_INC_RECV_LAYER && defined L_BUILD_FOR_MPI) || (!defined L_MEGA_DEBUG && defined L_BUILD_FOR_MPI)
				// Don't write out the receiver overlap in MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif				
				{

					// Write out rank
					litefile << MpiManager::my_rank << "\t";
				
					// Write out type
					litefile << LatTyp(i,j,k,YInd.size(),ZInd.size()) << "\t";

					// Write out X, Y, Z
					litefile << XPos[i] << "\t" << YPos[j] << "\t" << ZPos[k] << "\t";

					// Write out rho and u
					litefile << rho(i,j,k,YInd.size(),ZInd.size()) << "\t";
					for (v = 0; v < L_dims; v++) {
						litefile << u(i,j,k,v,YInd.size(),ZInd.size(),L_dims) << "\t";
					}
#if (L_dims != 3)
					litefile << 0 << "\t";
#endif

#ifdef L_MEGA_DEBUG
					// Write out F and Feq
					for (v = 0; v < L_nVels; v++) {
						litefile << f(i,j,k,v,YInd.size(),ZInd.size(),L_nVels) << "\t";
					}
					for (v = 0; v < L_nVels; v++) {
						litefile << feq(i,j,k,v,YInd.size(),ZInd.size(),L_nVels) << "\t";
					}
#endif
				
					// Write out time averaged rho and u
					litefile << rho_timeav(i,j,k,YInd.size(),ZInd.size()) << "\t";
					for (v = 0; v < L_dims; v++) {
						litefile << ui_timeav(i,j,k,v,YInd.size(),ZInd.size(),L_dims) << "\t";
					}
#if (L_dims != 3)
					litefile << 0 << "\t";
#endif

					// Write out time averaged u products
					litefile << uiuj_timeav(i,j,k,0,YInd.size(),ZInd.size(),(3*L_dims-3)) << "\t";
					litefile << uiuj_timeav(i,j,k,1,YInd.size(),ZInd.size(),(3*L_dims-3)) << "\t";
#if (L_dims == 3)
					litefile << uiuj_timeav(i,j,k,2,YInd.size(),ZInd.size(),(3*L_dims-3)) << "\t";
#else
					litefile << 0 << "\t";
#endif
#if (L_dims == 3)
					litefile << uiuj_timeav(i,j,k,3,YInd.size(),ZInd.size(),(3*L_dims-3)) << "\t";
					litefile << uiuj_timeav(i,j,k,4,YInd.size(),ZInd.size(),(3*L_dims-3)) << "\t";
					litefile << uiuj_timeav(i,j,k,5,YInd.size(),ZInd.size(),(3*L_dims-3)) << "\t";
#else
					litefile << uiuj_timeav(i,j,k,2,YInd.size(),ZInd.size(),(3*L_dims-3)) << "\t";
					litefile << 0 << "\t" << 0 << "\t";
#endif

					// New line
					litefile << std::endl;


				}

			}
		}
	}


#ifndef L_MEGA_DEBUG
	// Now do any sub-grids
	if (L_NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].io_lite(tval,"");
		}
	}
#endif

}

// ***************************************************************************//
// HDF5 writer which writes everything out as scalar arrays but is a fast, 
// portable format for easier post-processing
int GridObj::io_hdf5(double tval) {


	// Construct file name
	H5std_string FILE_NAME(GridUtils::path_str + "/hdf_R" + std::to_string(region_number) + "N" + std::to_string(level) + ".h5");
	const std::string time_string("/Time_" + std::to_string(static_cast<int>(tval)));

	// Get local grid sizes
	size_t N_lim = XPos.size();
	size_t M_lim = YPos.size();
	size_t K_lim = ZPos.size();

	// File ID
	int fid = level + region_number * L_NumLev;

	// Try block to detect exceptions raised by any of the calls inside it
	try	{

		// Turn off the auto-printing when failure occurs so that we can
		// handle the errors appropriately
		Exception::dontPrint();

		// Top level grid creates the file and the time step
		if (t == 0)	GridUtils::hdf_file.push_back(new H5File(FILE_NAME, H5F_ACC_TRUNC));
		GridUtils::hdf_file[fid]->createGroup(time_string);
		
		// Create dimension dataspace in file
		hsize_t dims_size[L_dims];
#if (L_dims == 3)
		dims_size[2] = K_lim; dims_size[1] = M_lim; dims_size[0] = N_lim;
#else
		dims_size[1] = M_lim; dims_size[0] = N_lim;
#endif
		DataSpace* dims_space = new DataSpace(L_dims, dims_size);

		// Create dataspace for positions (1D arrays)
		DataSpace* mem_pos_space = NULL;
		hsize_t dpos_size[1];

		// Create dataset handle
		DataSet* dataset = NULL;

		// Define hyperslab in memory space for velocity (flattened 3D arrays)
		hsize_t start[1];	// Start of hyperslab
		hsize_t stride[1];	// Stride of hyperslab
		hsize_t count[1];	// Block count
		hsize_t block[1];	// Block sizes
		start[0] = 0;
		stride[0] = L_dims;
		count[0] = N_lim * M_lim * K_lim;
		block[0] = 1;
		hsize_t d1_size[1];
		d1_size[0] = u.size();
		DataSpace* mem_space = new DataSpace(1, d1_size);
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

		// Define hyperslab in file for positions
		hsize_t start_p[L_dims];		// Start of hyperslab
		hsize_t stride_p[L_dims];		// Stride of hyperslab
		hsize_t count_p[L_dims];		// Block count
		hsize_t block_p[L_dims];		// Block sizes
		start_p[0] = 0; start_p[1] = 0;
		stride_p[0] = 1; stride_p[1] = 1;
		count_p[0] = 1; count_p[1] = 1;
		block_p[0] = 1; block_p[1] = 1;
#if (L_dims == 3)
		start_p[2] = 0; stride_p[2] = 1; count_p[2] = 1; block_p[2] = 1;
#endif

		// WRITE POSITION X
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/XPos", PredType::NATIVE_DOUBLE, *dims_space));
		dpos_size[0] = XPos.size();	mem_pos_space = new DataSpace(1, dpos_size);	// 1D array in memory
		count_p[0] = XPos.size();
		for (int j = 0; j < YPos.size(); j++) {
			start_p[1] = j;
#if (L_dims == 3)
			for (int k = 0; k < ZPos.size(); k++) {
				start_p[2] = k;
#else
			{
#endif
				dims_space->selectHyperslab(H5S_SELECT_SET, count_p, start_p, stride_p, block_p);	// Select slab in file
				dataset->write(&XPos[0], PredType::NATIVE_DOUBLE, *mem_pos_space, *dims_space);

			}
		}
		dataset->close();
		dims_space->selectNone();	// Cancel hyperslab selection
		count_p[0] = 1; start_p[1] = 0;	// Reset slab definitions for next one
#if (L_dims == 3)
		start_p[2] = 0;
#endif
		delete mem_pos_space;
		delete dataset;

		// WRITE POSITION Y
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/YPos", PredType::NATIVE_DOUBLE, *dims_space));
		dpos_size[0] = YPos.size();	mem_pos_space = new DataSpace(1, dpos_size);	// 1D array in memory
		count_p[1] = YPos.size();
		for (int i = 0; i < XPos.size(); i++) {
			start_p[0] = i;
#if (L_dims == 3)
			for (int k = 0; k < ZPos.size(); k++) {
				start_p[2] = k;
#else
				{
#endif
					dims_space->selectHyperslab(H5S_SELECT_SET, count_p, start_p, stride_p, block_p);	// Select slab in file
					dataset->write(&YPos[0], PredType::NATIVE_DOUBLE, *mem_pos_space, *dims_space);

				}
			}
		dataset->close();
		dims_space->selectNone();
		count_p[1] = 1; start_p[0] = 0;	// Reset slab definitions
#if (L_dims == 3)
		start_p[2] = 0;
#endif
		delete mem_pos_space;
		delete dataset;

		// WRITE POSITION Z
#if (L_dims)
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/ZPos", PredType::NATIVE_DOUBLE, *dims_space));
		dpos_size[0] = ZPos.size();	mem_pos_space = new DataSpace(1, dpos_size);	// 1D array in memory
		count_p[2] = ZPos.size();
		for (int i = 0; i < XPos.size(); i++) {
			start_p[0] = i;
			for (int j = 0; j < YPos.size(); j++) {
				start_p[1] = j;

				dims_space->selectHyperslab(H5S_SELECT_SET, count_p, start_p, stride_p, block_p);	// Select slab in file
				dataset->write(&ZPos[0], PredType::NATIVE_DOUBLE, *mem_pos_space, *dims_space);

				}
			}
		dataset->close();
		dims_space->selectNone();
		count_p[2] = 1; start_p[0] = 0; start_p[1] = 0;	// Reset slab definitions
		delete mem_pos_space;
		delete dataset;
#endif
		
		// WRITE LATTYP
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/LatTyp", PredType::NATIVE_INT, *dims_space));
		dataset->write(&LatTyp[0], PredType::NATIVE_INT);
		dataset->close();
		delete dataset;
		
		// WRITE RHO
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Rho", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&rho[0], PredType::NATIVE_DOUBLE);
		dataset->close();
		delete dataset;
	
		// WRITE UX
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Ux", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&u[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;

		// WRITE UY
		mem_space->selectNone();
		start[0] = 1;	// Shift by 1 element in memory
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Uy", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&u[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;

		// WRITE UZ
#if (L_dims == 3)
		mem_space->selectNone();
		start[0] = 2;	// Shift by 1 element in memory
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Uz", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&u[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;
#endif

		// WRITE RHO_TIMEAV
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Rho_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&rho_timeav[0], PredType::NATIVE_DOUBLE);
		dataset->close();
		delete dataset;

		// WRITE UX_TIMEAV
		mem_space->selectNone();
		start[0] = 0;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Ux_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&ui_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;

		// WRITE UY_TIMEAV
		mem_space->selectNone();
		start[0] = 1;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Uy_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&ui_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;

		// WRITE UZ_TIMEAV
#if (L_dims == 3)
		mem_space->selectNone();
		start[0] = 2;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/Uz_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&ui_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;
#endif

		// WRITE UXUX_TIMEAV
		delete mem_space;
		d1_size[0] = uiuj_timeav.size();	// New 1D mem_space defined for uiuj_timeav array
		mem_space = new DataSpace(1, d1_size);
		mem_space->selectNone();
		stride[0] = (3 * L_dims - 3);	// There are 3(d -1) velocity products
		start[0] = 0;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/UxUx_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;

		// WRITE UXUY_TIMEAV
		mem_space->selectNone();
		start[0] = 1;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/UxUy_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;


		// WRITE UYUY_TIMEAV
		mem_space->selectNone();
#if (L_dims == 3)
		start[0] = 3;
#else
		start[0] = 2;
#endif
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/UyUy_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;


#if (L_dims == 3)
		// WRITE UXUZ_TIMEAV
		mem_space->selectNone();
		start[0] = 2;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/UxUz_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;

		// WRITE UYUZ_TIMEAV
		mem_space->selectNone();
		start[0] = 4;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/UyUz_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;

		// WRITE UZUZ_TIMEAV
		mem_space->selectNone();
		start[0] = 5;
		mem_space->selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		dataset = new DataSet(
			GridUtils::hdf_file[fid]->createDataSet(time_string + "/UzUz_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
		dataset->close();
		delete dataset;
#endif

		// Clean up dataspaces
		delete dims_space;
		delete mem_space;


		// Call recursively on child grids
		if (L_NumLev > level) for (GridObj& g : subGrid) g.io_hdf5(tval);



	}  // End Try

	// Catch failure caused by the H5File operations
	catch (FileIException error)
	{
		error.printError();
		return -1;
	}

	// Catch failure caused by the DataSet operations
	catch (DataSetIException error)
	{
		error.printError();
		return -1;
	}

	// Catch failure caused by the DataSpace operations
	catch (DataSpaceIException error)
	{
		error.printError();
		return -1;
	}



	return 0;

}
// ***************************************************************************//
