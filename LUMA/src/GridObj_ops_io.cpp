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
		N_lim = static_cast<int>(XInd.size());
		M_lim = static_cast<int>(YInd.size());
		K_lim = static_cast<int>(ZInd.size());

		// Counters
		int i,j,k,v;

		// Write out global grid indices and then the values of f, u and rho
		for (k = 0; k < K_lim; k++) {
			for (j = 0; j < M_lim; j++) {
				for (i = 0; i < N_lim; i++) {
					
#ifdef L_BUILD_FOR_MPI
					// Don't write out the receiver layer sites to avoid duplication
					if (GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k])) continue;
#endif

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

		file.open("./input/restart_LBM.out", std::ios::in);
		if (!file.is_open()) {
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "Error opening LBM restart file. Exiting." << std::endl;
			exit(LUMA_FAILED);
		}
		// Counters, sizes and indices
		int i,j,k,v;
		int N_lim = static_cast<int>(XInd.size());
		int M_lim = static_cast<int>(YInd.size());
		int K_lim = static_cast<int>(ZInd.size());
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
	int M_lims = static_cast<int>(YInd.size()), K_lims = static_cast<int>(ZInd.size());

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
	litefile << "T = " << std::to_string(tval) << std::endl;
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
// HDF5 writer which writes quantities on grid out as scalar arrays
int GridObj::io_hdf5(double tval) {

#ifdef L_MPI_VERBOSE
	*MpiManager::logout << "Rank " << MpiManager::my_rank << ": Writing out Level " << level << ", Region " << region_number << std::endl;
#endif

	/***********************/
	/****** FILE SETUP *****/
	/***********************/

	// Construct filename
	std::string FILE_NAME(GridUtils::path_str + "/hdf_R" + std::to_string(region_number) + "N" + std::to_string(level) + ".h5");

	// Declarations
	hid_t file_id = NULL, plist_id = NULL, group_id = NULL;
	hid_t filespace = NULL; hsize_t dimsf[L_dims];
	hid_t memspace = NULL; hsize_t dimsm[1];
	hid_t attspace = NULL; hsize_t dimsa[1];
	hid_t dataset_id = NULL; hid_t attrib_id = NULL;
	herr_t status = 0;
	std::string variable_name;
	MpiManager::phdf5_struct p_data;

	// Turn auto error printing off
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	// Construct time string
	const std::string time_string("/Time_" + std::to_string(static_cast<int>(tval)));

	// Get local grid sizes
	int N_lim = static_cast<int>(XPos.size());
	int M_lim = static_cast<int>(YPos.size());
	int K_lim = static_cast<int>(ZPos.size());

	// Get modified local grid size (minus TL cells)
	int TL_thickness = 0;
	if (level != 0) TL_thickness = static_cast<int>(pow(2, level));
	int N_mod = N_lim - (2 * TL_thickness);
	int M_mod = M_lim - (2 * TL_thickness);
#if (L_dims == 3)
	int K_mod = K_lim - (2 * TL_thickness);
#else
	int K_mod = K_lim;
#endif


#ifdef L_BUILD_FOR_MPI

	///////////////////
	// PARALLEL CASE //
	///////////////////

	/* Cpp wrapper is not sufficient to access all the parallel IO function 
	 * and there is no tutorial on it so will have to implement in C for 
	 * parallel IO. Will leave the Cpp version of serial IO in though */


	// Retrieve writable data information to check whether this 
	// is a viable call to a parallel write.
	MpiManager* mpim = MpiManager::getInstance();
	for (MpiManager::phdf5_struct pd : mpim->p_data) {
		if (pd.level == level && pd.region == region_number) {
			p_data = pd;
			break;
		}
	}

	if (!p_data.writable_data_count) {
		*GridUtils::logfile << "Skipping HDF5 write as no writable data on this grid..." << std::endl;
		return -2;
	}

	// Create file parallel access property list
	MPI_Info info = MPI_INFO_NULL;
	plist_id = H5Pcreate(H5P_FILE_ACCESS);

	// Set communicator to be used
	if (level == 0)	{
		// Global communicator
		status = H5Pset_fapl_mpio(
			plist_id, MpiManager::getInstance()->world_comm, info);
	}

	// Else must be writable sub-grid so set communicator
	else {
		// Appropriate sub-grid communicator
		status = H5Pset_fapl_mpio(
			plist_id, MpiManager::getInstance()->subGrid_comm[(level - 1) + region_number * L_NumLev], info);
	}
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Set file access list failed: " << status << std::endl;

#else
	plist_id = H5P_DEFAULT;
#endif

	// Create/open file using the property list defined above
	if (t == 0) file_id = H5Fcreate(FILE_NAME.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	else file_id = H5Fopen(FILE_NAME.c_str(), H5F_ACC_RDWR, plist_id);
	if (file_id == NULL) *GridUtils::logfile << "HDF5 ERROR: Open file failed!" << std::endl;
	status = H5Pclose(plist_id);	 // Close access to property list now we have finished with it
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close file property list failed: " << status << std::endl;


	/***********************/
	/****** DATA SETUP *****/
	/***********************/

	// Create group
	group_id = H5Gcreate(file_id, time_string.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		
	// Compute dataspaces (file space then memory space)
	if (level == 0) {
		// L0 from definitions
		dimsf[0] = L_N;
		dimsf[1] = L_M;
#if (L_dims == 3)
		dimsf[2] = L_K;
#endif
	}
	else {
		// >L0 must get global sizes from the refinement region specification (ex. TL)
		dimsf[0] = (RefXend[level - 1][region_number] - RefXstart[level - 1][region_number] + 1) * 2 - (2 * TL_thickness);
		dimsf[1] = (RefYend[level - 1][region_number] - RefYstart[level - 1][region_number] + 1) * 2 - (2 * TL_thickness);
#if (L_dims == 3)
		dimsf[2] = (RefZend[level - 1][region_number] - RefZstart[level - 1][region_number] + 1) * 2 - (2 * TL_thickness);
#endif
	}
	filespace = H5Screate_simple(L_dims, dimsf, NULL);	// File space is globally sized

	// Memory space is always 1D scalar sized (ex. halo for MPI builds)
#ifdef L_BUILD_FOR_MPI
	dimsm[0] = p_data.writable_data_count;
#else
	dimsm[0] = N_mod * M_mod * K_mod;
#endif
	memspace = H5Screate_simple(1, dimsm, NULL);


	/***********************/
	/***** ATTRIBUTES ******/
	/***********************/

	if (t == 0) {

		// Create 1D attribute buffers
		int buffer_int_array[L_dims];
		int buffer_int = 0;
		buffer_int_array[0] = static_cast<int>(dimsf[0]);
		buffer_int_array[1] = static_cast<int>(dimsf[1]);
	#if (L_dims == 3)
		buffer_int_array[2] = static_cast<int>(dimsf[2]);
	#endif

		// Write Grid Size
		dimsa[0] = L_dims;
		attspace = H5Screate_simple(1, dimsa, NULL);
		attrib_id = H5Acreate(file_id, "GridSize", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int_array[0]);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
		status = H5Aclose(attrib_id);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

		if (level != 0) {
			// Write Sub-Grid Start
			buffer_int_array[0] = RefXstart[level - 1][region_number];
			buffer_int_array[1] = RefYstart[level - 1][region_number];
#if (L_dims == 3)
			buffer_int_array[2] = RefZstart[level - 1][region_number];
#endif
			attspace = H5Screate_simple(1, dimsa, NULL);
			attrib_id = H5Acreate(file_id, "RefinementStart", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int_array[0]);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

			// Write Sub-Grid End
			buffer_int_array[0] = RefXend[level - 1][region_number];
			buffer_int_array[1] = RefYend[level - 1][region_number];
#if (L_dims == 3)
			buffer_int_array[2] = RefZend[level - 1][region_number];
#endif
			attrib_id = H5Acreate(file_id, "RefinementEnd", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int_array[0]);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;
		}
		status = H5Sclose(attspace);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute space close failed: " << status << std::endl;
		
		// Write Timesteps
		buffer_int = L_Timesteps;
		dimsa[0] = 1;
		attspace = H5Screate_simple(1, dimsa, NULL);
		attrib_id = H5Acreate(file_id, "Timesteps", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
		status = H5Aclose(attrib_id);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

		// Write Out Frequency
		buffer_int = L_out_every;
		attrib_id = H5Acreate(file_id, "OutputFrequency", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
		status = H5Aclose(attrib_id);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

		// Write Levels
		buffer_int = L_NumLev + 1;
		attrib_id = H5Acreate(file_id, "NumberOfGrids", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
		status = H5Aclose(attrib_id);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

		// Write Regions
		buffer_int = L_NumReg;
		attrib_id = H5Acreate(file_id, "NumberOfRegions", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
		status = H5Aclose(attrib_id);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

		// Write Dimensions
		buffer_int = L_dims;
		attrib_id = H5Acreate(file_id, "Dimensions", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
		status = H5Aclose(attrib_id);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;
		status = H5Sclose(attspace);
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute space close failed: " << status << std::endl;

	}



	/***********************/
	/******* SCALARS *******/
	/***********************/

	// WRITE LATTYP
	variable_name = time_string + "/LatTyp";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eScalar, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &LatTyp[0], H5T_NATIVE_INT, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE RHO
	variable_name = time_string + "/Rho";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eScalar, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &rho[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE RHO_TIMEAV
	variable_name = time_string + "/Rho_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eScalar, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &rho_timeav[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;



	/***********************/
	/******* VECTORS *******/
	/***********************/

	// WRITE UX
	variable_name = time_string + "/Ux";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &u[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UY
	variable_name = time_string + "/Uy";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &u[1], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UZ
#if (L_dims == 3)
	variable_name = time_string + "/Uz";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod,  this, &u[2], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif

	// WRITE UX_TIMEAV
	variable_name = time_string + "/Ux_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &ui_timeav[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UY_TIMEAV
	variable_name = time_string + "/Uy_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &ui_timeav[1], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UZ_TIMEAV
#if (L_dims == 3)
	variable_name = time_string + "/Uz_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod,  this, &ui_timeav[2], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif


	/***********************/
	/*** PRODUCT VECTORS ***/
	/***********************/

	// WRITE UXUX_TIMEAV
	variable_name = time_string + "/UxUx_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &uiuj_timeav[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UXUY_TIMEAV
	variable_name = time_string + "/UxUy_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &uiuj_timeav[1], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UYUY_TIMEAV
	variable_name = time_string + "/UyUy_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#if (L_dims == 3)
	hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod,  this, &uiuj_timeav[3], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
#else
	hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &uiuj_timeav[2], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
#endif
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

#if (L_dims == 3)
	// WRITE UXUZ_TIMEAV
	variable_name = time_string + "/UxUz_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod,  this, &uiuj_timeav[2], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UYUZ_TIMEAV
	variable_name = time_string + "/UyUz_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod,  this, &uiuj_timeav[4], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE UZUZ_TIMEAV
	variable_name = time_string + "/UzUz_TimeAv";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod,  this, &uiuj_timeav[5], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif


	/***********************/
	/****** POSITIONS ******/
	/***********************/

	// WRITE POSITION X
	variable_name = time_string + "/XPos";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, ePosX, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &XPos[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

	// WRITE POSITION Y
	variable_name = time_string + "/YPos";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, ePosY, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod, this, &YPos[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;


	// WRITE POSITION Z
#if (L_dims == 3)
	variable_name = time_string + "/ZPos";
	dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hdf5_writeDataSet(memspace, filespace, dataset_id, ePosZ, N_lim, M_lim, K_lim, N_mod, M_mod, K_mod,  this, &ZPos[0], H5T_NATIVE_DOUBLE, TL_thickness, p_data);
	status = H5Dclose(dataset_id); // Close dataset
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif


	// Close memspace
	status = H5Sclose(memspace);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close memspace failed: " << status << std::endl;

	// Close filespace
	status = H5Sclose(filespace);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close filespace failed: " << status << std::endl;

	// Close group
	status = H5Gclose(group_id);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close group failed: " << status << std::endl;

	// Close file
	status = H5Fclose(file_id);
	if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close file failed: " << status << std::endl;



	// Call recursively on any present sub-grids
	if (level < L_NumLev) for (GridObj& g : subGrid) g.io_hdf5(tval);



	
/////////////////////////////// EVERYTHING BELOW HERE SHOULD BE SCRAPPED //////////////////////////////////////
	
//	/////////////////
//	// SERIAL CASE //
//	/////////////////
//
//	// Construct filename
//	H5std_string FILE_NAME(GridUtils::path_str + "/hdf_R" + std::to_string(region_number) + "N" + std::to_string(level) + ".h5");
//
//	// Try block_m to detect exceptions raised by any of the calls inside it
//	try	{
//
//		// Turn off the auto-printing when failure occurs so that we can
//		// handle the errors appropriately
//		Exception::dontPrint();
//		
//		// Create/open the file
//		H5File* hdf_file = NULL;
//		if (t == 0)	{
//
//			hdf_file = new H5File(FILE_NAME, H5F_ACC_TRUNC);	// New file, overwriting
//
//			// Create attribute data
//			int buffer_int_array[L_dims];
//			buffer_int_array[0] = N_lim;
//			buffer_int_array[1] = M_lim;
//#if (L_dims == 3)
//			buffer_int_array[2] = K_lim;
//#endif
//			int buffer_int = 0;
//
//			// Grid size
//			hsize_t attrib_size[1];
//			attrib_size[0] = L_dims;
//			DataSpace* d_attrib_space = new DataSpace(1, attrib_size);
//			Attribute *attrib = new Attribute(
//				hdf_file->createAttribute("GridSize", PredType::NATIVE_INT, *d_attrib_space));
//			attrib->write(PredType::NATIVE_INT, &buffer_int_array[0]);
//			attrib->close();
//			delete attrib;
//
//			// Only do on sub-grids
//			if (level != 0) {
//
//				// Grid extent start_m
//				buffer_int_array[0] = RefXstart[level - 1][region_number];
//				buffer_int_array[1] = RefYstart[level - 1][region_number];
//#if (L_dims == 3)
//				buffer_int_array[2] = RefZstart[level - 1][region_number];
//#endif
//				attrib = new Attribute(
//					hdf_file->createAttribute("RefinementStart", PredType::NATIVE_INT, *d_attrib_space));
//				attrib->write(PredType::NATIVE_INT, &buffer_int_array[0]);
//				attrib->close();
//				delete attrib;
//
//				// Grid extent end
//				buffer_int_array[0] = RefXend[level - 1][region_number];
//				buffer_int_array[1] = RefYend[level - 1][region_number];
//#if (L_dims == 3)
//				buffer_int_array[2] = RefZend[level - 1][region_number];
//#endif
//				attrib = new Attribute(
//					hdf_file->createAttribute("RefinementEnd", PredType::NATIVE_INT, *d_attrib_space));
//				attrib->write(PredType::NATIVE_INT, &buffer_int_array[0]);
//				attrib->close();
//				delete attrib;
//			}
//			d_attrib_space->close();
//			delete d_attrib_space;
//
//			// Time steps
//			buffer_int = L_Timesteps;
//			attrib_size[0] = 1;
//			d_attrib_space = new DataSpace(1, attrib_size);
//			attrib = new Attribute(
//				hdf_file->createAttribute("Timesteps", PredType::NATIVE_INT, *d_attrib_space));
//			attrib->write(PredType::NATIVE_INT, &buffer_int);
//			attrib->close();
//			delete attrib;
//
//			// Output frequency
//			buffer_int = L_out_every;
//			attrib = new Attribute(
//				hdf_file->createAttribute("OutputFrequency", PredType::NATIVE_INT, *d_attrib_space));
//			attrib->write(PredType::NATIVE_INT, &buffer_int);
//			attrib->close();
//			delete attrib;
//
//			// Number of levels
//			buffer_int = L_NumLev + 1;
//			attrib = new Attribute(
//				hdf_file->createAttribute("NumberOfLevels", PredType::NATIVE_INT, *d_attrib_space));
//			attrib->write(PredType::NATIVE_INT, &buffer_int);
//			attrib->close();
//			delete attrib;
//
//			// Number of regions
//			buffer_int = L_NumReg;
//			attrib = new Attribute(
//				hdf_file->createAttribute("NumberOfRegions", PredType::NATIVE_INT, *d_attrib_space));
//			attrib->write(PredType::NATIVE_INT, &buffer_int);
//			attrib->close();
//			delete attrib;
//
//			// Dimensions of the Problem
//			buffer_int = L_dims;
//			attrib = new Attribute(
//				hdf_file->createAttribute("Dimensions", PredType::NATIVE_INT, *d_attrib_space));
//			attrib->write(PredType::NATIVE_INT, &buffer_int);
//			attrib->close();
//			delete attrib;
//			d_attrib_space->close();
//			delete d_attrib_space;
//
//		}
//		else hdf_file = new H5File(FILE_NAME, H5F_ACC_RDWR);			// Existing file, open read-write
//
//		// Create time step group
//		hdf_file->createGroup(time_string);
//
//		// Create dataset handle
//		DataSet* dataset = NULL;
//				
//		// Create file space
//		hsize_t filespace_size[L_dims];
//		filespace_size[0] = N_mod;
//		filespace_size[1] = M_mod;
//#if (L_dims == 3)
//		filespace_size[2] = K_mod;
//#endif
//		DataSpace* file_space = new DataSpace(L_dims, filespace_size);
//
//
//		/***********************/
//		/****** POSITIONS ******/
//		/***********************/
//
//		// Initialise hyperslab in file for positions
//		hsize_t start_f[L_dims];		// Start of hyperslab
//		hsize_t stride_f[L_dims];		// Stride of hyperslab
//		hsize_t count_f[L_dims];		// Block count
//		hsize_t block_f[L_dims];		// Block sizes
//		start_f[0] = 0; start_f[1] = 0;
//		stride_f[0] = 1; stride_f[1] = 1;
//		count_f[0] = 1; count_f[1] = 1;
//		block_f[0] = 1; block_f[1] = 1;
//#if (L_dims == 3)
//		start_f[2] = 0;
//		stride_f[2] = 1;
//		count_f[2] = 1;
//		block_f[2] = 1;
//#endif
//		
//		// Create dataspace for positions (1D arrays)
//		DataSpace* mem_pos_space = NULL;
//		hsize_t dpos_size[1];
//
//		// WRITE POSITION X
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/XPos", PredType::NATIVE_DOUBLE, *file_space));
//		dpos_size[0] = N_mod;	mem_pos_space = new DataSpace(1, dpos_size);	// 1D array in memory covering the position vector
//		count_f[0] = N_mod;
//		for (int j = 0; j < YPos.size() - 2 * TL_thickness; j++) {
//			start_f[1] = j;
//#if (L_dims == 3)
//			for (int k = 0; k < ZPos.size() - 2 * TL_thickness; k++) {
//				start_f[2] = k;
//#else
//			{
//#endif
//				file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);	// Select slab in file
//				dataset->write(&XPos[TL_thickness], PredType::NATIVE_DOUBLE, *mem_pos_space, *file_space);
//
//			}
//		}
//		dataset->close();
//		file_space->selectNone();	// Cancel hyperslab selection
//		count_f[0] = 1; start_f[1] = 0;	// Reset slab definitions for next one
//#if (L_dims == 3)
//		start_f[2] = 0;
//#endif
//		delete mem_pos_space;
//		delete dataset;
//
//		// WRITE POSITION Y
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/YPos", PredType::NATIVE_DOUBLE, *file_space));
//		dpos_size[0] = M_mod;	mem_pos_space = new DataSpace(1, dpos_size);	// 1D array in memory
//		count_f[1] = M_mod;
//		for (int i = 0; i < XPos.size() - 2 * TL_thickness; i++) {
//			start_f[0] = i;
//#if (L_dims == 3)
//			for (int k = 0; k < ZPos.size() - 2 * TL_thickness; k++) {
//				start_f[2] = k;
//#else
//				{
//#endif
//					file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);	// Select slab in file
//					dataset->write(&YPos[TL_thickness], PredType::NATIVE_DOUBLE, *mem_pos_space, *file_space);
//
//				}
//			}
//		dataset->close();
//		file_space->selectNone();
//		count_f[1] = 1; start_f[0] = 0;	// Reset slab definitions
//#if (L_dims == 3)
//		start_f[2] = 0;
//#endif
//		delete mem_pos_space;
//		delete dataset;
//
//		// WRITE POSITION Z
//#if (L_dims == 3)
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/ZPos", PredType::NATIVE_DOUBLE, *file_space));
//		dpos_size[0] = K_mod;	mem_pos_space = new DataSpace(1, dpos_size);	// 1D array in memory
//		count_f[2] = K_mod;
//		for (int i = 0; i < XPos.size() - 2 * TL_thickness; i++) {
//			start_f[0] = i;
//			for (int j = 0; j < YPos.size() - 2 * TL_thickness; j++) {
//				start_f[1] = j;
//
//				file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);	// Select slab in file
//				dataset->write(&ZPos[TL_thickness], PredType::NATIVE_DOUBLE, *mem_pos_space, *file_space);
//
//				}
//			}
//		dataset->close();
//		file_space->selectNone();
//		count_f[2] = 1; start_f[0] = 0; start_f[1] = 0;	// Reset slab definitions
//		delete mem_pos_space;
//		delete dataset;
//#endif
//
//		/***********************/
//		/******* SCALARS *******/
//		/***********************/
//
//		// Define memory space for scalars
//		hsize_t mem_size[1];
//		mem_size[0] = N_lim * M_lim * K_lim; // Memory space includes TL
//		DataSpace* mem_space = new DataSpace(1, mem_size);
//
//		// Initialise hyperslab in memory space for scalars
//		hsize_t start_m[1];		// Start of hyperslab
//		hsize_t stride_m[1];	// Stride of hyperslab
//		hsize_t count_m[1];		// Block count
//		hsize_t block_m[1];		// Block sizes
//		// Memory offset due to not including TL
//#if (L_dims == 3)
//		start_m[0] = TL_thickness + TL_thickness * K_lim + TL_thickness * M_lim * K_lim;
//		block_m[0] = K_mod;
//		count_m[0] = M_mod;
//		stride_m[0] = block_m[0] + (2 * TL_thickness);
//#else
//		start_m[0] = TL_thickness + TL_thickness * M_lim;
//		block_m[0] = M_mod;
//		count_m[0] = N_mod;
//		stride_m[0] = block_m[0] + (2 * TL_thickness);
//#endif
//
//		// WRITE LATTYP
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/LatTyp", PredType::NATIVE_INT, *file_space));
//#if (L_dims == 3)
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = TL_thickness + TL_thickness * K_lim + i * M_lim * K_lim;
//#else
//		{
//#endif
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			dataset->write(&LatTyp[0], PredType::NATIVE_INT, *mem_space);
//			mem_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//		
//		// WRITE RHO
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/Rho", PredType::NATIVE_DOUBLE, *file_space));
//#if (L_dims == 3)
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = TL_thickness + TL_thickness * K_lim + i * M_lim * K_lim;
//#else
//		{
//#endif
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			dataset->write(&rho[0], PredType::NATIVE_DOUBLE, *mem_space);
//		}
//		dataset->close();
//		delete dataset;
//
//		// WRITE RHO_TIMEAV
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/Rho_TimeAv", PredType::NATIVE_DOUBLE, *file_space));
//#if (L_dims == 3)
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = TL_thickness + TL_thickness * K_lim + i * M_lim * K_lim;
//#else
//			{
//#endif
//				mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//				dataset->write(&rho_timeav[0], PredType::NATIVE_DOUBLE, *mem_space);
//			}
//			dataset->close();
//			delete dataset;
//
//		// Clear scalar memory space
//		mem_space->selectNone();
//		delete mem_space;
//
//
//		/***********************/
//		/******* VECTORS *******/
//		/***********************/
//
//		// Create memory space for D-dimensional vectors
//		mem_size[0] = N_lim * M_lim * K_lim * L_dims;
//		mem_space = new DataSpace(1, mem_size);
//
//		// Initialise hyperslab variables
//		start_f[0] = 0; start_f[1] = 0;
//		stride_f[0] = 1; stride_f[1] = 1;
//		count_f[0] = 1; count_f[1] = M_mod;
//		block_f[0] = 1; block_f[1] = 1;
//#if (L_dims == 3)
//		start_f[2] = 0; stride_f[2] = 1; count_f[2] = 1; block_f[2] = 1;
//#endif
//
//		// Initialise memory hyperslab variables for D-dimensional vectors
//#if (L_dims == 3)
//
//		// TODO
//#else
//		block_m[0] = 1;
//		count_m[0] = M_mod;
//		stride_m[0] = L_dims;
//#endif
//
//		// WRITE UX
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/Ux", PredType::NATIVE_DOUBLE, *file_space));
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			// Set memory offset
//			start_m[0] = 0 + TL_thickness * L_dims + i * L_dims * M_lim;
//			// Set file offset
//			start_f[0] = i - TL_thickness;
//			// Select memory hyperslab
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			// Select file hyperslab
//			file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);
//			// Write
//			dataset->write(&u[0], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//			// Clear selections
//			mem_space->selectNone();
//			file_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//		
//
//		// WRITE UY
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/Uy", PredType::NATIVE_DOUBLE, *file_space));
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = 0 + TL_thickness * L_dims + i * L_dims * M_lim;
//			start_f[0] = i - TL_thickness;
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);
//			dataset->write(&u[1], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//			mem_space->selectNone();
//			file_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//
//
//		// WRITE UZ
//#if (L_dims == 3)
//		// TODO
//#endif
//
//		// WRITE UX_TIMEAV
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/Ux_TimeAv", PredType::NATIVE_DOUBLE, *file_space));
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = 0 + TL_thickness * L_dims + i * L_dims * M_lim;
//			start_f[0] = i - TL_thickness;
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);
//			dataset->write(&ui_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//			mem_space->selectNone();
//			file_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//
//		// WRITE UY_TIMEAV
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/Uy_TimeAv", PredType::NATIVE_DOUBLE, *file_space));
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = 0 + TL_thickness * L_dims + i * L_dims * M_lim;
//			start_f[0] = i - TL_thickness;
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);
//			dataset->write(&ui_timeav[1], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//			mem_space->selectNone();
//			file_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//
//
//		// WRITE UZ_TIMEAV
//#if (L_dims == 3)
//		// TODO
//#endif
//
//		// Clear selections
//		mem_space->selectNone();
//		file_space->selectNone();
//		delete mem_space;
//
//
//		/***********************/
//		/****** PRODUCTS *******/
//		/***********************/
//
//		// Create memory space for D-dimensional vectors
//		mem_size[0] = N_lim * M_lim * K_lim * (3 * L_dims - 3);
//		mem_space = new DataSpace(1, mem_size);
//
//		// Set initial file hyperslab variables
//		start_f[0] = 0; start_f[1] = 0;
//		stride_f[0] = 1; stride_f[1] = 1;
//		count_f[0] = 1; count_f[1] = M_mod;
//		block_f[0] = 1; block_f[1] = 1;
//#if (L_dims == 3)
//		start_f[2] = 0; stride_f[2] = 1; count_f[2] = 1; block_f[2] = 1;
//#endif
//
//		// Set initial memory hyperslab variables for product vectors
//#if (L_dims == 3)
//
//		// TODO
//#else
//		block_m[0] = 1;
//		count_m[0] = M_mod;
//		stride_m[0] = 3 * L_dims - 3;
//#endif
//
//		// WRITE UXUX_TIMEAV
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/UxUx_TimeAv", PredType::NATIVE_DOUBLE, *file_space));
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = 0 + TL_thickness * (3 * L_dims - 3) + i * (3 * L_dims - 3) * M_lim;
//			start_f[0] = i - TL_thickness;
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);
//			dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//			mem_space->selectNone();
//			file_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//
//		// WRITE UXUY_TIMEAV
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/UxUy_TimeAv", PredType::NATIVE_DOUBLE, *file_space));
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = 0 + TL_thickness * (3 * L_dims - 3) + i * (3 * L_dims - 3) * M_lim;
//			start_f[0] = i - TL_thickness;
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);
//			dataset->write(&uiuj_timeav[1], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//			mem_space->selectNone();
//			file_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//
//
//		// WRITE UYUY_TIMEAV
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/UyUy_TimeAv", PredType::NATIVE_DOUBLE, *file_space));
//		for (int i = TL_thickness; i < N_lim - TL_thickness; i++) {
//			start_m[0] = 0 + TL_thickness * (3 * L_dims - 3) + i * (3 * L_dims - 3) * M_lim;
//			start_f[0] = i - TL_thickness;
//			mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//			file_space->selectHyperslab(H5S_SELECT_SET, count_f, start_f, stride_f, block_f);
//#if (L_dims == 3)
//			dataset->write(&uiuj_timeav[3], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//#else
//			dataset->write(&uiuj_timeav[2], PredType::NATIVE_DOUBLE, *mem_space, *file_space);
//#endif
//			
//			mem_space->selectNone();
//			file_space->selectNone();
//		}
//		dataset->close();
//		delete dataset;
//
//
//#if (L_dims == 3)
//		// WRITE UXUZ_TIMEAV
//		mem_space->selectNone();
//		start_m[0] = 2;
//		mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/UxUz_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
//		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
//		dataset->close();
//		delete dataset;
//
//		// WRITE UYUZ_TIMEAV
//		mem_space->selectNone();
//		start_m[0] = 4;
//		mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/UyUz_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
//		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
//		dataset->close();
//		delete dataset;
//
//		// WRITE UZUZ_TIMEAV
//		mem_space->selectNone();
//		start_m[0] = 5;
//		mem_space->selectHyperslab(H5S_SELECT_SET, count_m, start_m, stride_m, block_m);
//		dataset = new DataSet(
//			hdf_file->createDataSet(time_string + "/UzUz_TimeAv", PredType::NATIVE_DOUBLE, *dims_space));
//		dataset->write(&uiuj_timeav[0], PredType::NATIVE_DOUBLE, *mem_space, DataSpace::ALL);
//		dataset->close();
//		delete dataset;
//#endif
//
//		// Clean up dataspaces
//		delete file_space;
//		delete mem_space;
//
//
//		// Call recursively on child grids
//		if (L_NumLev > level) for (GridObj& g : subGrid) g.io_hdf5(tval);
//
//		// Close file
//		hdf_file->close();
//
//
//	}  // End Try
//
//	// Catch failure caused by the H5File operations
//	catch (FileIException error)
//	{
//		*GridUtils::logfile << "File I Exception: ";
//		error.printError();
//		return -51;
//	}
//
//	// Catch failure caused by the DataSet operations
//	catch (DataSetIException error)
//	{
//		*GridUtils::logfile << "Data Set I Exception: ";
//		error.printError();
//		return -52;
//	}
//
//	// Catch failure caused by the DataSpace operations
//	catch (DataSpaceIException error)
//	{
//		*GridUtils::logfile << "Data Space I Exception: ";
//		error.printError();
//		return -53;
//	}


	return 0;

}
// ***************************************************************************//

