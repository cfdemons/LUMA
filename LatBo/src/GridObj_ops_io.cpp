// Routines for reading and writing operations.

#include "../inc/stdafx.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "../inc/GridObj.h"
#include "../inc/MpiManager.h"
#include "../inc/ObjectManager.h"
#include "../inc/definitions.h"
#include "../inc/globalvars.h"


using namespace std;

// ***************************************************************************************************
// Writes all the contents of the class at time t and call recursviely for any subgrids.
// Writes to text file "Grids.out" by default.
void GridObj::io_textout(std::string output_tag) {

	// Get limits for current level
	size_t N_lim = XPos.size();
	size_t M_lim = YPos.size();
	size_t K_lim = ZPos.size();

	// Create stream and open text file
	ofstream gridoutput;
	gridoutput.precision(6);

	// Construct File Name
	string FNameG, N_str, M_str, K_str, ex_str, NumLev_str, NumReg_str, mpirank;
	N_str = to_string((int)N);
	M_str = to_string((int)M);
	K_str = to_string((int)K);
	NumLev_str = to_string(level);
	if (NumLev == 0) ex_str = to_string(0);
	else ex_str = to_string(CoarseLimsX[0]) + string("_") + to_string(CoarseLimsY[0]) + string("_") + to_string(CoarseLimsZ[0]);
	if (NumLev == 0) NumReg_str = to_string(0);
	else NumReg_str = to_string(region_number);
	mpirank = to_string(MpiManager::my_rank);
	// Build string
	FNameG = string(GridUtils::path_str + "/Grids")
			+ string("D") +  to_string(dims)
			+ string("x") + N_str
			+ string("y") + M_str
			+ string("z") + K_str
			+ string("Lev") + NumLev_str
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
		gridoutput << "L0 Grid Size = " << N << " x " << M << " x " << K << endl;
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
					gridoutput << "Local region # " << reg << " refinement = " << (((float)finex)*((float)finey)*((float)finez)*100) / (N*M*K) << "%" << endl;
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
		for (size_t v = 0; v < nVels; v++) {
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
						gridoutput << f(i,j,k,v,M_lim,K_lim,nVels) << "\t";

					}
				}
			}
		}

		gridoutput << "\n\nfeq Values";
		for (size_t v = 0; v < nVels; v++) {
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
						gridoutput << feq(i,j,k,v,M_lim,K_lim,nVels) << "\t";

					}
				}
			}
		}

		// Macroscopic (u, rho)
		gridoutput << "\n\nVelocity Values";
		for (size_t n = 0; n < dims; n++) {
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
						gridoutput << u(i,j,k,n,M_lim,K_lim,dims) << "\t";

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

// ***************************************************************************************************
// This routine writes/reads the current rank's data in the custom restart file format to the file
// whose handle is provided.
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
		unsigned int N_lim, M_lim, K_lim, min, minz;
		N_lim = XInd.size();
		M_lim = YInd.size();
		K_lim = ZInd.size();
		min = 0; minz = min;

		// If building for MPI then correct grid sizes to avoid writing out outer overlap
		if (level == 0) {

#ifdef BUILD_FOR_MPI
			min = 1; minz = min;
			N_lim = XInd.size()-1;
			M_lim = YInd.size()-1;
#if (dims == 3)
			K_lim = ZInd.size()-1;
#else
			K_lim = 1;
			minz = 0;
#endif
#endif

		}

		// Counters
		unsigned int i,j,k,v;

		// Write out global grid indices and then the values of f, u and rho
		for (k = minz; k < K_lim; k++) {
			for (j = min; j < M_lim; j++) {
				for (i = min; i < N_lim; i++) {

					// Grid level and region
					file << level << "\t" << region_number << "\t";

					// Global Index
					file << XInd[i] << "\t" << YInd[j] << "\t" << ZInd[k] << "\t";

					// f values
					for (v = 0; v < nVels; v++) {
						file << f(i,j,k,v,YInd.size(),ZInd.size(),nVels) << "\t";
					}

					// u values
					for (v = 0; v < dims; v++) {
						file << u(i,j,k,v,YInd.size(),ZInd.size(),dims) << "\t";
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

#ifdef IBM_ON

		ObjectManager::io_restart(IO_flag, level);

#endif



	} else {

		// Input stream
		std::ifstream file;
		*GridUtils::logfile << "Initialising grid level " << level << " region " << region_number << " from restart file..." << endl;


		//////////////////////
		// LBM Data -- READ //
		//////////////////////

		file.open(GridUtils::path_str + "/restart_LBM.out", std::ios::in);
		if (!file.is_open()) {
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "Error opening LBM restart file. Exiting." << std::endl;
			exit(EXIT_FAILURE);
		}
		// Counters, sizes and indices
		int i,j,k,v;
		unsigned int N_lim = XInd.size(), M_lim = YInd.size(), K_lim = ZInd.size();
		int gi, gj, gk;
		unsigned int in_level, in_regnum;

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

			// Global and local are the same for serial code or for lower grids
			i = gi; j = gj; k = gk;

			// For MPI cases, need to check on level 0 that site is on this rank before proceding
#ifdef BUILD_FOR_MPI
				if (level == 0) {

					// Check whether on overlap or core of this rank

					// Convert global to local indices
					if (gi == XInd[0]) { i = 0;	} else if (gi == XInd[XInd.size()-1]) { i = XInd.size()-1; } // Overlap
					else if (gi >= XInd[1] && gi <= XInd[XInd.size()-2]) { i = gi - XInd[1] + 1; } // Core
					else { continue; } // Not on the rank

					if (gj == YInd[0]) { j = 0;	} else if (gj == YInd[YInd.size()-1]) { j = YInd.size()-1; } // Overlap
					else if (gj >= YInd[1] && gj <= YInd[YInd.size()-2]) { j = gj - YInd[1] + 1; } // Core
					else { continue; } // Not on the rank
#if (dims == 3)
					if (gk == ZInd[0]) { k = 0;	} else if (gk == ZInd[ZInd.size()-1]) { k = ZInd.size()-1; } // Overlap
					else if (gk >= ZInd[1] && gk <= ZInd[ZInd.size()-2]) { k = gk - ZInd[1] + 1; } // Core
					else { continue; } // Not on the rank
#else
					k = gk; // In 2D no need to convert
#endif
				}
#endif


			// Read in f values
			for (v = 0; v < nVels; v++) {
				iss >> f(i,j,k,v,M_lim,K_lim,nVels);
			}

			// Read in u values
			for (v = 0; v < dims; v++) {
				iss >> u(i,j,k,v,M_lim,K_lim,dims);
			}

			// Read in rho value
			iss >> rho(i,j,k,M_lim,K_lim);

		}

		// Reached end of file so close file
		file.close();


		//////////////////////
		// IBM Data -- READ //
		//////////////////////

#ifdef IBM_ON

		ObjectManager::io_restart(IO_flag, level);

#endif


	}

	// Call recursively for subgrids
	if (level < NumLev && subGrid.size()) {
		for (unsigned int g = 0; g < subGrid.size(); g++) {
			subGrid[g].io_restart(IO_flag);
		}
	}


}
// ***************************************************************************************************
// Custom routine for writing out point probes or other high frequency, low volume data
void GridObj::io_probe_output() {

	// Declarations
	std::ofstream probefile;
	unsigned int i,j,d,i_local,j_local,k_local,i_global,j_global,k_global;
	unsigned int M_lims = YInd.size(), K_lims = ZInd.size();

	if (t == out_every_probe && MpiManager::my_rank == 0) {
		// Overwrite existing first time through
		probefile.open("./output/probe.out", std::ios::out);
	} else {
		// Append to existing
		probefile.open("./output/probe.out", std::ios::out | std::ios::app);
	}

	// Start a new line if first rank
	if (MpiManager::my_rank == 0) probefile << std::endl;


	// Probe spacing in each direction
	int pspace[dims];
	pspace[0] = abs(xProbeLims[1] - xProbeLims[0]) / (nProbes[0] - 1);
	pspace[1] = abs(yProbeLims[1] - yProbeLims[0]) / (nProbes[1] - 1);
#if (dims == 3)
	pspace[2] = abs(zProbeLims[1] - zProbeLims[0]) / (nProbes[2] - 1);
#endif

	// Loop over probe points
	for (i = 0; i < nProbes[0]; i++) {
		i_global = xProbeLims[0] + i*pspace[0];

		for (j = 0; j < nProbes[1]; j++) {
			j_global = yProbeLims[0] + j*pspace[1];

#if (dims == 3)
			for (unsigned int k = 0; k < nProbes[2]; k++) {
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
#ifdef BUILD_FOR_MPI
				i_local = i_global - XInd[1] + 1;
				j_local = j_global - YInd[1] + 1;
#if (dims == 3)
				k_local = k_global - ZInd[1] + 1;
#else
				k_local = k_global;
#endif

#else
				i_local = i_global; j_local = j_global; k_local = k_global;
#endif
				

				// Write out value and add tab
				for (d = 0; d < dims; d++) {
					probefile << u(i_local,j_local,k_local,d,M_lims,K_lims,dims) << "\t";
				}

			}

		}

	}

	probefile.close();


}
// ***************************************************************************************************
