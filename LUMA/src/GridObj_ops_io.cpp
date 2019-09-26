/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

// Routines for reading and writing operations.

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/ObjectManager.h"
#include "../inc/hdf5luma.h"

using namespace std;

// *****************************************************************************
/// \brief	Verbose ASCII writer.
///
///			Writes all the contents of the grid class at time t and call recursviely
///			for any sub-grids. Writes to text file "Grids.out" by default.
///
/// \param output_tag	text string added to top of output for identification.
void GridObj::io_textout(std::string output_tag) {

	int rank = GridUtils::safeGetRank();

	// Create stream and open text file
	ofstream gridoutput;
	gridoutput.precision(L_OUTPUT_PRECISION);

	// Construct File Name
	string FNameG, N_str, M_str, K_str, ex_str, L_NumLev_str, NumReg_str, mpirank;
	N_str = to_string((int)L_N);
	M_str = to_string((int)L_M);
	K_str = to_string((int)L_K);
	L_NumLev_str = to_string(level);
	if (L_NUM_LEVELS == 0) ex_str = to_string(0);
	else ex_str = to_string(CoarseLimsX[eMinimum]) + string("_") + to_string(CoarseLimsY[eMinimum]) + string("_") + to_string(CoarseLimsZ[eMinimum]);
	if (L_NUM_LEVELS == 0) NumReg_str = to_string(0);
	else NumReg_str = to_string(region_number);
	mpirank = to_string(rank);
	// Build string
	FNameG = string(GridUtils::path_str + "/Grids")
			+ string("D") +  to_string(L_DIMS)
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
		gridoutput << "Local Grid Size = " << N_lim << " x " << M_lim << " x " << K_lim << " (including any MPI overlap)" << std::endl;

		if (level == 0) {
			// If refined levels exist, print refinement ratio
			if (subGrid.size() != 0) {
				gridoutput << "Grid is refined." << endl;
				// Get size of regions
				for (size_t reg = 0; reg < subGrid.size(); reg++) {
					int finex = subGrid[reg]->CoarseLimsX[eMaximum] - subGrid[reg]->CoarseLimsX[eMinimum] + 1;
					int finey = subGrid[reg]->CoarseLimsY[eMaximum] - subGrid[reg]->CoarseLimsY[eMinimum] + 1;
					int finez = subGrid[reg]->CoarseLimsZ[eMaximum] - subGrid[reg]->CoarseLimsZ[eMinimum] + 1;
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
		for (size_t v = 0; v < L_NUM_VELS; v++) {
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
						gridoutput << f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) << "\t";

					}
				}
			}
		}

		gridoutput << "\n\nfeq Values";
		for (size_t v = 0; v < L_NUM_VELS; v++) {
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
						gridoutput << feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) << "\t";

					}
				}
			}
		}

		// Macroscopic (u, rho)
		gridoutput << "\n\nVelocity Values";
		for (size_t n = 0; n < L_DIMS; n++) {
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
						gridoutput << u(i,j,k,n,M_lim,K_lim,L_DIMS) << "\t";

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

				subGrid[reg]->io_textout(output_tag);

			}
		}


		// Close file
		gridoutput.close();

	} else {

		*GridUtils::logfile << "Cannot open file" << endl;

	}

}

// *****************************************************************************
/// \brief	.fga file writer.
///
///			Writes the components of the macroscopic velocity of the grid at time t and call
///         recursively for any sub-grid. Writes the data of each subgrid in a different .fga file. 
///			.fga is the ASCII file format used by Unreal Engine 4 to read the data that populates a 
///         VectorField object.
///         It doesn't do anything if the model is not 2D or 3D. Since .fga files can only store 3D data
///
void GridObj::io_fgaout() {
	
	// Write the file only if the model has the correct number of dimensions
	if ((L_DIMS == 2) || (L_DIMS == 3)){
		_io_fgaout(t);
	}
	else{
		cout << "Error: See Log File" << endl;
		*GridUtils::logfile << ".fga output can only be set for 2 or 3 dimensional models" << endl;
	}

}



// *****************************************************************************
/// \brief	.fga file writer.
///
///			Writes the components of the macroscopic velocity of the grid at time t and call
///         recursively for any sub-grid. Writes the data of each subgrid in a different .fga file. 
///			.fga is the ASCII file format used by Unreal Engine 4 to read the data that populates a 
///         VectorField object
///
void GridObj::_io_fgaout(int timeStepL0) {

	int rank = GridUtils::safeGetRank();
	
	// Create stream and open text file
	ofstream gridoutput;
	gridoutput.precision(L_OUTPUT_PRECISION);

	// Construct File Name
	string FNameG, N_str, M_str, K_str, ex_str, L_NumLev_str, NumReg_str, mpirank,L_timeStep_str;
	N_str = to_string((int)L_N);
	M_str = to_string((int)L_M);
	K_str = to_string((int)L_K);
	L_NumLev_str = to_string(level);
	L_timeStep_str = to_string(timeStepL0);
	//if (L_NUM_LEVELS == 0) ex_str = to_string(0);
	//else ex_str = to_string(CoarseLimsX[eMinimum]) + string("_") + to_string(CoarseLimsY[eMinimum]) + string("_") + to_string(CoarseLimsZ[eMinimum]);
	if (L_NUM_LEVELS == 0) NumReg_str = to_string(0);
	else NumReg_str = to_string(region_number);
	mpirank = to_string(rank);
	// Build string
	FNameG = string(GridUtils::path_str + "/Grids")
		+ string("D") + to_string(L_DIMS)
		+ string("x") + N_str
		+ string("y") + M_str
		+ string("z") + K_str
		+ string("dtNum") + L_timeStep_str
		+ string("Lev") + L_NumLev_str
		+ string("Reg") + NumReg_str
		//+ string("P") + ex_str
		+ string("Rnk") + mpirank
		+ string(".fga");
	// Get character pointer
	const char* FNameG_c = FNameG.c_str();

	gridoutput.open(FNameG_c, ios::out);

	if (gridoutput.is_open()) {

		// Print number of points (cells) in each direction
		gridoutput << N_lim << "," << M_lim << "," << K_lim << "," << endl;

		// Write the minimum coordinate of the grid. Convert the data to cm, since this is the distance unit used by UE4. 
		gridoutput << GridUnits::m2cm(dh / 2.0) << "," << GridUnits::m2cm(dh / 2.0) << "," << GridUnits::m2cm(dh / 2.0) << "," << endl;

		// Write the maximum coordinate of the grid
		gridoutput << GridUnits::m2cm(dh / 2.0 + (N_lim - 1) * dh) << "," << GridUnits::m2cm(dh / 2.0 + (M_lim - 1) * dh) << "," << GridUnits::m2cm(dh / 2.0 + (K_lim - 1) * dh) << "," << endl;

		// Auxiliary array to store the velocity data. This way I can use the same code with L_dim = 2 and L_dim = 3
		double v[3] = { 0.0, 0.0, 0.0 };

		for (size_t k = K_lim; k-- > 0;){

			for (size_t j = 0; j < M_lim; j++){

				for (size_t i = 0; i < N_lim; i++){

					for (size_t n = 0; n < L_DIMS; n++){
						// Fill the v array
						v[n] = u(i, j, k, n, M_lim, K_lim, L_DIMS);
					}
					// Write the data to the file. 
					gridoutput << GridUnits::m2cm(GridUnits::ulat2uphys(v[0], this)) << "," 
						       << GridUnits::m2cm(GridUnits::ulat2uphys(v[1], this)) << "," 
							   << GridUnits::m2cm(GridUnits::ulat2uphys(-v[2], this)) << "," 
							   << endl;
				}
			}
		}

		// Close file
		gridoutput.close();

		// Call recursively for all child subgrids
		size_t regions = subGrid.size();
		if (regions != 0) {
			for (size_t reg = 0; reg < regions; reg++) {

				subGrid[reg]->_io_fgaout(timeStepL0);
			}
		}

	}
	else {

		*GridUtils::logfile << "Cannot open file for .fga output" << endl;

	}

}

// *****************************************************************************
/// \brief	Restart file read-writer.
///
///			This routine writes/reads the current rank's particle distribution functions
///			in the custom restart file format. If the file already exists, data is appended. 
///			IB body data are also written out but no other body information at present. 
///			It writes the restart data following the procedure in: 
///			"Physically based Animation of Free Surface Flows with the Lattice Boltzmann Method" by N. Thuerey, 
///			section 6.1. 
///			This implementation allows to restart the simulation with a different time step value. 
///         The data written for each cell is: 
///			- Global position 
///         - Dimensionless velocity
///         - Density in LBM units
///         - Time-scaled non equilibrium distribution functions: ((f - f_eq) * omega) / (f_eq*dt)
///
/// \param IO_flag	flag to indicate whether a write or read
void GridObj::io_restart(eIOFlag IO_flag) {

	// Get GM Instance
	GridManager *gm = GridManager::getInstance();

	// Rank string
	std::string rnk_str = std::to_string(GridUtils::safeGetRank());


	if (IO_flag == eWrite) {

		// Restart file stream
		std::ofstream file;
		*GridUtils::logfile << "Writing grid level " << level << " region " << region_number << " to restart file..." << endl;


		///////////////////////
		// LBM Data -- WRITE //
		///////////////////////

		if (level == 0) {
			// New file
			file.open(GridUtils::path_str + "/restart_LBM_Rnk" + rnk_str + ".out", std::ios::out);
		}
		else {
			// Append
			file.open(GridUtils::path_str + "/restart_LBM_Rnk" + rnk_str + ".out", std::ios::out | std::ios::app);
		}

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

					// Global Position
					file << XPos[i] << "\t" << YPos[j] << "\t" << ZPos[k] << "\t";

					// Dimensionless u values
					for (v = 0; v < L_DIMS; v++) {
						file << GridUnits::ulbm2ud(u(i, j, k, v, M_lim, K_lim, L_DIMS), this) << "\t";
					}

					// rho in lbm units (rho does not depend on dt). 
					file << rho(i, j, k, M_lim, K_lim) << "\t";

					int id = k + j * K_lim + i * K_lim * M_lim;
					// time - scaled fneq values
					for (v = 0; v < L_NUM_VELS; v++) {
						double f_eq = _LBM_equilibrium_opt(id, v);
						double f_neq_restart = ((f(i, j, k, v, M_lim, K_lim, L_NUM_VELS) - f_eq) * omega) / (f_eq*dt);
						file << f_neq_restart << "\t";
					}

					file << std::endl;

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

		// Call write recursively for subgrids
		if (level < L_NUM_LEVELS && subGrid.size()) {
			for (size_t g = 0; g < subGrid.size(); g++) {
				subGrid[g]->io_restart(IO_flag);
			}
		}


	} else {

		// Input stream
		std::ifstream file;
		L_INFO("Initialising grids from restart file...", GridUtils::logfile);


		//////////////////////
		// LBM Data -- READ //
		//////////////////////

		file.open("./input/restart_LBM_Rnk" + rnk_str + ".out", std::ios::in);
		if (!file.is_open())
			L_ERROR("Error opening LBM restart file. Exiting.", GridUtils::logfile);

		// Counters, sizes and indices
		int i,j,k,v;
		double x, y, z;
		int in_level, in_regnum;
		std::vector<int> ind;

		// Read in one line of file at a time
		std::string line_in;	// String to store line
		std::istringstream iss;	// Buffer stream

		while(!file.eof())
		{

			// Get line and put in buffer
			std::getline(file,line_in,'\n');
			iss.str(line_in);
			iss.seekg(0); // Reset buffer position to start of buffer

			// Read in level and region
			iss >> in_level >> in_regnum;

			// Get grid
			GridObj *g = nullptr;
			GridUtils::getGrid(gm->Grids, in_level, in_regnum, g);
			if (!g) L_ERROR("Could not retrieve grid level " + std::to_string(in_level) + 
				" region " + std::to_string(in_regnum) + ". Is restart file correct?", GridUtils::logfile);

			// Read in positions
			iss >> x >> y >> z;

			// Get indices
			eLocationOnRank *loc = nullptr;
			std::vector<int> ijk;
			if (!GridUtils::isOnThisRank(x, y, z, loc, g, &ijk)) continue;
			i = ijk[0];
			j = ijk[1];
			k = ijk[2];
			int id = k + j * K_lim + i * K_lim * M_lim;

			// Read in u values and convert them to lbm units
			for (v = 0; v < L_DIMS; v++) {
				double u_temp;
				iss >> u_temp;
				g->u(i, j, k, v, g->M_lim, g->K_lim, L_DIMS) = GridUnits::ud2ulbm(u_temp, g);
			}

			// Read in rho value (it is already written in lbm units, since rho doesn't depend on dt)
			iss >> g->rho(i, j, k, g->M_lim, g->K_lim);

			// Read in f values and convert them to the new dt
			for (v = 0; v < L_NUM_VELS; v++) {
				double f_temp;
				double f_eq = _LBM_equilibrium_opt(id, v);
				iss >> f_temp;
				g->f(i, j, k, v, g->M_lim, g->K_lim, L_NUM_VELS) = f_eq*(1 + (g->dt*f_temp) / omega);
				g->fNew(i, j, k, v, g->M_lim, g->K_lim, L_NUM_VELS) = g->f(i, j, k, v, g->M_lim, g->K_lim, L_NUM_VELS);
			}

		}

		// Reached end of file so close file
		file.close();

		L_INFO("Restart complete.", GridUtils::logfile);


		//////////////////////
		// IBM Data -- READ //
		//////////////////////

#ifdef L_IBM_ON
			ObjectManager::getInstance()->io_restart(IO_flag, level);
#endif

	}

}

// *****************************************************************************
/// \brief	Probe writer.
///
///			This routine writes the quantities at the probe locations to a single 
///			file.
void GridObj::io_probeOutput() {

	int rank = GridUtils::safeGetRank();

	// Declarations
	std::ofstream probefile;
	int i, j, d;
	double x, y, z;

	if (t == 0 && rank == 0) {
		// Overwrite existing first time through
		probefile.open(GridUtils::path_str + "/probe.out", std::ios::out);
	} else {
		// Append to existing
		probefile.open(GridUtils::path_str + "/probe.out", std::ios::out | std::ios::app);
	}
	probefile.precision(L_OUTPUT_PRECISION);

	// Start a new line if first rank
	if (rank == 0 && t != 0) probefile << std::endl;

	// Declarations
	eLocationOnRank loc = eNone;
	GridObj *g = nullptr;
	std::vector<int> ijk;

	// Probe spacing in each direction
	double pspace[L_DIMS];
	if (cNumProbes[0] > 1)
		pspace[0] = abs(cProbeLimsX[1] - cProbeLimsX[0]) / (cNumProbes[0] - 1);
	if (cNumProbes[1] > 1)
		pspace[1] = abs(cProbeLimsY[1] - cProbeLimsY[0]) / (cNumProbes[1] - 1);
#if (L_DIMS == 3)
	if (cNumProbes[2] > 1)
		pspace[2] = abs(cProbeLimsZ[1] - cProbeLimsZ[0]) / (cNumProbes[2] - 1);
#endif

	// Loop over probe points to compute positions
	for (i = 0; i < cNumProbes[0]; i++) {
		x = cProbeLimsX[0] + i*pspace[0];

		for (j = 0; j < cNumProbes[1]; j++) {
			y = cProbeLimsY[0] + j*pspace[1];

#if (L_DIMS == 3)
			for (int k = 0; k < cNumProbes[2]; k++) {
				z = cProbeLimsZ[0] + k*pspace[2];
#else
			z = 0.0; {
#endif
				// Set written flag
				bool bProbeWritten = false;

				// Get each grid available on this rank in reverse order
				for (int lev = L_NUM_LEVELS; lev >= 0; --lev)
				{
					for (int reg = 0; reg < L_NUM_REGIONS; ++reg)
					{

						// Get grid if available
						GridUtils::getGrid(GridManager::getInstance()->Grids, lev, reg, g);
						if (g == nullptr) continue;

						// Determine if point is on this grid
						if (!GridUtils::isOnThisRank(x, y, z, &loc, g, &ijk)) continue;

						// Don't want it if on a halo as it will be duplicated
						if (loc == eHalo) continue;

						// As long as not on a TL to finer we can use it
						if (g->LatTyp(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim) == eTransitionToFiner) continue;

						// Write out velocity components and add tab
						for (d = 0; d < L_DIMS; d++)
						{
							probefile << g->u(ijk[0], ijk[1], ijk[2], d, g->M_lim, g->K_lim, L_DIMS) << "\t";
						}
#if (L_DIMS != 3)
						probefile << 0 << "\t";
#endif

						// Write out density
						probefile << g->rho(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim) << "\t";

						bProbeWritten = true;
						break;

					}

					if (bProbeWritten) break;
				}

			}
		}
	}

	probefile.close();


}

// *****************************************************************************
/// \brief	ASCII dump of grid data.
///
///			Generic ASCII writer for each rank to write out all grid data in rows 
///			into a single, unsorted file.
///
/// \param tval	time value being written out.
/// \param TAG	text identifier for the data.
void GridObj::io_lite(double tval, std::string TAG) {

	int rank = GridUtils::safeGetRank();
	std::ofstream litefile;

	// Filename
	std::string filename ("./" + GridUtils::path_str + "/io_lite.Lev" + std::to_string(level) + ".Reg" + std::to_string(region_number)
			+ ".Rnk" + std::to_string(rank) + "." + std::to_string((int)tval) + ".dat");

	// Create file
	litefile.open(filename, std::ios::out);

	// Set precision and force fixed formatting
	litefile.precision(L_OUTPUT_PRECISION);
	litefile.setf(std::ios::fixed);
	litefile.setf(std::ios::showpoint);
	
	// Indices
	size_t i,j,k,v;
		
	// Write out values
	for (k = 0; k < K_lim; k++) {
		for (j = 0; j < M_lim; j++) {
			for (i = 0; i < N_lim; i++) {


#if (defined L_INC_RECV_LAYER && defined L_BUILD_FOR_MPI)
				// Don't write out the receiver overlap in MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif				
				{

					// Write out rank
					litefile << rank << "\t";
				
					// Write out type
					litefile << LatTyp(i,j,k,M_lim,K_lim) << "\t";

					// Write out X, Y, Z
					litefile << XPos[i] << "\t" << YPos[j] << "\t" << ZPos[k] << "\t";

					// Write out rho and u
					litefile << rho(i,j,k,M_lim,K_lim) << "\t";
					for (v = 0; v < L_DIMS; v++) {
						litefile << u(i,j,k,v,M_lim,K_lim,L_DIMS) << "\t";
					}
#if (L_DIMS != 3)
					litefile << std::to_string(0.0) << "\t";
#endif

					// Write out F and Feq
					for (v = 0; v < L_NUM_VELS; v++) {
						litefile << f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) << "\t";
					}
					for (v = 0; v < L_NUM_VELS; v++) {
						litefile << fNew(i,j,k,v,M_lim,K_lim,L_NUM_VELS) << "\t";
					}
				
#ifdef L_COMPUTE_TIME_AVERAGED_QUANTITIES
					// Write out time averaged rho and u
					litefile << rho_timeav(i,j,k,M_lim,K_lim) << "\t";
					for (v = 0; v < L_DIMS; v++) {
						litefile << ui_timeav(i,j,k,v,M_lim,K_lim,L_DIMS) << "\t";
					}
#if (L_DIMS != 3)
					litefile << std::to_string(0.0) << "\t";
#endif

					// Write out time averaged u products
					litefile << uiuj_timeav(i,j,k,0,M_lim,K_lim,(3*L_DIMS-3)) << "\t";
					litefile << uiuj_timeav(i,j,k,1,M_lim,K_lim,(3*L_DIMS-3)) << "\t";
#if (L_DIMS == 3)
					litefile << uiuj_timeav(i,j,k,2,M_lim,K_lim,(3*L_DIMS-3)) << "\t";
#else
					litefile << std::to_string(0.0) << "\t";
#endif
#if (L_DIMS == 3)
					litefile << uiuj_timeav(i,j,k,3,M_lim,K_lim,(3*L_DIMS-3)) << "\t";
					litefile << uiuj_timeav(i,j,k,4,M_lim,K_lim,(3*L_DIMS-3)) << "\t";
					litefile << uiuj_timeav(i,j,k,5,M_lim,K_lim,(3*L_DIMS-3)) << "\t";
#else
					litefile << uiuj_timeav(i,j,k,2,M_lim,K_lim,(3*L_DIMS-3)) << "\t";
					litefile << std::to_string(0.0) << "\t" << std::to_string(0.0) << "\t";
#endif

#endif // L_COMPUTE_TIME_AVERAGED_QUANTITIES

					// New line
					litefile << std::endl;


				}

			}
		}
	}

	// Now do any sub-grids
	if (L_NUM_LEVELS > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg]->io_lite(tval, "");
		}
	}

}

// *****************************************************************************
/// \brief	HDF5 writer.
///
///			Useful grid quantities written out as scalar arrays. Creates one 
///			*.h5 file per grid and data is grouped into timesteps within each 
///			file. Should be used with the merge tool at post-processing to 
///			conver to sructured VTK output readable in paraview.
///
/// \param tval	time value being written out.
int GridObj::io_hdf5(double tval)
{

	// Get GM and lower edge information
	GridManager *gm = GridManager::getInstance();
	double minEdges[3];
	if (level == 0)
		for (int d = 0; d < 3; ++d) minEdges[d] = 0.0;
	else
	{
		minEdges[eXDirection] = gm->global_edges[eXMin][level + region_number * L_NUM_LEVELS];
		minEdges[eYDirection] = gm->global_edges[eYMin][level + region_number * L_NUM_LEVELS];
		minEdges[eZDirection] = gm->global_edges[eZMin][level + region_number * L_NUM_LEVELS];
	}

	// Get MPIM
#ifdef L_BUILD_FOR_MPI
	MpiManager *mpim = MpiManager::getInstance();

#ifdef L_MPI_VERBOSE
	*mpim->logout << "Writing out Level " << std::to_string(level) << ", Region " << std::to_string(region_number) << std::endl;
#endif

#endif

	/***********************/
	/****** FILE SETUP *****/
	/***********************/

	// Construct filename
	std::string FILE_NAME(GridUtils::path_str + 
		"/hdf_R" + std::to_string(region_number) + 
		"N" + std::to_string(level) + ".h5");

	// ID declarations
	hid_t file_id = static_cast<hid_t>(NULL);
	hid_t plist_id = static_cast<hid_t>(NULL);
	hid_t group_id = static_cast<hid_t>(NULL);
	hid_t filespace = static_cast<hid_t>(NULL);
	hid_t memspace = static_cast<hid_t>(NULL);
	hid_t attspace = static_cast<hid_t>(NULL);
	hid_t dataset_id = static_cast<hid_t>(NULL);
	hid_t attrib_id = static_cast<hid_t>(NULL);

	// Dimensions of file, memory and attribute spaces
	hsize_t dimsf[L_DIMS];
	hsize_t dimsm[1];
	hsize_t dimsa[1];

	// Others
	herr_t status = 0;
	std::string variable_name;
	HDFstruct p_data;
	int TL_thickness;
	bool TL_present[3];		// Access using eCartesianDirection

	// Turn auto error printing off
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	// Construct time string
	const std::string time_string("/Time_" + std::to_string(static_cast<int>(tval)));

	// Set TL thickness
	if (level == 0) {
		// TL is zero on L0 grid
		TL_thickness = 0;
		TL_present[eXDirection] = false;
		TL_present[eYDirection] = false;
		TL_present[eZDirection] = false;
	}
	else {
		// TL thickness is 2 cells on sub-grids
		TL_thickness = 2;
		TL_present[eXDirection] = gm->subgrid_tlayer_key[eXMin][level + region_number * L_NUM_LEVELS - 1];
		TL_present[eYDirection] = gm->subgrid_tlayer_key[eYMin][level + region_number * L_NUM_LEVELS - 1];
		TL_present[eZDirection] = gm->subgrid_tlayer_key[eZMin][level + region_number * L_NUM_LEVELS - 1];
	}



	// Retrieve writable data information for this grid from GM
	for (HDFstruct pd : gm->p_data) {
		if (pd.level == level && pd.region == region_number) {
			p_data = pd;
			break;
		}
	}

#ifdef L_BUILD_FOR_MPI

	///////////////////
	// PARALLEL CASE //
	///////////////////

	if (p_data.writable_data_count)
	{

		// Create file parallel access property list
		MPI_Info info = MPI_INFO_NULL;
		plist_id = H5Pcreate(H5P_FILE_ACCESS);

		// Set communicator to be used
		if (level == 0)	{
			// Global communicator
			status = H5Pset_fapl_mpio(
				plist_id, mpim->world_comm, info
				);
		}
		else {
			// Appropriate sub-grid communicator
			status = H5Pset_fapl_mpio(
				plist_id, mpim->subGrid_comm[(level - 1) + region_number * L_NUM_LEVELS], info
				);
		}
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Set file access list failed: " << status << std::endl;

#else

		// Simple serial property list
		plist_id = H5P_DEFAULT;

#endif // L_BUILD_FOR_MPI

		// Create/open file using the property list defined above
		if (t == 0) file_id = H5Fcreate(FILE_NAME.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
		else file_id = H5Fopen(FILE_NAME.c_str(), H5F_ACC_RDWR, plist_id);
		if (file_id == static_cast<hid_t>(NULL)) *GridUtils::logfile << "HDF5 ERROR: Open file failed!" << std::endl;
		status = H5Pclose(plist_id);	 // Close access to property list now we have finished with it
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close file property list failed: " << status << std::endl;

#ifdef L_BUILD_FOR_MPI
		// Synchronise after opening
		if (level == 0) MPI_Barrier(mpim->world_comm);
		else MPI_Barrier(mpim->subGrid_comm[(level - 1) + region_number * L_NUM_LEVELS]);
#endif


		/***********************/
		/****** DATA SETUP *****/
		/***********************/

		// Create group
		group_id = H5Gcreate(file_id, time_string.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		// Compute dataspaces (file space data in GM and ex. TL where appropriate)
		int idx = level + region_number * L_NUM_LEVELS;
		dimsf[0] = gm->global_size[eXDirection][idx];
		dimsf[1] = gm->global_size[eYDirection][idx];

		if (level > 0)
		{

#ifdef L_HDF_DEBUG
			// Write file space sizing parameters
			L_WARN("Global Size (inc. TL): Level " + std::to_string(level) + ", Region " + std::to_string(region_number) + ": "
				+ std::to_string(gm->global_size[eXDirection][idx]) + ","
				+ std::to_string(gm->global_size[eYDirection][idx]) + ","
				+ std::to_string(gm->global_size[eZDirection][idx]), GridUtils::logfile);
			L_INFO("TL Thickness: " + std::to_string(TL_thickness), GridUtils::logfile);
			L_INFO("TL keys: " + std::to_string(gm->subgrid_tlayer_key[eXMin][idx - 1]) + "," + std::to_string(gm->subgrid_tlayer_key[eXMax][idx - 1]) + "/" +
				std::to_string(gm->subgrid_tlayer_key[eYMin][idx - 1]) + "," + std::to_string(gm->subgrid_tlayer_key[eYMax][idx - 1]) + "/" +
				std::to_string(gm->subgrid_tlayer_key[eZMin][idx - 1]) + "," + std::to_string(gm->subgrid_tlayer_key[eZMax][idx - 1]),
				GridUtils::logfile);
#endif
			// Correct global size using knowledge of TLs that need clipping
			dimsf[0] -= 
				(gm->subgrid_tlayer_key[eXMin][idx - 1] * TL_thickness
				+ gm->subgrid_tlayer_key[eXMax][idx - 1] * TL_thickness);
			dimsf[1] -=
				(gm->subgrid_tlayer_key[eYMin][idx - 1] * TL_thickness
				+ gm->subgrid_tlayer_key[eYMax][idx - 1] * TL_thickness);
		}
#if (L_DIMS == 3)
		dimsf[2] = gm->global_size[eZDirection][idx];
		if (level > 0)
		{
			dimsf[2] -=
				(gm->subgrid_tlayer_key[eZMin][idx - 1] * TL_thickness
				+ gm->subgrid_tlayer_key[eZMax][idx - 1] * TL_thickness);
		}
#endif

#if (defined L_BUILD_FOR_MPI && defined L_HDF_DEBUG)
		// Check all members of communicator have same file size
		if (level != 0)	hdf_checkFileSpace(&dimsf[0], mpim->subGrid_comm[(level - 1) + region_number * L_NUM_LEVELS]);
#endif

		// File space is globally sized
		filespace = H5Screate_simple(L_DIMS, dimsf, NULL);

		// Write out file space to log file for reference
#ifdef L_HDF_DEBUG
		L_INFO("Level " + std::to_string(level) + ", Region " + std::to_string(region_number)
			+ ": Filespace size = "
			+ std::to_string(dimsf[eXDirection]) + " x "
			+ std::to_string(dimsf[eYDirection]) + " x "
#if (L_DIMS == 3)
			+ std::to_string(dimsf[eZDirection])
#else
			+ std::to_string(1)
#endif			
			, GridUtils::logfile);
#endif

		// Memory space is always 1D scalar sized (ex. TL and halo for MPI builds)
		dimsm[0] = p_data.writable_data_count;
		memspace = H5Screate_simple(1, dimsm, NULL);


		/***********************/
		/***** ATTRIBUTES ******/
		/***********************/

		if (t == 0)
		{

			// Create 1D attribute buffers
			int buffer_int_array[L_DIMS];
			int buffer_int = 0;
			double buffer_double = 0.0;
			buffer_int_array[0] = static_cast<int>(dimsf[0]);
			buffer_int_array[1] = static_cast<int>(dimsf[1]);
#if (L_DIMS == 3)
			buffer_int_array[2] = static_cast<int>(dimsf[2]);
#endif

			// Write Grid Size
			dimsa[0] = L_DIMS;
			attspace = H5Screate_simple(1, dimsa, NULL);
			attrib_id = H5Acreate(file_id, "GridSize", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int_array[0]);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;
			status = H5Sclose(attspace);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute space close failed: " << status << std::endl;

			// Write Timesteps
			buffer_int = L_TOTAL_TIMESTEPS;
			dimsa[0] = 1;
			attspace = H5Screate_simple(1, dimsa, NULL);
			attrib_id = H5Acreate(file_id, "Timesteps", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

			// Write Out Frequency
			buffer_int = L_GRID_OUT_FREQ;
			attrib_id = H5Acreate(file_id, "OutputFrequency", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

			// Write dh
			buffer_double = dh;
			attrib_id = H5Acreate(file_id, "Dx", H5T_NATIVE_DOUBLE, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_DOUBLE, &buffer_double);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

			// Write Levels
			buffer_int = L_NUM_LEVELS + 1;
			attrib_id = H5Acreate(file_id, "NumberOfGrids", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

			// Write Regions
			buffer_int = L_NUM_REGIONS;
			attrib_id = H5Acreate(file_id, "NumberOfRegions", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

			// Write MPI flag
#ifdef L_BUILD_FOR_MPI
			buffer_int = 1;
#else
			buffer_int = 0;
#endif
			attrib_id = H5Acreate(file_id, "Mpi", H5T_NATIVE_INT, attspace, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attrib_id, H5T_NATIVE_INT, &buffer_int);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute write failed: " << status << std::endl;
			status = H5Aclose(attrib_id);
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Attribute close failed: " << status << std::endl;

			// Write Dimensions
			buffer_int = L_DIMS;
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
		hdf5_writeDataSet(memspace, filespace, dataset_id, eScalar, this, &LatTyp[0], H5T_NATIVE_INT, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE RHO
		variable_name = time_string + "/Rho";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eScalar, this, &rho[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

#ifdef L_COMPUTE_TIME_AVERAGED_QUANTITIES

		// WRITE RHO_TIMEAV
		variable_name = time_string + "/Rho_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eScalar, this, &rho_timeav[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

#endif // L_COMPUTE_TIME_AVERAGED_QUANTITIES



		/***********************/
		/******* VECTORS *******/
		/***********************/

		// WRITE UX
		variable_name = time_string + "/Ux";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, this, &u[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UY
		variable_name = time_string + "/Uy";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, this, &u[1], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UZ
#if (L_DIMS == 3)
		variable_name = time_string + "/Uz";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, this, &u[2], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif

#ifdef L_COMPUTE_TIME_AVERAGED_QUANTITIES

		// WRITE UX_TIMEAV
		variable_name = time_string + "/Ux_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, this, &ui_timeav[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UY_TIMEAV
		variable_name = time_string + "/Uy_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, this, &ui_timeav[1], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UZ_TIMEAV
#if (L_DIMS == 3)
		variable_name = time_string + "/Uz_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eVector, this, &ui_timeav[2], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif


		/***********************/
		/*** PRODUCT VECTORS ***/
		/***********************/

		// WRITE UXUX_TIMEAV
		variable_name = time_string + "/UxUx_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, this, &uiuj_timeav[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UXUY_TIMEAV
		variable_name = time_string + "/UxUy_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, this, &uiuj_timeav[1], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UYUY_TIMEAV
		variable_name = time_string + "/UyUy_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#if (L_DIMS == 3)
		hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, this, &uiuj_timeav[3], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
#else
		hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, this, &uiuj_timeav[2], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
#endif
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

#if (L_DIMS == 3)
		// WRITE UXUZ_TIMEAV
		variable_name = time_string + "/UxUz_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, this, &uiuj_timeav[2], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UYUZ_TIMEAV
		variable_name = time_string + "/UyUz_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, this, &uiuj_timeav[4], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

		// WRITE UZUZ_TIMEAV
		variable_name = time_string + "/UzUz_TimeAv";
		dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hdf5_writeDataSet(memspace, filespace, dataset_id, eProductVector, this, &uiuj_timeav[5], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
		status = H5Dclose(dataset_id); // Close dataset
		if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif

#endif // L_COMPUTE_TIME_AVERAGED_QUANTITIES

		// Only write positions and block labels on first time step as these don't change
		if (t == 0)
		{

			/***********************/
			/***** BLOCK LABELS ****/
			/***********************/

#ifdef L_BUILD_FOR_MPI

			// Generate this data on the fly since all the same label and only done once
			std::vector<int> blockLabels(N_lim * M_lim * K_lim, mpim->my_rank);
			variable_name = time_string + "/MpiBlock";
			dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			hdf5_writeDataSet(memspace, filespace, dataset_id, eScalar, this, &blockLabels[0], H5T_NATIVE_INT, TL_present, TL_thickness, &minEdges[0], p_data);
			status = H5Dclose(dataset_id); // Close dataset
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif

			/***********************/
			/****** POSITIONS ******/
			/***********************/

			// WRITE POSITION X
			variable_name = time_string + "/XPos";
			dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			hdf5_writeDataSet(memspace, filespace, dataset_id, ePosX, this, &XPos[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
			status = H5Dclose(dataset_id); // Close dataset
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;

			// WRITE POSITION Y
			variable_name = time_string + "/YPos";
			dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			hdf5_writeDataSet(memspace, filespace, dataset_id, ePosY, this, &YPos[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
			status = H5Dclose(dataset_id); // Close dataset
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;


			// WRITE POSITION Z
#if (L_DIMS == 3)
			variable_name = time_string + "/ZPos";
			dataset_id = H5Dcreate(file_id, variable_name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			hdf5_writeDataSet(memspace, filespace, dataset_id, ePosZ, this, &ZPos[0], H5T_NATIVE_DOUBLE, TL_present, TL_thickness, &minEdges[0], p_data);
			status = H5Dclose(dataset_id); // Close dataset
			if (status != 0) *GridUtils::logfile << "HDF5 ERROR: Close dataset failed: " << status << std::endl;
#endif

		}

#ifdef L_BUILD_FOR_MPI
		// Synchronise before closing anything
		if (level == 0) MPI_Barrier(mpim->world_comm);
		else MPI_Barrier(mpim->subGrid_comm[(level - 1) + region_number * L_NUM_LEVELS]);
#endif

#ifdef L_HDF_DEBUG
		// Signal write completion
		L_INFO("Writing finished. Closing files...", GridUtils::logfile);
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

#ifdef L_BUILD_FOR_MPI
	}

	// No writable data
	else
	{
#ifdef L_MPI_VERBOSE
		L_INFO("Skipping HDF5 write as no writable data on L" + std::to_string(level) + " R" + std::to_string(region_number) + "...", mpim->logout);
#endif

#ifdef L_HDF_DEBUG
		// Check that the communicator was setup properly
		if (mpim->subGrid_comm[(level - 1) + region_number * L_NUM_LEVELS] != MPI_COMM_NULL)
		{
			L_ERROR("Communicator has a non-null value despite having no writable data: " +
				std::to_string(mpim->subGrid_comm[(level - 1) + region_number * L_NUM_LEVELS]),
				GridUtils::logfile);
		}
#endif	// L_HDF_DEBUG

	}
#endif	// L_BUILD_FOR_MPI

#ifdef L_HDF_DEBUG
	// Signal write completion
	L_WARN("Write out on this grid complete.", GridUtils::logfile);
#endif

	// Try call recursively on any present sub-grids
	for (GridObj *g : subGrid) g->io_hdf5(tval);	

	return 0;

}
// ***************************************************************************//

