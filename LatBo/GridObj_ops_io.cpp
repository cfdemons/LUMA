// Routines for reading and writing operations.

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "GridObj.h"
#include "definitions.h"
#include "globalvars.h"

using namespace std;

// ***************************************************************************************************
// Writes all the contents of the class at time t and call recursviely for any subgrids.
// Writes to text file "Grids.out" by default.
void GridObj::io_textout() {

	// Get limits for current level
	size_t N_lim = XPos.size();
	size_t M_lim = YPos.size();
	size_t K_lim = ZPos.size();

	// Create stream and open text file
	ofstream gridoutput;
	gridoutput.precision(4);

	// Construct File Name
	string FNameG, N_str, M_str, K_str, ex_str, NumLev_str, NumReg_str, mpirank;
	N_str = to_string((int)N);
	M_str = to_string((int)M);
	K_str = to_string((int)K);
	NumLev_str = to_string(NumLev);
	if (NumLev == 0) ex_str = to_string(0);
	else ex_str = to_string(RefXend[0]) + string("_") + to_string(RefYend[0]) + string("_") + to_string(RefZend[0]);
	if (NumLev == 0) NumReg_str = to_string(0);
	else NumReg_str = to_string(NumReg);
	mpirank = to_string(my_rank);
	// Build string
	FNameG = string("./Output/Grids")
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
	if (t == 0 && level == 0) gridoutput.open(FNameG_c, ios::out);
	else gridoutput.open(FNameG_c, ios::out |ios::app);

	if ( gridoutput.is_open() ) {

		// Draw a line to begin
		gridoutput << "\n-------------------------------------------------------------------------------------" << endl;
		gridoutput << "-----------------------------------START OF OUTPUT-----------------------------------" << endl;
		gridoutput << "-------------------------------------------------------------------------------------" << endl;

		if (level == 0) {
			// Print L0 Grid Size header
			gridoutput << "L0 Grid Size = " << N << " x " << M << " x " << K << endl;
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
				
				subGrid[reg].io_textout();

			}
		}
		

		// Close file
		gridoutput.close();

	} else {

		cout << "Cannot open file" << endl;

	}

}


// ***************************************************************************************************
// Routine to write out the coordinates of IBbodies at a given time step
void GridObj::io_write_body_pos() {

	// Get time step
	unsigned int timestep = t;

	for (size_t ib = 0; ib < iBody.size(); ib++) {
		

			// Open file for given time step
			std::ofstream jout;
			jout.open("./Output/Body_" + to_string(ib) + "_position_" + to_string(timestep) + "_rank" + std::to_string(my_rank) + ".out", std::ios::out);
			jout << "x" + to_string(timestep) + ", y" + to_string(timestep) + ", z" << std::endl;
	
			// Write out position
			for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
#if (dims == 3)
				jout << iBody[ib].markers[i].position[0] << ", " << iBody[ib].markers[i].position[1] << ", " << iBody[ib].markers[i].position[2] << std::endl;
#else
				jout << iBody[ib].markers[i].position[0] << ", " << iBody[ib].markers[i].position[1] << ", " << 0.0 << std::endl;
#endif
			}
			jout.close();

	}

}


// ***************************************************************************************************
// Routine to write out the coordinates of IBbodies at a given time step
void GridObj::io_write_lift_drag() {	

	// Get time step
	unsigned int timestep = t;

	for (size_t ib = 0; ib < iBody.size(); ib++) {
		

			// Open file for given time step
			std::ofstream jout;
			jout.open("./Output/Body_" + to_string(ib) + "_LD_" + to_string(timestep) + "_rank" + std::to_string(my_rank) + ".out", std::ios::out);
			jout << "L" + to_string(timestep) + ", D" + to_string(timestep) << std::endl;

			// Sum variables
			double Lsum = 0.0, Dsum = 0.0;

			// Compute lift and drag
			for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
				jout << iBody[ib].markers[i].force_xyz[0] << ", " << iBody[ib].markers[i].force_xyz[1] << std::endl;
				Lsum += iBody[ib].markers[i].force_xyz[0];
				Dsum += iBody[ib].markers[i].force_xyz[1];
			}

			jout << "Totals = " << std::endl;
			jout << Lsum << ", " << Dsum << std::endl;
			jout.close();

	}

}

// ***************************************************************************************************