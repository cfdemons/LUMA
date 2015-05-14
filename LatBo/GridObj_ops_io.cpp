// Routines for reading and writing operations.

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "GridObj.h"
#include "definitions.h"
#include "globalvars.h"

using namespace std;

string int2str(int number); // Forward declaration for int2str


// Writes all the contents of the class at time t and call recursviely for any subgrids.
// Writes to text file "Grids.out" by default.
void GridObj::lbm_write3(int t) {

	// Get limits for current level
	size_t N_lim = XPos.size();
	size_t M_lim = YPos.size();
	size_t K_lim = ZPos.size();

	// Create stream and open text file
	ofstream gridoutput;
	gridoutput.precision(4);

	// Construct File Name
	string FNameG, N_str, M_str, K_str, ex_str, NumLev_str, NumReg_str;
	N_str = int2str((int)N);
	M_str = int2str((int)M);
	K_str = int2str((int)K);
	NumLev_str = int2str(NumLev);
	if (NumLev == 0) ex_str = int2str(0);
	else ex_str = int2str(RefXend[0]) + string("_") + int2str(RefYend[0]) + string("_") + int2str(RefZend[0]);
	if (NumLev == 0) NumReg_str = int2str(0);
	else NumReg_str = int2str(NumReg);
	// Build string
	FNameG = string("./Output/Grids")
			+ string("D") +  int2str(dims)
			+ string("x") + N_str 
			+ string("y") + M_str 
			+ string("z") + K_str 
			+ string("Lev") + NumLev_str 
			+ string("Reg") + NumReg_str 
			+ string("P") + ex_str 
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
			if (NumLev != 0) {
				// Get size of regions
				for (int reg = 0; reg < NumReg; reg++) {
					int finex = subGrid[reg].CoarseLimsX[1] - subGrid[reg].CoarseLimsX[0] + 1;
					int finey = subGrid[reg].CoarseLimsY[1] - subGrid[reg].CoarseLimsY[0] + 1;
					int finez = subGrid[reg].CoarseLimsZ[1] - subGrid[reg].CoarseLimsZ[0] + 1;
					gridoutput << "Region # " << reg << " refinement = " << (((float)finex)*((float)finey)*((float)finez)*100) / (N*M*K) << "%" << endl;
				}
			}
			// Print time step
			string t_str = int2str(t);
			gridoutput << "Time Step = " << t << endl;
			gridoutput << "-------------------------------------------------------------------------------------" << endl;
		}

		// Print Grid Level
		string r_str = int2str(level);
		gridoutput << "Grid Level = " << r_str << endl;

		// Print region number
		string reg_str = int2str(region_number);
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

					// Index
					int idx = idxmap(i,j,k,M_lim,K_lim);
					gridoutput << LatTyp[idx] << "\t";

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

						// Index
						int idx = idxmap(i,j,k,v,M_lim,K_lim,nVels);
						gridoutput << f[idx] << "\t";

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

						// Index
						int idx = idxmap(i,j,k,v,M_lim,K_lim,nVels);
						gridoutput << feq[idx] << "\t";

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

						// Index
						int idx = idxmap(i,j,k,n,M_lim,K_lim,dims);
						gridoutput << u[idx] << "\t";

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

					// Index
					int idx = idxmap(i,j,k,M_lim,K_lim);
					gridoutput << rho[idx] << "\t";

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

				subGrid[reg].lbm_write3(t);

			}
		}
		

		// Close file
		gridoutput.close();

	} else {

		cout << "Cannot open file" << endl;

	}

}