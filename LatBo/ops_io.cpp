// Routines for reading and writing operations.

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "LBM_globalvars.h"

using namespace std;

string int2str(int number); // Forward declaration for int2str


// 3D write out all the contents of the Grid structure for level r at time t.
// Writes to text file "Grids.out".
void lbm_write3(int r, int t) {

	// Get limits for current level
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();

	// Create stream and open text file
	ofstream gridoutput;
	gridoutput.precision(4);

	// Construct File Name
	string FNameG, N_str, M_str, K_str, ex_str, Nref_str, PartOpt_str;
	N_str = int2str((int)N);
	M_str = int2str((int)M);
	K_str = int2str((int)K);
	ex_str = int2str(Ref_endX);
	Nref_str = int2str(Nref);
	FNameG = string("Grids")
			+ string("D") +  int2str(dims)
			+ string("x") + N_str 
			+ string("y") + M_str 
			+ string("z") + K_str 
			+ string("Nref") + Nref_str 
			+ string("P") + ex_str 
			+ string(".out");
	const char* FNameG_c = FNameG.c_str();

	gridoutput.open(FNameG_c, ios::out |ios::app);

	if ( gridoutput.is_open() ) {

		// Draw a line to begin
		gridoutput << "\n-------------------------------------------------------------------------------------" << endl;
		gridoutput << "-----------------------------------START OF OUTPUT-----------------------------------" << endl;
		gridoutput << "-------------------------------------------------------------------------------------" << endl;

		// Global Data
		// Identify L0 Grid Size
		gridoutput << "L0 Grid Size = " << N << " x " << M << " x " << K << endl;
		int finex = Grids[0].XInd.size(); int finey = Grids[0].YInd.size(); int finez = Grids[0].ZInd.size();
		gridoutput << "Refined = " << (((float)finex)*((float)finey)*((float)finez)*100) / (N*M*K) << "%" << endl;
		gridoutput << "-------------------------------------------------------------------------------------" << endl;

		// Identify Grid Level
		string r_str = int2str(r);
		gridoutput << "Grid Level = " << r_str << endl;
		
		// Identify time step
		string t_str = int2str(t);
		gridoutput << "Time Step = " << t << endl;

		// Now print omega
		gridoutput << "Omega = " << Grids[r].omega << endl;

		// Index Vectors
		gridoutput << "X Index: ";
		for (size_t i = 0; i < Grids[r].XInd.size(); i++) {
			gridoutput << Grids[r].XInd[i] << "\t";
		}
		gridoutput << "\nY Index: ";
		for (size_t j = 0; j < Grids[r].YInd.size(); j++) {
			gridoutput << Grids[r].YInd[j] << "\t";
		}
		gridoutput << "\nZ Index: ";
		for (size_t k = 0; k < Grids[r].ZInd.size(); k++) {
			gridoutput << Grids[r].ZInd[k] << "\t";
		}
	
		// Position Vectors
		gridoutput << "\nX Position: \t";
		for (size_t i = 0; i < Grids[r].XPos.size(); i++) {
			gridoutput << Grids[r].XPos[i] << "\t";
		}
		gridoutput << "\nY Position: \t";
		for (size_t j = 0; j < Grids[r].YPos.size(); j++) {
			gridoutput << Grids[r].YPos[j] << "\t";
		}
		gridoutput << "\nZ Position: \t";
		for (size_t k = 0; k < Grids[r].ZPos.size(); k++) {
			gridoutput << Grids[r].ZPos[k] << "\t";
		}

		// Typing Matrix
		gridoutput << "\n\nTyping Matrix";
		for (size_t k = 0; k < Grids[r].ZPos.size(); k++) {
			// New line with z-coordinate
			gridoutput << "\nz = " << Grids[r].ZPos[k] << "\n";

			for (size_t i = 0; i < Grids[r].XPos.size(); i++) {
				// New line
				gridoutput << "\n";
				for (size_t j = 0; j < Grids[r].YPos.size(); j++) {

					// Index
					int idx = idxmap(i,j,k,M_lim,K_lim);
					gridoutput << Grids[r].LatTyp[idx] << "\t";

				}
			}
		}

		// Populations (f, feq)
		gridoutput << "\n\nf Values";
		for (size_t v = 0; v < nVels; v++) {
			// Particular velocity
			gridoutput << "\nc = " << v+1;


			for (size_t k = 0; k < Grids[r].ZPos.size(); k++) {
				// New line with z-coordinate
				gridoutput << "\nz = " << Grids[r].ZPos[k] << "\n";

				for (size_t i = 0; i < Grids[r].XPos.size(); i++) {
					// New line
					gridoutput << "\n";
					for (size_t j = 0; j < Grids[r].YPos.size(); j++) {

						// Index
						int idx = idxmap(i,j,k,v,M_lim,K_lim,nVels);
						gridoutput << Grids[r].f[idx] << "\t";

					}
				}
			}
		}

		gridoutput << "\n\nfeq Values";
		for (size_t v = 0; v < nVels; v++) {
			// Particular velocity
			gridoutput << "\nc = " << v+1;


			for (size_t k = 0; k < Grids[r].ZPos.size(); k++) {
				// New line with z-coordinate
				gridoutput << "\nz = " << Grids[r].ZPos[k] << "\n";

				for (size_t i = 0; i < Grids[r].XPos.size(); i++) {
					// New line
					gridoutput << "\n";
					for (size_t j = 0; j < Grids[r].YPos.size(); j++) {

						// Index
						int idx = idxmap(i,j,k,v,M_lim,K_lim,nVels);
						gridoutput << Grids[r].feq[idx] << "\t";

					}
				}
			}
		}

		// Macroscopic (u, rho)
		gridoutput << "\n\nVelocity Values";
		for (size_t n = 0; n < dims; n++) {
			// Particular component
			gridoutput << "\nu(" << n+1 << ")";
			
			for (size_t k = 0; k < Grids[r].ZPos.size(); k++) {
				// New line with z-coordinate
				gridoutput << "\nz = " << Grids[r].ZPos[k] << "\n";

				for (size_t i = 0; i < Grids[r].XPos.size(); i++) {
					// New line
					gridoutput << "\n";
					for (size_t j = 0; j < Grids[r].YPos.size(); j++) {

						// Index
						int idx = idxmap(i,j,k,n,M_lim,K_lim,dims);
						gridoutput << Grids[r].u[idx] << "\t";

					}
				}
			}
		}

		gridoutput << "\n\nDensity";
		for (size_t k = 0; k < Grids[r].ZPos.size(); k++) {
			// New line with z-coordinate
			gridoutput << "\nz = " << Grids[r].ZPos[k] << "\n";

			for (size_t i = 0; i < Grids[r].XPos.size(); i++) {
				// New line
				gridoutput << "\n";
				for (size_t j = 0; j < Grids[r].YPos.size(); j++) {

					// Index
					int idx = idxmap(i,j,k,M_lim,K_lim);
					gridoutput << Grids[r].rho[idx] << "\t";

				}
			}
		}

		// Draw a line underneath
		gridoutput << "\n-------------------------------------------------------------------------------------" << endl;
		gridoutput << "------------------------------------END OF OUTPUT------------------------------------" << endl;
		gridoutput << "-------------------------------------------------------------------------------------" << endl;

		// Close file
		gridoutput.close();

	} else {

		cout << "Cannot open file" << endl;

	}

}