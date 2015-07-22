/* This file contains the routines necessary to write out to a vtk file
*/

#include "stdafx.h"
#include <sstream>
#include <iomanip>

#include "definitions.h"
#include "globalvars.h"
#include "GridObj.h"

using namespace std;

// Routine to write out the vtk for time step t
void GridObj::vtk_writer(int t, double tval)
{
	
	// Create file name then output file stream
	stringstream fileName;
	fileName << "./Output/vtk_out.Lev" << level << "Reg" << region_number << "." << t << ".vtk";
	
	ofstream fout;
	fout.open( fileName.str().c_str() );

	// Add header information
	fout << "# vtk DataFile Version 3.0f\n";
	fout << "LBM Output at time t = " << t << "\n";
	fout << "ASCII\n";

	// Grid information -- structured points for uniform lattice
	size_t ni = XPos.size();
	size_t nj = YPos.size();
#if (dims == 3)
	size_t nk = ZPos.size();
#else
	size_t nk = 1;
#endif
	fout << "DATASET STRUCTURED_POINTS\n";

	// Grid dimensions continued
	fout << "DIMENSIONS " << ni << " " << nj << " " << nk << "\n";
	fout << "SPACING " << dx << " " << dy << " " << dz << "\n";
	fout << "ORIGIN " << XPos[0] << " " << YPos[0] << " " << ZPos[0] << "\n";

	// Data set
	fout << "POINT_DATA " << ni * nj * nk << "\n";

	// Density
	fout << "SCALARS Density " << "float 1\n";	// Name of set = Density, type = float, components = 1
	fout << "LOOKUP_TABLE default"; // Required if not using a custom lookup table
	for (size_t k = 0; k < nk; k++) {
		for (size_t j = 0; j < nj; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < ni; i++) {
				fout << (float)rho(i,j,k,nj,nk) << " ";
			}
		}
	}

	// Velocity
	fout << "\nVECTORS Velocity " << "float";	// Name of set = Velocity, type = float
#if (dims == 3)
	for (size_t k = 0; k < nk; k++) {
		for (size_t j = 0; j < nj; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < ni; i++) {
				for (size_t dir = 0; dir < dims; dir++) {
					fout << (float)u(i,j,k,dir,nj,nk,dims) << " ";
				}
			}
		}
	}
#else
	for (size_t k = 0; k < nk; k++) {
		for (size_t j = 0; j < nj; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < ni; i++) {
				size_t dir = 0;
				fout << (float)u(i,j,k,dir,nj,nk,dims) << " ";
				dir = 1;
				fout << (float)u(i,j,k,dir,nj,nk,dims) << " ";
				// Set z as zero
				fout << 0.0 << " ";
			}
		}
	}
#endif


	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].vtk_writer(t, tval);
		}
	}

	return;

}
