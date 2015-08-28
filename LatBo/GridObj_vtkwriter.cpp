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
	fileName << "./Output/vtk_out.Lev" << level << ".Reg" << region_number << ".Rnk" << my_rank << "." << t << ".vtk";
	
	ofstream fout;
	fout.open( fileName.str().c_str() );

	// Add header information
	fout << "# vtk DataFile Version 3.0f\n";
	fout << "LBM Output at time t = " << t << "\n";
	fout << "ASCII\n";

	// Grid information -- structured points for uniform lattice
	size_t ni, nj, nk, ni_corrected, nj_corrected, nk_corrected, startx, starty, startz, endx, endy, endz;

	// Actual grid sizes
	ni = XPos.size();
	nj = YPos.size();
#if (dims == 3)
	nk = ZPos.size();
#else
	nk = 1;
#endif

	// If using MPI correct grid size and starting and end loop 
	// indices to avoid writing out the outer buffers
#ifdef BUILD_FOR_MPI
	ni_corrected = ni - 2;
	nj_corrected = nj - 2;
	startx = 1; endx = ni_corrected + 1;
	starty = 1; endy = nj_corrected + 1;
#if (dims == 3)
	nk_corrected = nk - 2;
	startz = 1; endz = nk_corrected + 1;
#else
	nk_corrected = nk;
	startz = 0; endz = nk_corrected;
#endif

#else
	ni_corrected = ni;
	nj_corrected = nj;
	nk_corrected = nk;
	startx = 0; endx = ni_corrected;
	starty = 0; endy = nj_corrected;
	startz = 0; endz = nk_corrected;

#endif

	fout << "DATASET STRUCTURED_POINTS\n";

	// Grid dimensions continued
	fout << "DIMENSIONS " << ni_corrected << " " << nj_corrected << " " << nk_corrected << "\n";
	fout << "SPACING " << dx << " " << dy << " " << dz << "\n";

	// Even refined grids when using MPI need to start at zero
	if ( level == 0 ) fout << "ORIGIN " << XPos[startx] << " " << YPos[starty] << " " << ZPos[startz] << "\n";
	else fout << "ORIGIN " << XPos[0] << " " << YPos[0] << " " << ZPos[0] << "\n";

	// Data set
	fout << "POINT_DATA " << ni_corrected * nj_corrected * nk_corrected << "\n";

	// Density
	fout << "SCALARS Density " << "float 1\n";	// Name of set = Density, type = float, components = 1
	fout << "LOOKUP_TABLE default"; // Required if not using a custom lookup table
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {
				fout << (float)rho(i,j,k,nj,nk) << " ";
			}
		}
	}

	// Velocity
	fout << "\nVECTORS Velocity " << "float";	// Name of set = Velocity, type = float
#if (dims == 3)
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {
				for (size_t dir = 0; dir < dims; dir++) {
					fout << (float)u(i,j,k,dir,nj,nk,dims) << " ";
				}
			}
		}
	}
#else
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {
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
