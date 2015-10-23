/* This file contains the routines necessary to write out to a vtk file
*/

#include "../inc/stdafx.h"
#include <sstream>
#include <iomanip>

#include "../inc/definitions.h"
#include "../inc/globalvars.h"
#include "../inc/GridObj.h"

using namespace std;

// Routine to write out the vtk for time step t
void GridObj::vtk_writer(double tval)
{

	// Create file name then output file stream
	stringstream fileName;
	fileName << "./output/vtk_out.Lev" << level << ".Reg" << region_number << ".Rnk" << my_rank << "." << (int)tval << ".vtk";

	ofstream fout;
	fout.open( fileName.str().c_str() );

	// Add header information
	fout << "# vtk DataFile Version 3.0f\n";
	fout << "LBM Output at time t = " << (int)tval << "\n";
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

	// Time-Averaged Density
	fout << "\nSCALARS TimeAveragedDensity " << "float 1\n";	// Name of set = TimeAveragedDensity, type = float, components = 1
	fout << "LOOKUP_TABLE default"; // Required if not using a custom lookup table
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {
				fout << (float)rho_timeav(i,j,k,nj,nk) << " ";
			}
		}
	}

	// Time-Averaged UiUj (dims-1 components)
	fout << "\nSCALARS TimeAveragedUiUj " << "float " << to_string(2*dims-3) << "\n";	// Components = 2*dims-3
	fout << "LOOKUP_TABLE default"; // Required if not using a custom lookup table
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {

				// Write out element 1 in 2D (u_01)
				fout << (float)uiuj_timeav(i,j,k,1,nj,nk,(3*dims-3)) << " ";
#if (dims == 3)
				// Also write out element 2 and 4 in 3D (u_02 and u_12 respectively)
				fout << (float)uiuj_timeav(i,j,k,2,nj,nk,(3*dims-3)) << " ";
				fout << (float)uiuj_timeav(i,j,k,4,nj,nk,(3*dims-3)) << " ";
#endif

			}
		}
	}


	// Velocity
	fout << "\nVECTORS Velocity " << "float";	// Name of set = Velocity, type = float
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {

				for (size_t dir = 0; dir < dims; dir++) {
					fout << (float)u(i,j,k,dir,nj,nk,dims) << " ";
				}

#if (dims != 3)
				// Need to add a z-component of zero to complete file in 2D
				fout << 0.0 << " ";
#endif

			}
		}
	}

	// Time-Averaged Ui
	fout << "\nVECTORS TimeAveragedUi " << "float";	// Name of set = TimeAveragedUi, type = float
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {

				for (size_t dir = 0; dir < dims; dir++) {
					fout << (float)ui_timeav(i,j,k,dir,nj,nk,dims) << " ";
				}

#if (dims != 3)
				// Need to add a z-component of zero to complete file in 2D
				fout << 0.0 << " ";
#endif

			}
		}
	}

	// Time-Averaged UiUi
	fout << "\nVECTORS TimeAveragedUiUi " << "float";
	for (size_t k = startz; k < endz; k++) {
		for (size_t j = starty; j < endy; j++) {
			fout << "\n"; // New line for each row
			for (size_t i = startx; i < endx; i++) {

				// Write out element 0 in both 2D and 3D (u_00)
				fout << (float)uiuj_timeav(i,j,k,0,nj,nk,(3*dims-3)) << " ";

#if (dims == 3)
				// Also write out elements 3 and 5 in 3D (u_11 and u_22 respectively)
				fout << (float)uiuj_timeav(i,j,k,3,nj,nk,(3*dims-3)) << " ";
				fout << (float)uiuj_timeav(i,j,k,5,nj,nk,(3*dims-3)) << " ";
#else
				// Also write out element 2 in 2D (u_11)
				fout << (float)uiuj_timeav(i,j,k,2,nj,nk,(3*dims-3)) << " ";
				// Need to add a z-component of zero to complete file in 2D
				fout << 0.0 << " ";
#endif

			}
		}
	}


	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].vtk_writer(tval);
		}
	}

	return;

}


// ***************************************************************************************************
// Routine to write out the vtk (position) for each IB body at time step t (current capability is for unclosed objects only)
void GridObj::vtk_IBwriter(double tval) {

    // Loop through each iBody
    for (size_t ib = 0; ib < iBody.size(); ib++) {

        // Create file name then output file stream
        stringstream fileName;
        fileName << "./output/vtk_IBout.Body" << ib << "." << (int)tval << ".vtk";

        ofstream fout;
        fout.open( fileName.str().c_str() );

        // Add header information
        fout << "# vtk DataFile Version 3.0f\n";
        fout << "IB Output for body ID " << ib << " at time t = " << (int)tval << "\n";
        fout << "ASCII\n";
        fout << "DATASET POLYDATA\n";


        // Write out the positions of each Lagrange marker
        fout << "POINTS " << iBody[ib].markers.size() << " float\n";
        for (size_t i = 0; i < iBody[ib].markers.size(); i++) {

#if (dims == 3)
				fout << iBody[ib].markers[i].position[0] << " " << iBody[ib].markers[i].position[1] << " " << iBody[ib].markers[i].position[2] << std::endl;
#else
				fout << iBody[ib].markers[i].position[0] << " " << iBody[ib].markers[i].position[1] << " " << 1.0 << std::endl; // z = 1.0 as fluid ORIGIN is at z = 1.0
#endif
        }


        // Write out the connectivity of each Lagrange marker
        size_t nLines = iBody[ib].markers.size() - 1;       // Non-closed surface only (for now)
        fout << "LINES " << nLines << " " << 3 * nLines << std::endl;

        for (size_t i = 0; i < nLines; i++) {
            fout << 2 << " " << i << " " << i + 1 << std::endl;
        }

        fout.close();
    }
}
