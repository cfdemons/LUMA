/* This file contains the routines necessary to write out to a vtk file
*/

#include "../inc/stdafx.h"
#include <sstream>
#include <iomanip>
#include "../inc/definitions.h"
#include "../inc/globalvars.h"
#include "../inc/GridObj.h"
#include "../inc/MpiManager.h"

using namespace std;

// Routine to write out the vtk for time step t
void GridObj::io_vtkwriter(double tval)
{

	// Create file name then output file stream
	stringstream fileName;
	fileName << GridUtils::path_str + "/vtk_out.Lev" << level << ".Reg" << region_number << 
		".Rnk" << MpiManager::my_rank << "." << (int)tval << ".vtk";

	ofstream fout;
	fout.open( fileName.str().c_str() );

	// Add header information
	fout << "# vtk DataFile Version 3.0f\n";
	fout << "LBM Output at time t = " << (int)tval << "\n";
	fout << "ASCII\n";

	// Grid information -- structured points for uniform lattice
	size_t ni, nj, nk;
	size_t ni_corrected = 0, nj_corrected = 0, nk_corrected = 0;
	size_t start_x = 0, start_y = 0, start_z = 0;

	// Actual grid sizes
	ni = XPos.size();
	nj = YPos.size();
#if (dims == 3)
	nk = ZPos.size();
#else
	nk = 1;
	start_z = 0;
	nk_corrected = 1;
#endif

	// If using MPI, correct grid size to avoid writing out the receiving buffers
#ifdef BUILD_FOR_MPI

	// Loop over positions to determine size of grid without overlap
	for (size_t i = 0; i < ni; i++) {
		if(!GridUtils::isOnRecvLayer(XPos[i],'x',"min") && !GridUtils::isOnRecvLayer(XPos[i],'x',"max")) {
			ni_corrected++;
			if (ni_corrected == 1) {
				start_x = i;
			}
		}
	}
	for (size_t j = 0; j < nj; j++) {
		if(!GridUtils::isOnRecvLayer(YPos[j],'y',"min") && !GridUtils::isOnRecvLayer(YPos[j],'y',"max")) {
			nj_corrected++;
			if (nj_corrected == 1) {
				start_y = j;
			}
		}
	}
#if (dims == 3)
	for (size_t k = 0; k < nk; k++) {
		if(!GridUtils::isOnRecvLayer(ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(ZPos[k],'z',"max")) {
			nk_corrected++;
			if (nk_corrected == 1) {
				start_z = k;
			}
		}
	}
#endif

#else
	// Not using MPI so no overlaps to consider
	ni_corrected = XInd.size();
	nj_corrected = YInd.size();

#if (dims == 3)
	nk_corrected = ZInd.size();
#endif

#endif

	fout << "DATASET STRUCTURED_POINTS\n";

	// Grid dimensions continued
	fout << "DIMENSIONS " << ni_corrected << " " << nj_corrected << " " << nk_corrected << "\n";
	fout << "SPACING " << dx << " " << dy << " " << dz << "\n";

	// Starting point
	fout << "ORIGIN " << XPos[start_x] << " " << YPos[start_y] << " " << ZPos[start_z] << "\n";

	// Data set
	fout << "POINT_DATA " << ni_corrected * nj_corrected * nk_corrected << "\n";

	// Density
	fout << "SCALARS Density " << "float 1\n";	// Name of set = Density, type = float, components = 1
	fout << "LOOKUP_TABLE default"; // Required if not using a custom lookup table
	for (size_t k = 0; k < ZInd.size(); k++) {
		for (size_t j = 0; j < YInd.size(); j++ ) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < XInd.size(); i++) {
#ifdef BUILD_FOR_MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif
				{
					fout << (float)rho(i,j,k,nj,nk) << " ";
				}
			}
		}
	}

	// Time-Averaged Density
	fout << "\nSCALARS TimeAveragedDensity " << "float 1\n";	// Name of set = TimeAveragedDensity, type = float, components = 1
	fout << "LOOKUP_TABLE default"; // Required if not using a custom lookup table
	for (size_t k = 0; k < ZInd.size(); k++) {
		for (size_t j = 0; j < YInd.size(); j++ ) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < XInd.size(); i++) {
#ifdef BUILD_FOR_MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif
				{
					fout << (float)rho_timeav(i,j,k,nj,nk) << " ";
				}
			}
		}
	}

	// Time-Averaged UiUj (dims-1 components)
	fout << "\nSCALARS TimeAveragedUiUj " << "float " << to_string(2*dims-3) << "\n";	// Components = 2*dims-3
	fout << "LOOKUP_TABLE default"; // Required if not using a custom lookup table
	for (size_t k = 0; k < ZInd.size(); k++) {
		for (size_t j = 0; j < YInd.size(); j++ ) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < XInd.size(); i++) {
#ifdef BUILD_FOR_MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif
				{
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
	}


	// Velocity
	fout << "\nVECTORS Velocity " << "float";	// Name of set = Velocity, type = float
	for (size_t k = 0; k < ZInd.size(); k++) {
		for (size_t j = 0; j < YInd.size(); j++ ) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < XInd.size(); i++) {
#ifdef BUILD_FOR_MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif
				{
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
	}

	// Time-Averaged Ui
	fout << "\nVECTORS TimeAveragedUi " << "float";	// Name of set = TimeAveragedUi, type = float
	for (size_t k = 0; k < ZInd.size(); k++) {
		for (size_t j = 0; j < YInd.size(); j++ ) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < XInd.size(); i++) {
#ifdef BUILD_FOR_MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif
				{

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
	}

	// Time-Averaged UiUi
	fout << "\nVECTORS TimeAveragedUiUi " << "float";
	for (size_t k = 0; k < ZInd.size(); k++) {
		for (size_t j = 0; j < YInd.size(); j++ ) {
			fout << "\n"; // New line for each row
			for (size_t i = 0; i < XInd.size(); i++) {
#ifdef BUILD_FOR_MPI
				if (!GridUtils::isOnRecvLayer(XPos[i],YPos[j],ZPos[k]))
#endif
				{

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
	}


	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].io_vtkwriter(tval);
		}		
	}

	return;

}
