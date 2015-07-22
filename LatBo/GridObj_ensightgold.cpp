#include "stdafx.h"
#include <sstream>
#include <iostream>
#include <iomanip>

#include "definitions.h"
#include "globalvars.h"
#include "GridObj.h"

using namespace std;

// Routine to generate the case file
void GridObj::ensight_gen_case(int nsteps)
{
	
	// Create file stream
	ofstream fout;
	
	// Create stringstream
	stringstream fileName;
	fileName.str("./Output/LBM.case");

	fout.open( fileName.str().c_str() );

	// Format section
	fout << "FORMAT" << endl;
	fout << "type: ensight gold" <<endl;

	// Geometry section
	fileName.str("geom.geo");
	fout << "GEOMETRY" << endl ;
	fout << "model: " << fileName.str().c_str() << endl;
	
	// Variable section
	// Only velocity (vector) and density (scalar) exported
	fout << "VARIABLE" << endl;
	fileName.str("velocity.*****.vec");
	fout << "vector per node: 1 velocity " << fileName.str().c_str() << endl;
	fileName.str("rho.*****.sca");
	fout << "scalar per node: 1 density  " << fileName.str().c_str() << endl;

	// Time section
	fout << "TIME" << endl;
	fout << "time set: 1" << endl;
	fout << "number of steps: " << nsteps / out_every << endl;	
	fout << "filename start number: 00000" << endl;
	fout << "filename increment: 1" << endl;
	fout << "time values: ";
	// Put in the time values
	for( int i = 0; i < nsteps / out_every; i++ ) {
		fout << 0 + i*dt*out_every << endl;
	}

	fout << endl ;

	fout.close() ;

}



// Routine to generate the geometry file
void GridObj::ensight_gen_geometry( )
{

	// Open file stream
	ofstream fout;

	char fileName[50];
	char buf[80]; // Character buffer
	sprintf(fileName,(char*)"./Output/geom.geo");
	
	if (level == 0) fout.open(fileName, ios::out);
	else fout.open(fileName, ios::out|ios::app); // Append

	// Only write the following header the first time round
	if (level == 0) {
		// Output format not required for ASCII out
			
		// Description line 1
		fout << "desc1";
	
		// Description line 2
		fout << "\ndesc2";
	
		// Tell it to assign node IDs for me
		fout << "\nnode id assign";

		// Tell it to assign element IDs for me
		fout << "\nelement id assign";
	}
	
	// Part header
	sprintf(buf, (char*)"\npart");
	fout << buf;
	// Unique part ID based on region and level
	int gridNum = level + region_number*NumLev;
	sprintf(buf, "\n%10d", gridNum); fout << buf;
	string buffer_str = "\nDesc: Grid Level = " + to_string(level) + 
							" Region = " + to_string(region_number);
	fout << buffer_str;
	
	// Dimensions of the grid
	sprintf(buf, (char*)"\nblock uniform iblanked"); // Uniform grid with blanking
	fout << buf;
	int ni = XPos.size();
	int nj = YPos.size();
#if (dims == 3)
	int nk = ZPos.size();
#else
	int nk = 1;
#endif

	// Print grid size (3 numbers on same line)
	sprintf(buf, "\n%10d ", ni); fout << buf;
	sprintf(buf, "%10d ", nj); fout << buf;
	sprintf(buf, "%10d", nk); fout << buf;

	// Origin location
	float x_orig = (float)XPos[0], y_orig = (float)YPos[0], z_orig = (float)ZPos[0];

	sprintf(buf, "\n%12.5e", (float)x_orig); fout << buf;
	sprintf(buf, "\n%12.5e", (float)y_orig); fout << buf;
	sprintf(buf, "\n%12.5e", (float)z_orig); fout << buf;

	// Grid spacing
	float x_delta = (float)dx, y_delta = (float)dy, z_delta = (float)dz;

	sprintf(buf, "\n%12.5e", (float)x_delta); fout << buf;
	sprintf(buf, "\n%12.5e", (float)y_delta); fout << buf;
	sprintf(buf, "\n%12.5e", (float)z_delta); fout << buf;

	// Node blanking flags
	for (int k = 0; k < nk; k++) {
		for (int j = 0; j < nj; j++) {
			for (int i = 0; i < ni; i++) {
				if (LatTyp(i,j,k,nj,nk) == 1 || LatTyp(i,j,k,nj,nk) == 4) {
					sprintf(buf, "\n%10d", 1); fout << buf;
				} else { 
					sprintf(buf, "\n%10d", 0); fout << buf;
				}
			}
		}
	}

	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].ensight_gen_geometry();
		}
	}

}


// Routine to generate the vectors file
void GridObj::ensight_gen_vector(int fileNum)
{

	// Open file stream
	ofstream fout;
	char fileName[50] ;
	char buf[80]; // Character buffer
	sprintf(fileName, (char*)"./Output/velocity.%05d.vec", fileNum);

	if (level == 0) fout.open(fileName, ios::out);
	else fout.open(fileName, ios::out|ios::app); // Append

	if (level == 0) {
		// Description line
		sprintf(buf, (char*)"required description");
		fout << buf;
	}
	
	// Part header
	sprintf(buf, (char*)"\npart");
	fout << buf;

	int gridNum = level + region_number*NumLev;
	sprintf(buf, "\n%10d", gridNum); fout << buf;
	
	sprintf(buf, (char*)"\nblock");
	fout << buf;

	// Data for part
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();

	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				int dir = 0;
				sprintf(buf, "\n%12.5e", (float)u(i,j,k,dir,M_lim,K_lim,dims)); fout << buf;
			}
		}
	}

	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				int dir = 1;
				sprintf(buf, "\n%12.5e", (float)u(i,j,k,dir,M_lim,K_lim,dims)); fout << buf;
			}
		}
	}

#if (dims == 3)
	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				int dir = 2;
				sprintf(buf, "\n%12.5e", (float)u(i,j,k,dir,M_lim,K_lim,dims)); fout << buf;
			}
		}
	}
#else
	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				sprintf(buf, "\n%12.5e", 0.0); fout << buf;
			}
		}
	}
#endif

	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].ensight_gen_vector(fileNum);
		}
	}

}


// Routine to generate the scalars file
void GridObj::ensight_gen_scalar(int fileNum)
{

	// File stream
	ofstream fout;
	char fileName[50];
	char buf[80];
	sprintf(fileName, (char*)"./Output/rho.%05d.sca", fileNum);

	if (level == 0 ) fout.open(fileName, ios::out);
	else fout.open(fileName, ios::out|ios::app);

	if (level == 0) {
		// Description
		sprintf(buf, (char*)"required description");
		fout << buf;
	}
	
	// Part header
	sprintf(buf, (char*)"\npart");
	fout << buf;

	int gridNum = level + region_number*NumLev;
	sprintf(buf, "\n%10d", gridNum); fout << buf;

	sprintf(buf, (char*)"\nblock");
	fout << buf;

	// Data
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();

	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				sprintf(buf, "\n%12.5e", (float)rho(i,j,k,M_lim,K_lim)); fout << buf;
			}
		}
	}

	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].ensight_gen_scalar(fileNum);
		}
	}

}
