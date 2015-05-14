#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <iomanip>

#include "definitions.h"
#include "globalvars.h"
#include "GridObj.h"

using namespace std;

// Routine to generate the case file
void GridObj::genCase(int nsteps, int saveEvery)
{
	
	// Create file stream
	ofstream fout;
	
	char fileName[50];
	sprintf_s(fileName, (char*)"./Output/LBM.case");

	fout.open(fileName);

	// Format section
	fout << "FORMAT" << endl;
	fout << "type: ensight gold" <<endl;

	// Geometry section
	sprintf_s(fileName, (char*)"geom.geo");
	fout << "GEOMETRY" << endl ;
	fout << "model: " << fileName << endl;
	
	// Varible section
	// Only velocity (vector) and density (scalar) exported
	fout << "VARIABLE" << endl;
	sprintf_s(fileName, (char*)"velocity.*****.vec");
	fout << "vector per element: 1 velocity " << fileName << endl;
	sprintf_s(fileName, (char*)"rho.*****.sca");
	fout << "scalar per element: 1 density  " << fileName << endl;

	// Time section
	fout << "TIME" << endl;
	fout << "time set: 1" << endl;
	fout << "number of steps: " << nsteps << endl;	
	fout << "filename start number: 00000" << endl;
	fout << "filename increment: 1" << endl;
	fout << "time values: ";
	// Put in the time values
	for( int i = 0; i < nsteps; i++ ) {
		fout << 1+i*saveEvery << endl;
	}

	fout << endl ;

	fout.close() ;

}



// Routine to generate the geometry file
void GridObj::genGeo( )
{

	// Open file stream
	ofstream fout;

	char fileName[50];
	char buf[80]; // Character buffer
	sprintf_s(fileName,(char*)"./Output/geom.geo");
	
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
	sprintf_s(buf, (char*)"\npart");
	fout << buf;
	// Unique part ID based on region and level
	int gridNum = level + region_number*NumLev;
	sprintf_s(buf, "\n%10d", gridNum); fout << buf;
	string buffer_str = "\nDesc: Grid Level = " + to_string(level) + 
							" Region = " + to_string(region_number);
	fout << buffer_str;
	
	// Dimensions of the grid
	sprintf_s(buf, (char*)"\nblock uniform"); // Uniform grid
	fout << buf;
	int ni = XPos.size()+1;
	int nj = YPos.size()+1;
#if (dims == 3)
	int nk = ZPos.size()+1;
#else
	int nk = 1;
#endif

	// Print grid size (3 numbers on same line)
	sprintf_s(buf, "\n%10d ", ni); fout << buf;
	sprintf_s(buf, "%10d ", nj); fout << buf;
	sprintf_s(buf, "%10d", nk); fout << buf;

	// Origin location
	float x_orig = (float)XPos[0], y_orig = (float)YPos[0], z_orig = (float)ZPos[0];

	sprintf_s(buf, "\n%12.5e", (float)x_orig); fout << buf;
	sprintf_s(buf, "\n%12.5e", (float)y_orig); fout << buf;
	sprintf_s(buf, "\n%12.5e", (float)z_orig); fout << buf;

	float x_delta = (float)dx, y_delta = (float)dy, z_delta = (float)dz;

	sprintf_s(buf, "\n%12.5e", (float)x_delta); fout << buf;
	sprintf_s(buf, "\n%12.5e", (float)y_delta); fout << buf;
	sprintf_s(buf, "\n%12.5e", (float)z_delta); fout << buf;

	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].genGeo();
		}
	}

}


// Routine to generate the vectors file
void GridObj::genVec(int fileNum)
{

	// Open file stream
	ofstream fout;
	char fileName[50] ;
	char buf[80]; // Character buffer
	sprintf_s(fileName, (char*)"./Output/velocity.%05d.vec", fileNum);

	if (level == 0) fout.open(fileName, ios::out);
	else fout.open(fileName, ios::out|ios::app); // Append

	if (level == 0) {
		// Description line
		sprintf_s(buf, (char*)"required description");
		fout << buf;
	}
	
	// Part header
	sprintf_s(buf, (char*)"\npart");
	fout << buf;

	int gridNum = level + region_number*NumLev;
	sprintf_s(buf, "\n%10d", gridNum); fout << buf;
	
	sprintf_s(buf, (char*)"\nblock");
	fout << buf;

	// Data for part
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();
	int ct;
	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				int dir = 0;
				ct = idxmap(i,j,k,dir,M_lim,K_lim,dims);
				sprintf_s(buf, "\n%12.5e", (float)u[ct]); fout << buf;
			}
		}
	}

	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				int dir = 1;
				ct = idxmap(i,j,k,dir,M_lim,K_lim,dims);
				sprintf_s(buf, "\n%12.5e", (float)u[ct]); fout << buf;
			}
		}
	}

#if (dims == 3)
	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				int dir = 2;
				ct = idxmap(i,j,k,dir,M_lim,K_lim,dims);
				sprintf_s(buf, "\n%12.5e", (float)u[ct]); fout << buf;
			}
		}
	}
#else
	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				sprintf_s(buf, "\n%12.5e", 0.0); fout << buf;
			}
		}
	}
#endif

	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].genVec(fileNum);
		}
	}

}


// Routine to generate the scalars file
void GridObj::genScal(int fileNum)
{

	// File stream
	ofstream fout;
	char fileName[50];
	char buf[80];
	sprintf_s(fileName, (char*)"./Output/rho.%05d.sca", fileNum);

	if (level == 0 ) fout.open(fileName, ios::out);
	else fout.open(fileName, ios::out|ios::app);

	if (level == 0) {
		// Description
		sprintf_s(buf, (char*)"required description");
		fout << buf;
	}
	
	// Part header
	sprintf_s(buf, (char*)"\npart");
	fout << buf;

	int gridNum = level + region_number*NumLev;
	sprintf_s(buf, "\n%10d", gridNum); fout << buf;

	sprintf_s(buf, (char*)"\nblock");
	fout << buf;

	// Data
	int N_lim = XPos.size();
	int M_lim = YPos.size();
	int K_lim = ZPos.size();
	int ct;
	for (int k = 0; k < K_lim; k++) {
		for (int j = 0; j < M_lim; j++) {
			for (int i = 0; i < N_lim; i++) {
				ct = idxmap(i,j,k,M_lim,K_lim);
				sprintf_s(buf, "\n%12.5e", (float)rho[ct]); fout << buf;
			}
		}
	}

	fout.close();

	// Now do the rest of the grids
	if (NumLev > level) {
		for (size_t reg = 0; reg < subGrid.size(); reg++) {
			subGrid[reg].genScal(fileNum);
		}
	}

}
