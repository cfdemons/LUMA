#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <iomanip>

#include "ensightGold.h"
#include "definitions.h"
#include "globalvars.h"

using namespace std;

// Routine to generate the case file
void EnsightGold::genCase(int nsteps, int saveEvery)
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

	// Time section (required for transient cases)
	fout << "TIME" << endl;
	fout << "time set: 1" << endl;
	fout << "number of steps: " << nsteps << endl;	
	fout << "filename start number: 00001" << endl;
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
void EnsightGold::genGeo(int r)
{
	// Open file stream
	char fileName[50];
	sprintf_s(fileName, (char*)"./Output/geom.geo");
	ofstream fout(fileName, ios::out|ios::binary);

	// Output format
	char buf[80] ;
	sprintf_s(buf, (char*)"C Binary");
	fout.write(buf, 80);
	
	// Description line 1
	sprintf_s(buf, (char*)"desc1" );
	fout.write(buf, 80 );
	
	// Description line 2
	sprintf_s(buf, (char*)"desc2");
	fout.write(buf, 80);
	
	// Tell it to assign node IDs for me
	sprintf_s(buf, (char*)"node id assign");
	fout.write(buf, 80);

	// Tell it to assign element IDs for me
	sprintf_s(buf, (char*)"element id assign");
	fout.write(buf, 80);
	
	// Necessary lines but do not apply to LBM data
	sprintf_s(buf, (char*)"part");
	fout.write(buf, 80);
	int partNum = 1;
	// Use reinterpret_cast to convert integrer to character
	fout.write(reinterpret_cast<char*>(&partNum), sizeof(int));
	sprintf_s(buf, (char*)"part desc");
	fout.write(buf, 80);

	// Dimensions of the grid
	sprintf_s(buf, (char*)"block uniform"); // Uniform grid
	fout.write(buf, 80);
	int ni = Grids[r].XPos.size()+1;
	int nj = Grids[r].YPos.size()+1;
#if (dims == 3)
	int nk = Grids[r].ZPos.size()+1;
#else
	int nk = 1;
#endif

	fout.write(reinterpret_cast<char*>(&ni), sizeof(int));
	fout.write(reinterpret_cast<char*>(&nj), sizeof(int));
	fout.write(reinterpret_cast<char*>(&nk), sizeof(int));


	// Origin location -- will be specific to a given grid
	float x_orig = 0.0, y_orig = 0.0, z_orig = 0.0;

	fout.write(reinterpret_cast<char*>(&x_orig), sizeof(float));
	fout.write(reinterpret_cast<char*>(&y_orig), sizeof(float));
	fout.write(reinterpret_cast<char*>(&z_orig), sizeof(float));

	float x_delta = (float)Grids[r].dx, y_delta = (float)Grids[r].dy, z_delta = (float)Grids[r].dz;

	fout.write(reinterpret_cast<char*>(&x_delta), sizeof(float));
	fout.write(reinterpret_cast<char*>(&y_delta), sizeof(float));
	fout.write(reinterpret_cast<char*>(&z_delta), sizeof(float));

	fout.close();

}


// Routine to generate the vectors file
void EnsightGold::genVec(int fileNum, int r)
{
	// Open file stream
	char fileName[50] ;
	sprintf_s(fileName, (char*)"./Output/velocity.%05d.vec", fileNum);

	ofstream fout;
	fout.open(fileName, ios::out|ios::binary);

	// Descriprion line
	char buf[80];
	sprintf_s(buf, (char*)"required description");
	fout.write(buf, 80);
	
	// Part header
	sprintf_s(buf, (char*)"part");
	fout.write(buf,80);
	int partNum = 1;
	fout.write(reinterpret_cast<char*>(&partNum), sizeof(int));
	sprintf_s(buf, (char*)"block");
	fout.write(buf, 80);

	// Data for part
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();
	int c;
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				int dir = 0;
				c = idxmap(i,j,k,dir,M_lim,K_lim,dims);
				float ux = (float)Grids[r].u[c];
				fout.write(reinterpret_cast<char*>(&ux), sizeof(float));
			}
		}
	}

	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				int dir = 1;
				c = idxmap(i,j,k,dir,M_lim,K_lim,dims);
				float uy = (float)Grids[r].u[c];
				fout.write(reinterpret_cast<char*>(&uy), sizeof(float));
			}
		}
	}

#if (dims == 3)
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				int dir = 2;
				c = idxmap(i,j,k,dir,M_lim,K_lim,dims);
				float uz = (float)Grids[r].u[c];
				fout.write(reinterpret_cast<char*>(&uz), sizeof(float));
			}
		}
	}
#else
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				float uz = 0.0;
				fout.write(reinterpret_cast<char*>(&uz), sizeof(float));
			}
		}
	}
#endif

	fout.close();

}


// Routine to generate the scalars file
void EnsightGold::genScal(int fileNum, int r)
{
	// File stream
	char fileName[50];
	sprintf_s(fileName, (char*)"./Output/rho.%05d.sca", fileNum);

	ofstream fout;
	fout.open(fileName, ios::out|ios::binary);

	// Description
	char buf[80];
	sprintf_s(buf, (char*)"required description");
	fout.write(buf, 80);
	
	// Part header
	sprintf_s(buf, (char*)"part");
	fout.write(buf, 80);
	int partNum = 1;
	fout.write(reinterpret_cast<char*>(&partNum), sizeof(int));
	sprintf_s(buf, (char*)"block");
	fout.write(buf, 80);

	// Data
	int N_lim = Grids[r].XPos.size();
	int M_lim = Grids[r].YPos.size();
	int K_lim = Grids[r].ZPos.size();
	int c;
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				c = idxmap(i,j,k,M_lim,K_lim);
				float rho = (float)Grids[r].rho[c];
				fout.write(reinterpret_cast<char*>(&rho), sizeof(float));
			}
		}
	}

	fout.close();

}