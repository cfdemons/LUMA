// Header guard
#ifndef GLOBAL_H
#define GLOBAL_H

// Include definitions
#include "LBM_definitions.h"
 
// Define GridData structure
struct GridData {

	// 1D arrays
	std::vector<int> XInd;
	std::vector<int> YInd;
	std::vector<int> ZInd;
	std::vector<double> XPos;
	std::vector<double> YPos;
	std::vector<double> ZPos;

	// Flattened 4D arrays (i,j,k,vel)
	std::vector<double> f;
	std::vector<double> feq;
	std::vector<double> u;

	// Flattened 3D arrays (i,j,k)
	std::vector<double> rho;
	std::vector<int> LatTyp;

	// Scalars
	double omega;
	double dx;
	double dy;
	double dz;
	double dt;

};


// Global variable references
extern const int c[3][nVels];
extern const double w[nVels];
extern const double cs;
extern unsigned int PartOpt;
extern unsigned int GridOut;
extern GridData Grids[Nref+1];
 
#endif