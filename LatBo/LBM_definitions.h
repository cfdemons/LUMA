// Header guard
#ifndef HEADERLBM_H
#define HEADERLBM_H

// Declarations here
#include <time.h>		// Timing functionality
#include <iostream>		// IO functionality
#include <fstream>		// File functionality
#include <vector>		// Vector template access
#include <iomanip>		// Output precision control
#include <math.h>		// Mathematics
#include <string>		// String template access

// Forward declaration of functions

// Generic functions
std::vector<int> onespace(int min, int max);						// Function: onespace
std::vector<double> linspace(double min, double max, int n);		// Function: linspace
double vecnorm(double vec[2]);										// Function: vecnorm + 3 overloads
double vecnorm(double val1, double val2);
double vecnorm(double val1, double val2, double val3);
double vecnorm(double vec[3]);
std::string int2str(int number);									// Function: int2str
void lbm_write3(int r, int t);										// Function: lbm_write

// Mapping functions
int idxmap (int i, int j, int k, int vel, int M, int K, int nVels);	// Function: idxmap + 2 overloads
int idxmap (int i, int j, int vel, int M, int nVels);
int idxmap (int i, int j, int M);
double posmapref (double coarse_pos, int fine_level, char direction, char plusminus);
																	// Function: posmapref
std::vector<int> indmapref(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start);
																	// Function: indmapref

// Initialisation functions
void LBM_init_vel(int r);		// Function: LBM_init_vel
void LBM_init_rho(int r);		// Function: LBM_init_rho
void LBM_init_multi();			// Function: LBM_init_multi

// LBM operations
void LBM_multi(int r);							// Function: LBM_multi
void LBM_collide(int r, bool core_flag);		// Function: LBM_collide + 1 overload
double LBM_collide(int i, int j, int k, int v, int r);
void LBM_stream(int r);							// Function: LBM_stream
void LBM_macro(int r);							// Function: LBM_macro
void LBM_boundary (int r, int bc_type_flag);	// Function: LBM_boundary

// Multi-grid operations
void LBM_explode(int r);			// Function: LBM_explode
void LBM_coalesce(int r);			// Function: LBM_coalesce

/*	
***************************************************************************************************************
************************************** Global configuration data **********************************************
***************************************************************************************************************
*/
#define PI 3.14159265358979323846

/*	
***************************************************************************************************************
********************************************** Time data ******************************************************
***************************************************************************************************************
*/
#define T 50		// End time of simulation
#define deltat 1	// Time step size

/*	
***************************************************************************************************************
******************************************* Domain dimensions *************************************************
***************************************************************************************************************
*/
#define dims 2	// Number of dimensions to the problem
#define N 6	// Number of x lattice sites
#define M 7	// Number of y lattice sites
#define K 8	// Number of z lattice sites
#define a_x 0	// Start of domain-x
#define b_x 6	// End of domain-x
#define a_y 0	// Start of domain-y
#define b_y 7	// End of domain-y
#define a_z 0	// Start of domain-z
#define b_z 8	// End of domain-z

/*	
***************************************************************************************************************
*********************************************** Fluid data ****************************************************
***************************************************************************************************************
*/
#define u_0x .2		// Initial x-velocity
#define u_0y 0		// Initial y-velocity
#define u_0z 0		// Initial z-velocity
#define rho_in 1	// Initial density
#define nu .02		// Kinematic viscosity
#define kn 1		// Vortices in x direction domain /2
#define km 1		// Vortices in y direction domain /2
#define kk 1		// Vortices in z direction domain /2

/*	
***************************************************************************************************************
******************************************** Multi-grid data **************************************************
***************************************************************************************************************
*/
#define Nref 1		// Levels of refinement
// Lattice indices for refined region on level L0 start numbering at 0
#define Ref_startX 1
#define Ref_endX 4
#define Ref_startY 1
#define Ref_endY 5
#define Ref_startZ 1
#define Ref_endZ 6

/*	
***************************************************************************************************************
************************************** Clean-up -- no need to edit ********************************************
***************************************************************************************************************
*/

// Set default options if using 2D
#if dims == 3
	#define nVels 19	// Use D3Q19
#else
	#define nVels 9		// Use D2Q9
	
	#undef a_z
	#define a_z 0

	#undef b_z
	#define b_z 2

	#undef K
	#define K 1

	#undef Ref_startZ
	#define Ref_startZ 1

	#undef Ref_endZ
	#define Ref_endZ 1
#endif




#endif