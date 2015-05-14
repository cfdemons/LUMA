// Header guard
#ifndef LBM_DEFINITIONS_H
#define LBM_DEFINITIONS_H

// Declarations here
#include <time.h>		// Timing functionality
#include <iostream>		// IO functionality
#include <fstream>		// File functionality
#include <vector>		// Vector template access
#include <iomanip>		// Output precision control
#include <math.h>		// Mathematics
#include <string>		// String template access

#include "ops_generic.h"
#include "ops_mapping.h"

/*	
***************************************************************************************************************
************************************** Global configuration data **********************************************
***************************************************************************************************************
*/
#define PI 3.14159265358979323846
//#define TEXTOUT
#define ENSIGHTGOLD

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
#define dims 3	// Number of dimensions to the problem
#define N 32	// Number of x lattice sites
#define M 32	// Number of y lattice sites
#define K 32	// Number of z lattice sites
#define a_x 0	// Start of domain-x
#define b_x 32	// End of domain-x
#define a_y 0	// Start of domain-y
#define b_y 32	// End of domain-y
#define a_z 0	// Start of domain-z
#define b_z 32	// End of domain-z

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
#define NumLev 1		// Levels of refinement
#define NumReg 1		// Number of refined regions (can be arbitrary if NumLev = 0)

#if NumLev != 0
// Lattice indices for refined region on level L0 start numbering at 0

	#if NumReg == 2 // Inlcuded for testing purposes so I don't have to keep re-commenting bits
	static size_t RefXstart[NumReg]		= {1, 5};
	static size_t RefXend[NumReg]		= {4, 8};
	static size_t RefYstart[NumReg]		= {1, 5};
	static size_t RefYend[NumReg]		= {4, 8};
	// If doing 2D, these can be arbitrary values
	static size_t RefZstart[NumReg]		= {1, 1};
	static size_t RefZend[NumReg]		= {4, 4};

	#elif NumReg == 1
	static size_t RefXstart[NumReg]		= {1};
	static size_t RefXend[NumReg]		= {4};
	static size_t RefYstart[NumReg]		= {1};
	static size_t RefYend[NumReg]		= {4};
	static size_t RefZstart[NumReg]		= {1};
	static size_t RefZend[NumReg]		= {4};
	#endif

#endif

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
	
	// Set Z limits for 2D
	#undef a_z
	#define a_z 0

	#undef b_z
	#define b_z 2

	#undef K
	#define K 1
#endif

#if NumLev == 0
// Set region info to default as no refinement
#undef NumReg
#define NumReg 1
static size_t RefXstart[NumReg]		= {0};
static size_t RefXend[NumReg]		= {0};
static size_t RefYstart[NumReg]		= {0};
static size_t RefYend[NumReg]		= {0};
static size_t RefZstart[NumReg]		= {0};
static size_t RefZend[NumReg]		= {0};
#endif


#endif