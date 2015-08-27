/*
	**************************************************************************************************************
	**************************************************************************************************************
	**																											**
	**											LatBo Defintions File											**
	**											 (For user editing)												**
	**																											**
	**************************************************************************************************************
	**************************************************************************************************************
*/


// Header guard
#ifndef LBM_DEFINITIONS_H
#define LBM_DEFINITIONS_H

// Declarations here
#include <time.h>			// Timing functionality
#include <iostream>			// IO functionality
#include <fstream>			// File functionality
#include <vector>			// Vector template access
#include <iomanip>			// Output precision control
#include <math.h>			// Mathematics
#include <string>			// String template access
#include <mpi.h>			// Enable MPI (MSMPI on Windows)
#include "generic_ops.h"	// Forward declarations of generic functions



/*	
***************************************************************************************************************
************************************** Global configuration data **********************************************
***************************************************************************************************************
*/

// Numbers
#define PI 3.14159265358979323846

// Using MPI?
#define BUILD_FOR_MPI

// Output Options
#define out_every 10			// How many timesteps before output
// Types of output
//#define TEXTOUT
#define VTK_WRITER
//#define MPI_VERBOSE

// Gravity (acts in +x direction)
//#define GRAVITY_ON
// Expression for the gravity force
#define grav_force ( 3 * vecnorm(u_0x,u_0y,u_0z) * nu / pow(abs(b_y - a_y),2) )
#define grav_direction 0	// Gravity direction (0 = x, 1 = y, 2 = z)

// Initialisation
//#define NO_FLOW			// Initialise the domain with no flow


/*	
***************************************************************************************************************
********************************************** Time data ******************************************************
***************************************************************************************************************
*/

#define T 10		// End time of simulation (if each time step increments by physical dt = dx)

/*	
***************************************************************************************************************
******************************************* Domain dimensions *************************************************
***************************************************************************************************************
*/

// MPI Data
#define Xcores 2
#define Ycores 2
#define Zcores 2	// This gets set to 1 later if doing a 2D problem

// Lattice properties (in lattice units)
#define dims 2		// Number of dimensions to the problem
#define N 300		// Number of x lattice sites
#define M 100		// Number of y lattice sites
#define K 100		// Number of z lattice sites

// Physical dimensions (dictates scaling)
#define a_x 0		// Start of domain-x
#define b_x 6.0		// End of domain-x
#define a_y 0		// Start of domain-y
#define b_y 2.0		// End of domain-y
#define a_z 0		// Start of domain-z
#define b_z 2.0		// End of domain-z


/*	
***************************************************************************************************************
*********************************************** Fluid data ****************************************************
***************************************************************************************************************
*/

// Data in lattice units
#define u_0x 0.05	// Initial x-velocity
#define u_0y 0		// Initial y-velocity
#define u_0z 0		// Initial z-velocity
#define rho_in 1	// Initial density
#define Re 50		// Desired Reynolds number
// nu computed based on above selections


/*	
***************************************************************************************************************
******************************************* Immersed Boundary *************************************************
***************************************************************************************************************
*/

// Master IBM switches //
//#define IBM_ON						// Turn on IBM
#define IBM_DEBUG						// Write IBM body and matrix data out to text files
#define IBBODY_TRACER					// Write out IBbody positions
#define LD_OUT							// Write out lift and drag (sum x and y forces on Lagrange markers of body)
#define STOP_EPSILON_RECOMPUTE			// Prevent recomputing of epsilon in an attempt to save time
#define CHEAP_NEAREST_NODE_DETECTION	// Perform a nearest-neighbour-type nearest node operation for IBM support calculation

// Switches for inserting certain bodies (enable only one at once!)
//#define INSERT_CIRCLE_SPHERE
//#define INSERT_RECTANGLE_CUBOID
//#define INSERT_BOTH
//#define INSERT_FILAMENT
//#define INSERT_FILARRAY
//#define _2D_RIGID_PLATE_IBM
//#define _2D_PLATE_WITH_FLAP
//#define _3D_RIGID_PLATE_IBM
//#define _3D_PLATE_WITH_FLAP

// Global properties
#define num_markers 10		// Number of Lagrange points (approximately)
#define ibb_deform false	// Default deformable property of body to be built

// Physical dimensions of rigid IB body or flexible plate
#define ibb_x 2.0		// x Position of body centre
#define ibb_y 1.0		// y Position of body centre
#define ibb_z 0.5		// z Position of body centre
#define ibb_w 0.5		// width (x) of IB body
#define ibb_l 0.25		// length (y) of IB body
#define ibb_d 0.25		// depth (z) of IB body
#define ibb_r .5		// radius of IB body

// Physical dimensions of flexible IB filament
#define ibb_length 0.5		// length of filament
#define ibb_start_x 2.50	// start x position of the filament
#define ibb_start_y 1.0		// start y position of the filament
#define ibb_start_z 0.5		// start z position of the filament

// Angles of filament or plate
#define ibb_angle_vert 20	// Inclination of filament in xy plane
#define ibb_angle_horz 0	// Inclination of filament in xz plane

// Boundary conditions of flexible filament or flexible plate
#define start_BC 2			// Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
#define end_BC 0			// Type of boundary condition at filament end:		0 == free; 1 = simply supported; 2 == clamped

// Mechanical properties of filament
#define ibb_delta_rho 1.5	// Difference in density (lattice units) between solid and fluid
#define ibb_EI .025			// Flexural rigidity (lattice units) of filament


/*	
***************************************************************************************************************
********************************************** Wall data ******************************************************
***************************************************************************************************************
*/

// Switches
#define SOLID_BLOCK_ON			// Turn on solid object (bounce-back) specified below
#define WALLS_ON				// Turn on no-slip walls (default is top, bottom, front, back unless WALLS_ON_2D is used)
//#define WALLS_ON_2D				// Limit no-slip walls to top and bottom no-slip walls
#define INLET_ON				// Turn on inlet boundary (assumed left-hand wall for now - default Zou-He)
#define INLET_DO_NOTHING		// Specify the inlet to be a do-nothing inlet condition
//#define INLET_REGULARISED		// Specify the inlet to be a regularised inlet condition (Latt & Chopard)
#define OUTLET_ON				// Turn on outlet boundary (assumed right-hand wall for now)

#ifdef SOLID_BLOCK_ON
// Wall labelling routine implements this
// Specified in lattice units (i.e. by index)
#define obj_x_min 80		// Index of start of object/wall in x-direction
#define obj_x_max 120		// Index of end of object/wall in x-direction
#define obj_y_min 40		// Index of start of object/wall in y-direction
#define obj_y_max 60		// Index of end of object/wall in y-direction
#define obj_z_min 40		// Index of start of object/wall in z-direction
#define obj_z_max 60		// Index of end of object/wall in z-direction
#endif


/*	
***************************************************************************************************************
******************************************** Multi-grid data **************************************************
***************************************************************************************************************
*/

#define NumLev 1		// Levels of refinement (can't use with IBM yet and won't span MPI ranks)
#define NumReg 2		// Number of refined regions (can be arbitrary if NumLev = 0)

#if NumLev != 0
// Global lattice indices for refined region on level L0 start numbering at 0

	#if NumReg == 2 // Inlcuded for testing purposes so I don't have to keep re-commenting bits
	static size_t RefXstart[NumReg]		= {40, 160};
	static size_t RefXend[NumReg]		= {60, 180};
	static size_t RefYstart[NumReg]		= {60, 20};
	static size_t RefYend[NumReg]		= {80, 40};
	// If doing 2D, these can be arbitrary values
	static size_t RefZstart[NumReg]		= {24, 24};
	static size_t RefZend[NumReg]		= {36, 36};

	#elif NumReg == 1
	static size_t RefXstart[NumReg]		= {64};
	static size_t RefXend[NumReg]		= {72};
	static size_t RefYstart[NumReg]		= {24};
	static size_t RefYend[NumReg]		= {36};
	static size_t RefZstart[NumReg]		= {24};
	static size_t RefZend[NumReg]		= {36};
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
	
	#define MPI_dir 26	// 3D MPI

#else
	#define nVels 9		// Use D2Q9

	// MPI config to 2D
	#undef Zcores
	#define Zcores 1
	#define MPI_dir 8
	
	// Set Z limits for 2D
	#undef a_z
	#define a_z 0

	#undef b_z
	#define b_z 2

	#undef K
	#define K 1

	// Set object limits for 2D
	#undef obj_z_min
	#define obj_z_min 0

	#undef obj_z_max
	#define obj_z_max 0

	#undef ibb_d
	#define ibb_d 0

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