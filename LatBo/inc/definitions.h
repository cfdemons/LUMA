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
#include <mpi.h>			// Enable MPI



/*
***************************************************************************************************************
************************************** Global configuration data **********************************************
***************************************************************************************************************
*/

// Numbers
#define PI 3.14159265358979323846

// Using MPI?
//#define BUILD_FOR_MPI

// Output Options
#define out_every 10			// How many timesteps before whole grid output

// Types of output
//#define TEXTOUT
#define VTK_WRITER
//#define TECPLOT
//#define MPI_VERBOSE

// High frequency output options
//#define PROBE_OUTPUT
#define out_every_probe 250
const static unsigned int nProbes[3] = {3, 3, 3};		// Number of probes in each direction
// Start and End points for planes of probes
const static int xProbeLims[2] = {90, 270};
const static int yProbeLims[2] = {15, 45};
const static int zProbeLims[2] = {30, 120};


// Gravity
//#define GRAVITY_ON
// Expression for the gravity force
#define grav_force ( 3 * gUtils.vecnorm(u_0x,u_0y,u_0z) * nu / pow(fabs(b_y - a_y),2) )
#define grav_direction 0	// Gravity direction (0 = x, 1 = y, 2 = z)

// Initialisation
#define NO_FLOW			// Initialise the domain with no flow
//#define RESTARTING		// Initialise the GridObj with quantities read from a restart file
#define restart_out_every 500000

// LBM configuration
//#define USE_MRT

#if (dims == 3)
// MRT relaxation times (D3Q19)
#define mrt_relax {1.0, 1.19, 1.4, 1.0, 1.2, 1.0, 1.2, 1.0, 1.2, omega, 1.4, omega, 1.4, omega, omega, omega, 1.98, 1.98, 1.98}
#else
// MRT relaxation times (D2Q9)
#define mrt_relax {1.0, 1.4, 1.4, 1.0, 1.2, 1.0, 1.2, omega, omega}
#endif

/*
***************************************************************************************************************
********************************************** Time data ******************************************************
***************************************************************************************************************
*/

#define T 10000		// Number of time steps


/*
***************************************************************************************************************
******************************************* Domain dimensions *************************************************
***************************************************************************************************************
*/

// MPI Data
#define Xcores 1
#define Ycores 1
#define Zcores 1	// Set to 1 if doing a 2D problem when using custom MPI sizes

//#define USE_CUSTOM_MPI_SIZES

// MPI local grid sizes (Cartesian topolgy numbered in z, y then x directions)
#ifdef USE_CUSTOM_MPI_SIZES
const static size_t xRankSize[Xcores*Ycores*Zcores]		= {50, 50, 50, 50, 350, 350, 350, 350};
const static size_t yRankSize[Xcores*Ycores*Zcores]		= {20, 20, 130, 130, 20, 20, 130, 130};
// The following can be arbitrary if doing a 2D problem
const static size_t zRankSize[Xcores*Ycores*Zcores]		= {20, 30, 20, 30, 20, 30, 20, 30};
#endif


// Lattice properties (in lattice units)
#define dims 2		// Number of dimensions to the problem
#define N 150		// Number of x lattice sites
#define M 150		// Number of y lattice sites
#define K 30		// Number of z lattice sites


// Physical dimensions (dictates scaling)
#define a_x 0		// Start of domain-x
#define b_x 150.0	// End of domain-x
#define a_y 0		// Start of domain-y
#define b_y 150.0	// End of domain-y
#define a_z 0		// Start of domain-z
#define b_z 8		// End of domain-z


/*
***************************************************************************************************************
*********************************************** Fluid data ****************************************************
***************************************************************************************************************
*/

// Fluid data in lattice units
#define USE_INLET_PROFILE
#define u_ref 0.04		// Reference velocity for scaling (mean inlet velocity)
#define u_max 0.06		// Max velocity of profile

// If not using an inlet profile, specify values or expressions here
#define u_0x u_max*(1 - pow( ( (YPos[j] - ((b_y-a_y-dy)/2)) ) / ((b_y-a_y-dy)/2) ,2) )	// Initial x-velocity
#define u_0y 0			// Initial y-velocity
#define u_0z 0			// Initial z-velocity

#define rho_in 1		// Initial density
#define Re 10			// Desired Reynolds number

// nu computed based on above selections


/*
***************************************************************************************************************
******************************************* Immersed Boundary *************************************************
***************************************************************************************************************
*/

// Master IBM switches //
#define IBM_ON						// Turn on IBM
//#define IBM_DEBUG						// Write IBM body and matrix data out to text files
//#define IBBODY_TRACER					// Write out IBbody positions
//#define LD_OUT						// Write out lift and drag (sum x and y forces on Lagrange markers of body)
//#define STOP_EPSILON_RECOMPUTE			// Prevent recomputing of epsilon in an attempt to save time
#define CHEAP_NEAREST_NODE_DETECTION	// Perform a nearest-neighbour-type nearest node operation for IBM support calculation

// Switches for inserting certain bodies (enable only one at once!)
//#define INSERT_CIRCLE_SPHERE
//#define INSERT_RECTANGLE_CUBOID
//#define INSERT_BOTH
#define INSERT_FILAMENT
//#define INSERT_FILARRAY
#define _2D_RIGID_PLATE_IBM
//#define _2D_PLATE_WITH_FLAP
//#define _3D_RIGID_PLATE_IBM
//#define _3D_PLATE_WITH_FLAP

// Global properties
#define num_markers 31		// Number of Lagrange points (approximately)
#define ibb_deform true	// Default deformable property of body to be built

// Physical dimensions of rigid IB body or flexible plate
#define ibb_x 50.0		// x Position of body centre
#define ibb_y 10.0		// y Position of body centre
#define ibb_z 0.0		// z Position of body centre
#define ibb_w 10.0		// width (x) of IB body
#define ibb_l 10.0		// length (y) of IB body
#define ibb_d 0.0		// depth (z) of IB body
#define ibb_r 5.0		// radius of IB body

// Physical dimensions of flexible IB filament
#define ibb_length 30.0		// length of filament
#define ibb_start_x 75.0	// start x position of the filament
#define ibb_start_y 60.0	// start y position of the filament
#define ibb_start_z 0.0		// start z position of the filament

// Angles of filament or plate
#define ibb_angle_vert 90	// Inclination of filament in xy plane
#define ibb_angle_horz 0	// Inclination of filament in xz plane

// Boundary conditions of flexible filament or flexible plate
#define start_BC 2			// Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
#define end_BC 0			// Type of boundary condition at filament end:		0 == free; 1 = simply supported; 2 == clamped

// Mechanical properties of filament
#define ibb_delta_rho 1.5	// Difference in density (lattice units) between solid and fluid
#define ibb_EI 0.001			// Flexural rigidity (lattice units) of filament


/*
***************************************************************************************************************
********************************************** Wall data ******************************************************
***************************************************************************************************************
*/

// Switches
//#define SOLID_BLOCK_ON			// Turn on solid object (bounce-back) specified below
#define WALLS_ON				// Turn on no-slip walls (default is top, bottom, front, back unless WALLS_ON_2D is used)
#define WALLS_ON_2D				// Limit no-slip walls to top and bottom no-slip walls only
#define INLET_ON				// Turn on inlet boundary (assumed left-hand wall for now - default Zou-He)
//#define INLET_DO_NOTHING		// Specify the inlet to be a do-nothing inlet condition (overrides other options)
#define INLET_REGULARISED		// Specify the inlet to be a regularised inlet condition (Latt & Chopard)
//#define UNIFORM_INLET			// Make the inlet a uniform inlet
#define OUTLET_ON				// Turn on outlet boundary (assumed right-hand wall for now)
//#define PERIODIC_BOUNDARIES		// Turn on periodic boundary conditions (only applies to fluid-fluid interfaces)

#ifdef SOLID_BLOCK_ON
// Wall labelling routine implements this
// Specified in lattice units (i.e. by index)
#define obj_x_min 90		// Index of start of object/wall in x-direction
#define obj_x_max 120		// Index of end of object/wall in x-direction
#define obj_y_min 1			// Index of start of object/wall in y-direction
#define obj_y_max 31		// Index of end of object/wall in y-direction
#define obj_z_min 105		// Index of start of object/wall in z-direction
#define obj_z_max 135		// Index of end of object/wall in z-direction
#endif


/*
***************************************************************************************************************
******************************************** Multi-grid data **************************************************
***************************************************************************************************************
*/

#define NumLev 0		// Levels of refinement (can't use with IBM yet and won't span MPI ranks)
#define NumReg 1		// Number of refined regions (can be arbitrary if NumLev = 0)

#if NumLev != 0
// Global lattice indices for refined region on level L0 start numbering at 0

	#if NumReg == 2 // Makes it easier to set up multi-region cases
	const static size_t RefXstart[NumReg]	= {10, 50};
	const static size_t RefXend[NumReg]		= {30, 70};
	const static size_t RefYstart[NumReg]	= {10, 40};
	const static size_t RefYend[NumReg]		= {20, 50};
	// If doing 2D, these can be arbitrary values
	static size_t RefZstart[NumReg]		= {5, 25};
	static size_t RefZend[NumReg]		= {15, 35};

	#elif NumReg == 1
	const static size_t RefXstart[NumReg]	= {10};
	const static size_t RefXend[NumReg]		= {30};
	const static size_t RefYstart[NumReg]	= {5};
	const static size_t RefYend[NumReg]		= {15};
	static size_t RefZstart[NumReg]			= {5};
	static size_t RefZend[NumReg]			= {15};
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
const static size_t RefXstart[NumReg]		= {0};
const static size_t RefXend[NumReg]			= {0};
const static size_t RefYstart[NumReg]		= {0};
const static size_t RefYend[NumReg]			= {0};
static size_t RefZstart[NumReg]				= {0};
static size_t RefZend[NumReg]				= {0};
#endif

// Clean up for using profiled inlet
#ifdef USE_INLET_PROFILE
#undef u_0x
#define u_0x ux_in[j]
#undef u_0y
#define u_0y uy_in[j]
#undef u_0z
#define u_0z uz_in[j]
#endif

#endif
