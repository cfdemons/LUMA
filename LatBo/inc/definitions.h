/*
	**************************************************************************
	**************************************************************************
	**																		**
	**							LatBo Defintions File						**
	**							  (for user editing)						**
	**																		**
	**************************************************************************
	**************************************************************************
*/
// Definitions File Format ## 0.5-31 ## //


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
*******************************************************************************
**************************** Debugging Options ********************************
*******************************************************************************
*/


//#define MEGA_DEBUG				// Debug F, Feq, Macropcopic all in one file -- Warning: Heavy IO which kills performance
//#define DEBUG_STREAM				// Writes out the number and type of streaming operations used to test streaming exclusions
//#define MPI_VERBOSE				// Write out the buffers used by MPI plus more setup data
//#define IBM_DEBUG					// Write IBM body and matrix data out to text files
//#define IBBODY_TRACER				// Write out IBBody positions
//#define LD_OUT					// Write out lift and drag (sum x and y forces on Lagrange markers of IBBody)
//#define BFL_DEBUG					// Write out BFL marker positions and Q values out to files


/*
*******************************************************************************
************************* Global configuration data ***************************
*******************************************************************************
*/


// Numbers
#define PI 3.14159265358979323846

// Using MPI?
//#define BUILD_FOR_MPI

// Output Options
#define out_every 1			// How many timesteps before whole grid output
#define output_precision 16		// Precision of output


// Types of output
#define TEXTOUT
#define VTK_WRITER
//#define TECPLOT
//#define IO_LITE

// High frequency output options
//#define PROBE_OUTPUT
#define out_every_probe 250
const static int nProbes[3] = {3, 3, 3};		// Number of probes in each direction
// Start and End points for planes of probes
const static int xProbeLims[2] = {90, 270};
const static int yProbeLims[2] = {15, 45};
const static int zProbeLims[2] = {30, 120};


// Gravity
#define GRAVITY_ON
// Expression for the gravity force
#define grav_force 1e-10	//( 3 * gUtils.vecnorm(u_0x,u_0y,u_0z) * nu / pow(fabs(b_y - a_y),2) )
#define grav_direction 0	// Gravity direction (0 = x, 1 = y, 2 = z)

// Initialisation
//#define NO_FLOW			// Initialise the domain with no flow
//#define RESTARTING		// Initialise the GridObj with quantities read from a restart file
#define restart_out_every 500000

// LBM configuration
//#define USE_MRT

#if (dims == 3)
// MRT relaxation times (D3Q19) -- (see Stiebler 2011 paper for some improvements)
#define mrt_relax {1.0, 1.19, 1.4, 1.0, 1.2, 1.0, 1.2, 1.0, 1.2, omega, 1.4, omega, 1.4, omega, omega, omega, 1.98, 1.98, 1.98}
#else
// MRT relaxation times (D2Q9)
#define mrt_relax {1.0, 1.4, 1.4, 1.0, 1.2, 1.0, 1.2, omega, omega}
#endif

/*
*******************************************************************************
******************************** Time data ************************************
*******************************************************************************
*/

#define T 250	// Number of time steps


/*
*******************************************************************************
**************************** Domain Dimensions ********************************
*******************************************************************************
*/

// MPI Data
#define Xcores 2
#define Ycores 2
#define Zcores 2	// Set to 1 if doing a 2D problem when using custom MPI sizes

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
#define N 16		// Number of x lattice sites
#define M 16		// Number of y lattice sites
#define K 30		// Number of z lattice sites


// Physical dimensions (dictates scaling)
#define a_x 0		// Start of domain-x
#define b_x 0.61	// End of domain-x
#define a_y 0		// Start of domain-y
#define b_y 0.61	// End of domain-y
#define a_z 0		// Start of domain-z
#define b_z 8		// End of domain-z


/*
*******************************************************************************
******************************** Fluid Data ***********************************
*******************************************************************************
*/

// Fluid data in lattice units
//#define USE_INLET_PROFILE
#define u_ref 2.80583613916949e-05	// Reference velocity for scaling (mean inlet velocity)
#define u_max 0.06		// Max velocity of profile

// If not using an inlet profile, specify values or expressions here
#define u_0x 0			//u_ref //u_max*(1 - pow( ( (YPos[j] - ((b_y-a_y-dy)/2)) ) / ((b_y-a_y-dy)/2) ,2) )	// Initial x-velocity
#define u_0y 0			// Initial y-velocity
#define u_0z 0			// Initial z-velocity

#define rho_in 1		// Initial density
#define Re 1			// Desired Reynolds number

// nu computed based on above selections


/*
*******************************************************************************
****************************** Immersed Boundary ******************************
*******************************************************************************
*/

// Master IBM switches //
//#define IBM_ON						// Turn on IBM

//#define STOP_EPSILON_RECOMPUTE		// Prevent recomputing of epsilon in an attempt to save time
#define CHEAP_NEAREST_NODE_DETECTION	// Perform a nearest-neighbour-type nearest node operation for IBM support calculation

// Switches for inserting certain bodies (enable only one at once!)
//#define INSERT_CIRCLE_SPHERE
//#define INSERT_RECTANGLE_CUBOID
//#define INSERT_BOTH
#define INSERT_FILAMENT
//#define INSERT_FILARRAY
//#define _2D_RIGID_PLATE_IBM
//#define _2D_PLATE_WITH_FLAP
//#define _3D_RIGID_PLATE_IBM
//#define _3D_PLATE_WITH_FLAP

// Global properties
#define num_markers 19		// Number of Lagrange points (approximately)
#define ibb_deform false	// Default deformable property of body to be built

// Physical dimensions of rigid IB body or flexible plate
#define ibb_x 75.0		// x Position of body centre
#define ibb_y 75.0		// y Position of body centre
#define ibb_z 0.0		// z Position of body centre
#define ibb_w 10.0		// width (x) of IB body
#define ibb_l 10.0		// length (y) of IB body
#define ibb_d 0.0		// depth (z) of IB body
#define ibb_r 10.0		// radius of IB body

// Physical dimensions of flexible IB filament
#define ibb_length 0.2		// length of filament
#define ibb_start_x 0.3	// start x position of the filament
#define ibb_start_y 0.0	// start y position of the filament
#define ibb_start_z 0.0		// start z position of the filament

// Angles of filament or plate
#define ibb_angle_vert 90	// Inclination of filament in xy plane
#define ibb_angle_horz 0	// Inclination of filament in xz plane

// Boundary conditions of flexible filament or flexible plate
#define start_BC 2			// Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
#define end_BC 0			// Type of boundary condition at filament end:		0 == free; 1 = simply supported; 2 == clamped

// Mechanical properties of filament
#define ibb_delta_rho 1.0	// Difference in density (lattice units) between solid and fluid
#define ibb_EI 2.0			// Flexural rigidity (lattice units) of filament


/*
*******************************************************************************
********************************** Wall Data **********************************
*******************************************************************************
*/

// Inlets
//#define INLET_ON				// Turn on inlet boundary (assumed left-hand wall for now - default Zou-He)
//#define INLET_DO_NOTHING		// Specify the inlet to be a do-nothing inlet condition (overrides other options)
//#define INLET_REGULARISED		// Specify the inlet to be a regularised inlet condition (Latt & Chopard)
//#define UNIFORM_INLET			// Make the inlet a uniform inlet
//#define INLET_NRBC				// Turn on NRBC at inlet


// Outlets
//#define OUTLET_ON				// Turn on outlet boundary (assumed right-hand wall for now)
//#define OUTLET_NRBC				// Turn on NRBC at outlet


// Periodicity
#define	PERIODIC_BOUNDARIES

// Solids
#define WALLS_ON				// Turn on no-slip walls (default is top, bottom, front, back unless WALLS_ON_2D is used)
#define WALLS_ON_2D				// Limit no-slip walls to top and bottom no-slip walls only
//#define WALLS_ON_FLOOR_ONLY		// Limit no-slip walls to bottom no-slip wall only
#define wall_thickness	1		// Thickness of walls in coarsest lattice units



/*
*******************************************************************************
********************************* Object Data *********************************
*******************************************************************************
*/

// Bounce-back solids
//#define SOLID_BLOCK_ON			// Turn on solid object (bounce-back) specified below

#ifdef SOLID_BLOCK_ON
	#define block_on_grid_lev 0		// Provide grid level on which block should be added 
	#define block_on_grid_reg 0		// Provide grid region on which block should be added 
	// Wall labelling routine implements this
	// Specified in lattice units (i.e. by index) local to the chosen grid level
	#define obj_x_min 0		// Index of start of object/wall in x-direction
	#define obj_x_max 60		// Index of end of object/wall in x-direction
	#define obj_y_min 0		// Index of start of object/wall in y-direction
	#define obj_y_max 4		// Index of end of object/wall in y-direction
	#define obj_z_min 105		// Index of start of object/wall in z-direction
	#define obj_z_max 135		// Index of end of object/wall in z-direction
#endif


// Bounce-back objects from point clouds
//#define SOLID_FROM_FILE

#ifdef SOLID_FROM_FILE
	#define object_on_grid_lev 0		// Provide grid level on which object should be added 
	#define object_on_grid_reg 0		// Provide grid region on which object should be added
	// Following specified in lattice units (i.e. by index) local to the chosen grid level
	#define start_object_x 100
	#define start_object_y 1
	#define start_object_z 10
	#define object_length_x 60			// The object input is scaled based on this dimension
#endif



// BFL objects
//#define BFL_ON

#ifdef BFL_ON
	#define bfl_on_grid_lev 1		// Provide grid level on which BFL body should be added 
	#define bfl_on_grid_reg 0		// Provide grid region on which BFL body should be added
	// Following specified in lattice units (i.e. by index) local to the chosen grid level
	#define start_bfl_x 20
	#define start_bfl_y 10
	#define start_bfl_z 10
	#define bfl_length_x 40		// The BFL object input is scaled based on this dimension
#endif



/*
*******************************************************************************
****************************** Multi-grid Data ********************************
*******************************************************************************
*/

#define NumLev 0		// Levels of refinement (can't use with IBM yet)
#define NumReg 1		// Number of refined regions (can be arbitrary if NumLev = 0)

#if NumLev != 0
// Global lattice indices (in terms of each grid level) for each refined region specified on each level


// Following options are only here to making testing different grid combinations easier
#if (NumReg == 2 && NumLev == 2) 
	const static size_t RefXstart[NumLev][NumReg]	= { {5, 5}, {2, 2} };
	const static size_t RefXend[NumLev][NumReg]		= { {25, 25}, {20, 10} };
	const static size_t RefYstart[NumLev][NumReg]	= { {5, 14}, {5, 2} };
	const static size_t RefYend[NumLev][NumReg]		= { {12, 25}, {10, 10} };
	// If doing 2D, these can be arbitrary values
	static size_t RefZstart[NumLev][NumReg]		= { {5, 10}, {2, 2} };
	static size_t RefZend[NumLev][NumReg]		= { {20, 15}, {10, 10} };

#elif (NumReg == 1 && NumLev == 1)
	const static size_t RefXstart[NumLev][NumReg]	= { 10 };
	const static size_t RefXend[NumLev][NumReg]		= { 80 };
	const static size_t RefYstart[NumLev][NumReg]	= { 0 };
	const static size_t RefYend[NumLev][NumReg]		= { 30 };
	// If doing 2D, these can be arbitrary values
	static size_t RefZstart[NumLev][NumReg]		= { 5 };
	static size_t RefZend[NumLev][NumReg]		= { 35 };

#elif (NumReg == 1 && NumLev == 2)
	const static size_t RefXstart[NumLev][NumReg]	= { {226}, {12} };
	const static size_t RefXend[NumLev][NumReg]		= { {382}, {228} };
	const static size_t RefYstart[NumLev][NumReg]	= { {0}, {0} };
	const static size_t RefYend[NumLev][NumReg]		= { {48}, {84} };
	// If doing 2D, these can be arbitrary values
	static size_t RefZstart[NumLev][NumReg]		= { {78}, {12} };
	static size_t RefZend[NumLev][NumReg]		= { {162}, {156} };


#endif

#endif


/*
*******************************************************************************
************************* Clean-up: NOT FOR EDITING ***************************
*******************************************************************************
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

	// Set BFL start for 2D
	#undef start_object_z
	#define start_object_z 0

	// Set Object start for 2D
	#undef start_bfl_z
	#define start_bfl_z 0

#endif

#if NumLev == 0
	// Set region info to default as no refinement
	const static size_t RefXstart[1][1]		= {0};
	const static size_t RefXend[1][1]		= {0};
	const static size_t RefYstart[1][1]		= {0};
	const static size_t RefYend[1][1]		= {0};
	static size_t RefZstart[1][1]			= {0};
	static size_t RefZend[1][1]				= {0};
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
