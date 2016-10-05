/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) 2015, 2016
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * distribution without written consent.
 *
 */

/// LUMA version
#define LUMA_VERSION "1.2.0-alpha"


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


//#define L_MEGA_DEBUG				///< Debug F, Feq, Macroscopic all in one file -- Warning: Heavy IO which kills performance
//#define L_INC_RECV_LAYER			///< Flag to include writing out receiver layer sites in MPI builds
//#define L_DEBUG_STREAM			///< Writes out the number and type of streaming operations used to test streaming exclusions
//#define L_MPI_VERBOSE				///< Write out the buffers used by MPI plus more setup data
//#define L_IBM_DEBUG				///< Write IBM body and matrix data out to text files
//#define L_IBBODY_TRACER			///< Write out IBBody positions
//#define L_BFL_DEBUG				///< Write out BFL marker positions and Q values out to files
//#define L_CLOUD_DEBUG				///< Write out to a file the cloud that has been read in
//#define L_LOG_TIMINGS				///< Write out the initialisation, time step and mpi timings to an output file
//#define L_HDF_DEBUG				///< Write some HDF5 debugging information
//#define L_TEXTOUT					///< Verbose ASCII output of grid information


/*
*******************************************************************************
************************* Global configuration data ***************************
*******************************************************************************
*/


// Numbers
#define L_PI 3.14159265358979323846		///< PI definition

// Using MPI?
#define L_BUILD_FOR_MPI				///< Enable MPI features in build

// Output Options
#define L_out_every 100				///< How many timesteps before whole grid output
#define L_out_every_forces 100		///< Specific output frequency of body forces
#define L_output_precision 8		///< Precision of output (for text writers)

// Types of output
#define L_IO_LITE					///< ASCII dump on output
#define L_HDF5_OUTPUT				///< HDF5 dump on output
//#define L_LD_OUT					///< Write out lift and drag (all bodies)

// High frequency output options
//#define L_PROBE_OUTPUT						///< Turn on probe output
#define L_out_every_probe 250					///< Write out frequency of probe output
const static int nProbes[3] = {3, 3, 3};		///< Number of probes in each direction (x, y, z)
const static int xProbeLims[2] = {90, 270};		///< Limits of X plane for array of probes
const static int yProbeLims[2] = {15, 45};		///< Limits of Y plane for array of probes
const static int zProbeLims[2] = {30, 120};		///< Limits of Z plane for array of probes


// Gravity
//#define GRAVITY_ON						///< Turn on gravity force
/// Expression for the gravity force
#define L_grav_force 0.0001
#define L_grav_direction eXDirection		///< Gravity direction (specify using enumeration)

// Initialisation
//#define L_NO_FLOW							///< Initialise the domain with no flow
#define L_RESTARTING						///< Initialise the GridObj with quantities read from a restart file
#define L_restart_out_every 1000			///< Frequency of write out of restart file

// LBM configuration
#define L_USE_KBC_COLLISION					///< Use KBC collision operator instead of LBGK by default


/*
*******************************************************************************
******************************** Time data ************************************
*******************************************************************************
*/

#define L_Timesteps 100		///< Number of time steps to run simulation for


/*
*******************************************************************************
**************************** Domain Dimensions ********************************
*******************************************************************************
*/

// MPI Data
#define L_Xcores 4		///< Number of MPI ranks to divide domain into in X direction
#define L_Ycores 2		///< Number of MPI ranks to divide domain into in Y direction
/// Number of MPI ranks to divide domain into in Z direction.
/// Set to 1 if doing a 2D problem when using custom MPI sizes
#define L_Zcores 2

//#define L_USE_CUSTOM_MPI_SIZES		///< Define to use custom decomposition otherwise decomposition will be uniform

// MPI local grid sizes (Cartesian topolgy numbered in z, y then x directions)
#ifdef L_USE_CUSTOM_MPI_SIZES
	/// Number of sites in X direction for each custom rank
	const static size_t xRankSize[L_Xcores*L_Ycores*L_Zcores]		= {50, 50, 50, 50, 350, 350, 350, 350};
	/// Number of sites in Y direction for each custom rank
	const static size_t yRankSize[L_Xcores*L_Ycores*L_Zcores]		= {20, 20, 130, 130, 20, 20, 130, 130};
	/// Number of sites in Z direction for each custom rank.
	/// The following can be arbitrary if doing a 2D problem
	const static size_t zRankSize[L_Xcores*L_Ycores*L_Zcores]		= {20, 30, 20, 30, 20, 30, 20, 30};
#endif


// Lattice properties (in lattice units)
#define L_dims 2		///< Number of dimensions to the problem
#define L_N 100			///< Number of x lattice sites
#define L_M 100			///< Number of y lattice sites
#define L_K 100			///< Number of z lattice sites


// Physical dimensions (dictates scaling)
#define L_a_x 0			///< Start of domain-x
#define L_b_x 1		///< End of domain-x
#define L_a_y 0			///< Start of domain-y
#define L_b_y 1		///< End of domain-y
#define L_a_z 0			///< Start of domain-z
#define L_b_z 1		///< End of domain-z


/*
*******************************************************************************
******************************** Fluid Data ***********************************
*******************************************************************************
*/

// Fluid data in lattice units
//#define L_USE_INLET_PROFILE		///< Use an inlet profile
//#define L_PARABOLIC_INLET		///< Use analytic expression for inlet profile - if not then ASCII file is read (requires L_USE_INLET_PROFILE)
#define L_u_ref 0.04			///< Reference velocity for scaling, can be mean inelt velocity
#define L_u_max L_u_ref*1.5		///< Max velocity of inlet profile

// If not using an inlet profile, specify values or expressions here
#define L_u_0x L_u_ref		///< Initial/inlet x-velocity
#define L_u_0y 0			///< Initial/inlet y-velocity
#define L_u_0z 0			///< Initial/inlet z-velocity

#define L_rho_in 1			///< Initial density
#define L_Re 100			///< Desired Reynolds number

// nu computed based on above selections


/*
*******************************************************************************
****************************** Immersed Boundary ******************************
*******************************************************************************
*/

// Master IBM switches //
//#define L_IBM_ON						///< Turn on IBM
#define L_IB_Lev 0					///< Grid level for immersed boundary object (0 if no refined regions, -1 if no IBM)
#define L_IB_Reg 0					///< Grid region for immersed boundary object (0 if no refined regions, -1 if no IBM)

//#define L_STOP_EPSILON_RECOMPUTE		///< Prevent recomputing of epsilon in an attempt to save time
#define L_VTK_BODY_WRITE				///< Write out the bodies to a VTK file

// Read in IB Body from File
//#define L_IBB_FROM_FILE			///< Build immersed bodies from a point cloud file

	#define L_ibb_on_grid_lev L_IB_Lev		///< Provide grid level on which object should be added
	#define L_ibb_on_grid_reg L_IB_Reg		///< Provide grid region on which object should be added
	// Following specified in physical distances
	#define L_start_ibb_x 0.3		///< Start X of object bounding box
	#define L_start_ibb_y 0.2		///< Start Y of object bounding box
	#define L_centre_ibb_z 0.5		///< Centre of object bounding box in Z direction
	#define L_ibb_length 0.5		///< The object input is scaled based on this dimension
	#define L_ibb_scale_direction eXDirection	///< Scale in this direction (specify as enumeration)
	#define L_ibb_length_ref 0.5	///< Reference length to be used in the definition of Reynolds number

// Default global properties
#define L_num_markers 120		///< Number of Lagrange points to use when building a prefab body (approximately)
#define L_ibb_deform false		///< Default deformable property of body to be built (whether it moves or not)
#define L_ibb_flex_rigid false	///< Whether a structural calculation needs to be performed on the body


// Switches for inserting certain bodies (enable only one at once!)
//#define L_INSERT_CIRCLE_SPHERE
//#define L_INSERT_RECTANGLE_CUBOID
//#define L_INSERT_BOTH
#define L_INSERT_FILAMENT
//#define L_INSERT_FILARRAY
//#define L_2D_RIGID_PLATE_IBM
//#define L_2D_PLATE_WITH_FLAP
//#define L_3D_RIGID_PLATE_IBM
//#define L_3D_PLATE_WITH_FLAP

// Physical dimensions of rigid IB body or flexible plate
#define L_ibb_x 0.2		///< X Position of body centre
#define L_ibb_y 0.5		///< Y Position of body centre
#define L_ibb_z 0.0		///< Z Position of body centre
#define L_ibb_w 0.5		///< Width (x) of IB body
#define L_ibb_l 0.5		///< Length (y) of IB body
#define L_ibb_d 0.5		///< Depth (z) of IB body
#define L_ibb_r 0.25	///< Radius of IB body

// Physical dimensions of flexible IB filament
#define L_ibb_filament_length 0.5		///< Length of filament
#define L_ibb_filament_start_x 0.2		///< Start X position of the filament
#define L_ibb_filament_start_y 0.5		///< Start Y position of the filament
#define L_ibb_filament_start_z 0.5		///< Start Z position of the filament

// Angles of filament or plate
#define L_ibb_angle_vert 90		///< Inclination of filament in XY plane
#define L_ibb_angle_horz 0		///< Inclination of filament in XZ plane

// Boundary conditions of flexible filament or flexible plate
#define L_start_BC 2		///< Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
#define L_end_BC 0			///< Type of boundary condition at filament end:	0 == free; 1 = simply supported; 2 == clamped

// Mechanical properties of filament
#define L_ibb_delta_rho 1.0		///< Difference in density (lattice units) between solid and fluid
#define L_ibb_EI 2.0			///< Flexural rigidity (lattice units) of filament


/*
*******************************************************************************
********************************** Wall Data **********************************
*******************************************************************************
*/

// Virtual Wind Tunnels
//#define L_UPSTREAM_TUNNEL			///< Adds an inlet to all faces except exit
//#define L_FREESTREAM_TUNNEL			///< Adds a inlet to all faces


// Inlets
#define L_INLET_ON				///< Turn on inlet boundary (assumed left-hand wall - default Do Nothing)
//#define L_INLET_REGULARISED	///< Specify the inlet to be a regularised inlet condition (Latt & Chopard)
//#define L_INLET_NRBC			///< Turn on NRBC at inlet


// Outlets
#define L_OUTLET_ON				///< Turn on outlet boundary (assumed right-hand wall -- default First Order Extrap.)
//#define L_OUTLET_NRBC			///< Turn on NRBC at outlet


// Periodicity
//#define L_PERIODIC_BOUNDARIES		///< Turn on periodic boundary conditions (only applies to fluid-fluid interfaces)


// Solids
#define L_WALLS_ON				///< Turn on no-slip walls (default is top, bottom, front, back unless L_WALLS_ON_2D is used)
//#define L_WALLS_ON_2D				///< Limit no-slip walls to top and bottom no-slip walls only
#define L_wall_thickness_bottom 1		///< Thickness of walls in coarsest lattice units
#define L_wall_thickness_top 1			///< Thickness of top walls in coarsest lattice units
#define L_wall_thickness_front 1		///< Thickness of front (3D) walls in coarsest lattice units
#define L_wall_thickness_back 1			///< Thickness of back (3D) walls in coarsest lattice units



/*
*******************************************************************************
********************************* Object Data *********************************
*******************************************************************************
*/

// Bounce-back solids
#define L_SOLID_BLOCK_ON			///< Add solid block to the domain

	#define L_block_on_grid_lev 0		///< Provide grid level on which block should be added 
	#define L_block_on_grid_reg 0		///< Provide grid region on which block should be added 
	// Wall labelling routine implements this
	// Specified in lattice units (i.e. by index) local to the chosen grid level
	#define L_block_x_min 30		///< Index of start of object/wall in x-direction
	#define L_block_x_max 60		///< Index of end of object/wall in x-direction
	#define L_block_y_min 30		///< Index of start of object/wall in y-direction
	#define L_block_y_max 60		///< Index of end of object/wall in y-direction
	#define L_block_z_min 30		///< Index of start of object/wall in z-direction
	#define L_block_z_max 60		///< Index of end of object/wall in z-direction


// Bounce-back objects from point clouds
//#define L_SOLID_FROM_FILE			///< Build solid body from point cloud file

	#define L_object_on_grid_lev 0		///< Provide grid level on which object should be added 
	#define L_object_on_grid_reg 0		///< Provide grid region on which object should be added
	// Following specified in lattice units (i.e. by index) local to the chosen grid level
	#define L_start_object_x 30			///< Index for start of object bounding box in X direction
	#define L_start_object_y 30			///< Index for start of object bounding box in Y direction
	#define L_centre_object_z 50		///< Index for cetnre of object bounding box in Z direction
	#define L_object_length 40			///< The object input is scaled based on this dimension
	#define L_object_scale_direction eXDirection	///< Scale in this direction (specify as enumeration)
	#define L_object_length_ref 40		///< Reference length to be used in the definition of Reynolds number


// BFL objects
//#define L_BFL_ON					///< Build BFL body from point cloud

	#define L_bfl_on_grid_lev 0		///< Provide grid level on which BFL body should be added 
	#define L_bfl_on_grid_reg 0		///< Provide grid region on which BFL body should be added
	// Following specified in lattice units (i.e. by index) local to the chosen grid level
	#define L_start_bfl_x 30		///< Index for start of object bounding box in X direction
	#define L_start_bfl_y 30		///< Index for start of object bounding box in Y direction
	#define L_centre_bfl_z 50		///< Index for cetnre of object bounding box in Z direction
	#define L_bfl_length 40			///< The BFL object input is scaled based on this dimension
	#define L_bfl_scale_direction eXDirection	///< Scale in this direction (specify as enumeration)
	#define L_bfl_length_ref 40		///< Reference length to be used in the definition of Reynolds number



/*
*******************************************************************************
****************************** Multi-grid Data ********************************
*******************************************************************************
*/

#define L_NumLev 0		///< Levels of refinement (0 = coarse grid only
#define L_NumReg 1		///< Number of refined regions (can be arbitrary if L_NumLev = 0)

#if L_NumLev != 0
// Global lattice indices (in terms of each grid level) for each refined region specified on each level


// Following options are only here to making testing different grid combinations easier
#if (L_NumReg == 2 && L_NumLev == 2) 
	const static int RefXstart[L_NumLev][L_NumReg]	= { {5, 5}, {2, 2} };
	const static int RefXend[L_NumLev][L_NumReg]	= { {25, 25}, {20, 10} };
	const static int RefYstart[L_NumLev][L_NumReg]	= { {5, 14}, {5, 2} };
	const static int RefYend[L_NumLev][L_NumReg]	= { {12, 25}, {10, 10} };
	// If doing 2D, these can be arbitrary values
	static int RefZstart[L_NumLev][L_NumReg]		= { {5, 10}, {2, 2} };
	static int RefZend[L_NumLev][L_NumReg]			= { {20, 15}, {10, 10} };

#elif (L_NumReg == 1 && L_NumLev == 1)
	const static int RefXstart[L_NumLev][L_NumReg]	= { 5 };
	const static int RefXend[L_NumLev][L_NumReg]	= { 110 };
	const static int RefYstart[L_NumLev][L_NumReg]	= { 4 };
	const static int RefYend[L_NumLev][L_NumReg]	= { 38 };
	// If doing 2D, these can be arbitrary values
	static int RefZstart[L_NumLev][L_NumReg]		= { 4 };
	static int RefZend[L_NumLev][L_NumReg]			= { 28 };

#elif (L_NumReg == 1 && L_NumLev == 2)
	const static int RefXstart[L_NumLev][L_NumReg]	= { {5}, {5} };
	const static int RefXend[L_NumLev][L_NumReg]	= { {110}, {150} };
	const static int RefYstart[L_NumLev][L_NumReg]	= { {4}, {4} };
	const static int RefYend[L_NumLev][L_NumReg]	= { {38}, {61} };
	// If doing 2D, these can be arbitrary values
	static int RefZstart[L_NumLev][L_NumReg]		= { {20}, {5} };
	static int RefZend[L_NumLev][L_NumReg]			= { {40}, {35} };

#elif (NumReg == 1 && NumLev == 3)
	const static size_t RefXstart[NumLev][NumReg]	= { {8},	{4},	{8} };
	const static size_t RefXend[NumLev][NumReg]		= { {34},	{48},	{80} };
	const static size_t RefYstart[NumLev][NumReg]	= { {9},	{4},	{8} };
	const static size_t RefYend[NumLev][NumReg]		= { {23},	{24},	{32} };
	// If doing 2D, these can be arbitrary values
	static size_t RefZstart[NumLev][NumReg]		= { {2},	{4},	{8} };
	static size_t RefZend[NumLev][NumReg]		= { {30},	{52},	{88} };

#elif (L_NumReg == 1 && L_NumLev == 3)
	const static int RefXstart[L_NumLev][L_NumReg]	= { {190},	{10},	{10} };
	const static int RefXend[L_NumLev][L_NumReg]	= { {270},	{90},	{90} };
	const static int RefYstart[L_NumLev][L_NumReg]	= { {240},	{10},	{10} };
	const static int RefYend[L_NumLev][L_NumReg]	= { {270},	{50},	{70} };
	// If doing 2D, these can be arbitrary values
	static int RefZstart[L_NumLev][L_NumReg]		= { {2},	{4},	{8} };
	static int RefZend[L_NumLev][L_NumReg]			= { {30},	{52},	{88} };

#elif (L_NumReg == 1 && L_NumLev == 4)
	const static int RefXstart[L_NumLev][L_NumReg]	= { {10},	{5},	{10},	{20} };
	const static int RefXend[L_NumLev][L_NumReg]	= { {50},	{70},	{110},	{160} };
	const static int RefYstart[L_NumLev][L_NumReg]	= { {0},	{0},	{0},	{0} };
	const static int RefYend[L_NumLev][L_NumReg]	= { {15},	{25},	{40},	{60} };
	// If doing 2D, these can be arbitrary values
	static int RefZstart[L_NumLev][L_NumReg]		= { {15},	{5},	{10},	{20} };
	static int RefZend[L_NumLev][L_NumReg]			= { {45},	{55},	{90},	{140} };

#endif

#endif


/*
*******************************************************************************
************************* Clean-up: NOT FOR EDITING ***************************
*******************************************************************************
*/

// Set default value for level and region for IB body if no subgrids
#if (defined L_IBM_ON && L_NumLev == 0)
	#undef L_IB_Lev
	#undef L_IB_Reg
	#define L_IB_Lev 0				// Grid level for immersed boundary object (0 if no refined regions)
	#define L_IB_Reg 0				// Grid region for immersed boundary object (0 if no refined regions)
#endif

// Set default value for level and region for IB body if no subgrids
#ifndef L_IBM_ON
	#undef L_IB_Lev
	#undef L_IB_Reg
	#define L_IB_Lev -1				// Grid level for immersed boundary object (-1 if no IBM)
	#define L_IB_Reg -1				// Grid region for immersed boundary object (-1 if no IBM)
#endif

// Set dependent options
#if L_dims == 3
	#define L_nVels 27		///< Number of lattice velocities

	#define L_MPI_dir 26	///< Number of MPI directions

#else
	#define L_nVels 9		// Use D2Q9

	// MPI config to 2D
	#define L_MPI_dir 8

	// Set Z limits for 2D
	#undef L_a_z
	#define L_a_z 0

	#undef L_b_z
	#define L_b_z 2

	#undef L_K
	#define L_K 1

	// Set object limits for 2D
	#undef L_block_z_min
	#define L_block_z_min 0

	#undef L_block_z_max
	#define L_block_z_max 0

	#undef L_ibb_d
	#define L_ibb_d 0

	// Set BFL start for 2D
	#undef L_centre_object_z
	#define L_centre_object_z 0

	// Set Object start for 2D
	#undef L_centre_bfl_z
	#define L_centre_bfl_z 0

	// Set IBB start for 2D
	#undef L_centre_ibb_z
	#define L_centre_ibb_z 0

	// Set z inlet velocity
	#undef L_u_0z
	#define L_u_0z 0

#endif

#if L_NumLev == 0
	// Set region info to default as no refinement
	const static int RefXstart[1][1]	= {0};
	const static int RefXend[1][1]		= {0};
	const static int RefYstart[1][1]	= {0};
	const static int RefYend[1][1]		= {0};
	static int RefZstart[1][1]			= {0};
	static int RefZend[1][1]			= {0};

	#undef L_NumReg
	#define L_NumReg 1
#endif

// Clean up for using profiled inlet
#ifdef L_USE_INLET_PROFILE
	#undef L_u_0x
	#define L_u_0x ux_in[j]
	#undef L_u_0y
	#define L_u_0y uy_in[j]
	#undef L_u_0z
	#define L_u_0z uz_in[j]
#endif

#endif
