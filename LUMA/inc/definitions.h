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
#define LUMA_VERSION "1.4.0-alpha"


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
//#define L_INIT_VERBOSE			///< Write out initialisation information such as refinement mappings
//#define L_MPI_VERBOSE				///< Write out the buffers used by MPI plus more setup data
//#define L_MPI_WRITE_LOAD_BALANCE	///< Write out the load balancing information based on active cell count
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
#define L_OUT_EVERY 100			///< How many timesteps before whole grid output
#define L_OUT_EVERY_FORCES 100		///< Specific output frequency of body forces
#define L_OUTPUT_PRECISION 5		///< Precision of output (for text writers)

// Types of output
//#define L_IO_LITE					///< ASCII dump on output
#define L_HDF5_OUTPUT				///< HDF5 dump on output
//#define L_LD_OUT					///< Write out lift and drag (all bodies)
//#define L_IO_FGA                  ///< Write the components of the macroscopic velocity in a .fga file. (To be used in Unreal Engine 4).
//#define L_COMPUTE_TIME_AVERAGED_QUANTITIES

// High frequency output options
//#define L_PROBE_OUTPUT						///< Turn on probe output
#define L_PROBE_OUT_FREQ 250					///< Write out frequency of probe output
const static int cNumProbes[3] = {3, 3, 3};		///< Number of probes in each direction (x, y, z)
const static double cProbeLimsX[2] = {0.1, 0.2};	///< Limits of X plane for array of probes
const static double cProbeLimsY[2] = {0.1, 0.2};	///< Limits of Y plane for array of probes
const static double cProbeLimsZ[2] = {0.1, 0.2};	///< Limits of Z plane for array of probes


// Gravity
//#define L_GRAVITY_ON						///< Turn on gravity force
/// Expression for the gravity force
#define L_GRAVITY_FORCE 0.0001
#define L_GRAVITY_DIRECTION eXDirection		///< Gravity direction (specify using enumeration)

// Initialisation
#define L_NO_FLOW							///< Initialise the domain with no flow
//#define L_INIT_VELOCITY_FROM_FILE			///< Read initial velocity from file
//#define L_RESTARTING						///< Initialise the GridObj with quantities read from a restart file
#define L_RESTART_OUT_FREQ 50000			///< Frequency of write out of restart file

// LBM configuration
//#define L_USE_KBC_COLLISION					///< Use KBC collision operator instead of LBGK by default
//#define L_USE_BGKSMAG
#define L_CSMAG 0.07


/*
*******************************************************************************
******************************** Time data ************************************
*******************************************************************************
*/

#define L_TOTAL_TIMESTEPS 1000		///< Number of time steps to run simulation for


/*
*******************************************************************************
**************************** Domain Dimensions ********************************
*******************************************************************************
*/

// MPI Data
#define L_MPI_XCORES 2		///< Number of MPI ranks to divide domain into in X direction
#define L_MPI_YCORES 2		///< Number of MPI ranks to divide domain into in Y direction
/// Number of MPI ranks to divide domain into in Z direction.
#define L_MPI_ZCORES 2

/*
*******************************************************************************
****************************** Physical Data **********************************
*******************************************************************************
*/

// Lattice properties
#define L_DIMS 2			///< Number of dimensions to the problem
#define L_RESOLUTION 100		///< Number of coarse lattice sites per unit length
#define L_TIMESTEP 0.001		///< The timestep in non-dimensional units

// Non-dimensional domain dimensions
#define L_BX 2.2		///< End of domain in X (non-dimensional units)
#define L_BY 0.43		///< End of domain in Y (non-dimensional units)
#define L_BZ 1.0		///< End of domain in Z (non-dimensional units)

// Physical velocity
#define L_PHYSICAL_U 0.2		///< Reference velocity of the real fluid to model [m/s]


/*
*******************************************************************************
******************************** Fluid Data ***********************************
*******************************************************************************
*/

// Fluid data in lattice units
//#define L_USE_INLET_PROFILE	   ///< Use an inlet profile
//#define L_PARABOLIC_INLET	   ///< Use analytic expression for inlet profile - if not then ASCII file is read (requires L_USE_INLET_PROFILE)

// If not using an inlet profile, specify values or expressions here
#define L_UX0 0.04			///< Initial/inlet x-velocity
#define L_UY0 0.0			///< Initial/inlet y-velocity
#define L_UZ0 0.0			///< Initial/inlet z-velocity

#define L_RHOIN 1			///< Initial density. In lattice units. 
//#define L_NU 0            ///< Dimensionless kinematic viscosity L_NU = 1/Re. Comment it to use L_RE instead.  
#define L_RE 20			///< Desired Reynolds number


/*
*******************************************************************************
****************************** Immersed Boundary ******************************
*******************************************************************************
*/

// Master IBM switches //
#define L_IBM_ON						///< Turn on IBM

//#define L_STOP_EPSILON_RECOMPUTE		///< Prevent recomputing of epsilon in an attempt to save time
#define L_VTK_BODY_WRITE				///< Write out the bodies to a VTK file
#define L_GEOMETRY_FILE					///< If defined LUMA will read for geometry config file


/*
*******************************************************************************
********************************** Wall Data **********************************
*******************************************************************************
*/

// Virtual Wind Tunnels
//#define L_UPSTREAM_TUNNEL			///< Adds an inlet to all faces except exit
//#define L_FREESTREAM_TUNNEL		///< Adds a inlet to all faces


// Inlets
#define L_INLET_ON				///< Turn on inlet boundary (assumed left-hand wall - default Do Nothing)
//#define L_INLET_REGULARISED	///< Specify the inlet to be a regularised inlet condition (Latt & Chopard)
//#define L_INLET_NRBC			///< Turn on NRBC at inlet


// Outlets
#define L_OUTLET_ON				///< Turn on outlet boundary (assumed right-hand wall -- default Do Nothing)
//#define L_OUTLET_NRBC			///< Turn on NRBC at outlet


// Solids
#define L_WALLS_ON				///< Turn on no-slip walls (default is top, bottom, front, back unless L_WALLS_ON_2D is used)
#define L_WALLS_ON_2D							///< Limit no-slip walls to top and bottom no-slip walls only
#define L_WALL_THICKNESS_BOTTOM (static_cast<double>(L_BX)/static_cast<double>(L_N))		///< Thickness of wall
#define L_WALL_THICKNESS_TOP (static_cast<double>(L_BX)/static_cast<double>(L_N))			///< Thickness of top wall
#define L_WALL_THICKNESS_FRONT (static_cast<double>(L_BX)/static_cast<double>(L_N))		///< Thickness of front (3D) wall
#define L_WALL_THICKNESS_BACK (static_cast<double>(L_BX)/static_cast<double>(L_N))		///< Thickness of back (3D) wall



/*
*******************************************************************************
********************************* Object Data *********************************
*******************************************************************************
*/

// Bounce-back solids
//#define L_SOLID_BLOCK_ON			///< Add solid block to the domain

	#define L_BLOCK_ON_GRID_LEV 2		///< Provide grid level on which block should be added 
	#define L_BLOCK_ON_GRID_REG 0		///< Provide grid region on which block should be added 
	// Wall labelling routine implements this
	// Specified in lattice units (i.e. by index) local to the chosen grid level
	#define L_BLOCK_MIN_X 0.9		///< Start of object/wall in x-direction
	#define L_BLOCK_MAX_X 1.1		///< End of object/wall in x-direction
	#define L_BLOCK_MIN_Y 0.4		///< Start of object/wall in y-direction
	#define L_BLOCK_MAX_Y 0.6		///< End of object/wall in y-direction
	#define L_BLOCK_MIN_Z 0.3		///< Start of object/wall in z-direction
	#define L_BLOCK_MAX_Z 0.7		///< End of object/wall in z-direction


// Bounce-back objects from point clouds
//#define L_SOLID_FROM_FILE			///< Build solid body from point cloud file

	#define L_OBJECT_ON_GRID_LEV 2		///< Provide grid level on which object should be added 
	#define L_OBJECT_ON_GRID_REG 0		///< Provide grid region on which object should be added
	// Following specified in lattice units (i.e. by index) local to the chosen grid level
	#define L_START_OBJECT_X 0.7		///< Start of object bounding box in X direction
	#define L_START_OBJECT_Y 0.4		///< Start of object bounding box in Y direction
	#define L_CENTRE_OBJECT_Z 0.5		///< Centre of object bounding box in Z direction
	#define L_OBJECT_LENGTH 0.2			///< The object input is scaled based on this dimension
	#define L_OBJECT_SCALE_DIRECTION eXDirection	///< Scale in this direction (specify as enumeration)
	#define L_OBJECT_REF_LENGTH 0.2		///< Reference length to be used in the definition of Reynolds number


/*
*******************************************************************************
****************************** Multi-grid Data ********************************
*******************************************************************************
*/

#define L_NUM_LEVELS 1		///< Levels of refinement (0 = coarse grid only)
#define L_NUM_REGIONS 1		///< Number of refined regions (can be arbitrary if L_NUM_LEVELS = 0)

#if L_NUM_LEVELS != 0
// Position of each refined region

static double cRefStartX[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 0.1 }
};
static double cRefEndX[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 0.5 }
};
static double cRefStartY[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 0.11 }
};
static double cRefEndY[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 0.31 }
};
static double cRefStartZ[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 0.0 }
};
static double cRefEndZ[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 6.5 }
};

#endif


/*
*******************************************************************************
************************* Clean-up: NOT FOR EDITING ***************************
*******************************************************************************
*/

#define L_N static_cast<int>(L_BX * L_RESOLUTION)
#define L_M static_cast<int>(L_BY * L_RESOLUTION)
#define L_K static_cast<int>(L_BZ * L_RESOLUTION)

// Set dependent options
#if (L_DIMS == 3)

	#ifdef L_USE_KBC_COLLISION
		#define L_NUM_VELS 27		///< Number of lattice velocities
	#else
		#define L_NUM_VELS 19		///< Number of lattice velocities
	#endif

	#define L_MPI_DIRS 26	///< Number of MPI directions

#else
	#define L_NUM_VELS 9		// Use D2Q9

	// MPI config to 2D
	#define L_MPI_DIRS 8

	// Set Z limits for 2D
	#undef L_BZ
	#define L_BZ 2

	#undef L_K
	#define L_K 1

	#undef L_MPI_ZCORES
	#define L_MPI_ZCORES 1

	// Set object limits for 2D
	#undef L_BLOCK_MIN_Z
	#define L_BLOCK_MIN_Z 0.0

	#undef L_BLOCK_MAX_Z
	#define L_BLOCK_MAX_Z 0.0

	// Set BBB start for 2D
	#undef L_CENTRE_OBJECT_Z
	#define L_CENTRE_OBJECT_Z 0.0

	// Set z inlet velocity
	#undef L_UZ0
	#define L_UZ0 0.0

#endif

#if L_NUM_LEVELS == 0
	// Set region info to default as no refinement
	static double cRefStartX[1][1]	= {0.0};
	static double cRefEndX[1][1]	= {0.0};
	static double cRefStartY[1][1]	= {0.0};
	static double cRefEndY[1][1]	= {0.0};
	static double cRefStartZ[1][1]	= {0.0};
	static double cRefEndZ[1][1]	= {0.0};

	#undef L_NUM_REGIONS
	#define L_NUM_REGIONS 1
#endif

// Clean up for using profiled inlet
#ifdef L_USE_INLET_PROFILE
	#undef L_UX0
	#define L_UX0 ux_in[j]
	#undef L_UY0
	#define L_UY0 uy_in[j]
	#undef L_UZ0
	#define L_UZ0 uz_in[j]
#endif

#endif
