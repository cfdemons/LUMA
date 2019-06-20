/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2018 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

/// LUMA version
#define LUMA_VERSION "1.7.7"


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
**************************** DO NOT EDIT!!!!!! ********************************
*******************************************************************************
*/
#define L_N static_cast<int>((L_BX) * L_RESOLUTION)	///< Number of coarse cells in X-direction
#define L_M static_cast<int>((L_BY) * L_RESOLUTION)	///< Number of coarse cells in Y-direction
#define L_K static_cast<int>((L_BZ) * L_RESOLUTION)	///< Number of coarse cells in Z-direction
/// Width of a coarse cell in dimensionless units
#define L_COARSE_SITE_WIDTH (1.0 / static_cast<double>(L_RESOLUTION))


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
//#define L_MOMEX_DEBUG				///< Debug momentum exchange by writing out F contributions verbosely
#define L_SHOW_TIME_TO_COMPLETE		///< Write the estimated time to completion to the terminal


/*
*******************************************************************************
************************* Global configuration data ***************************
*******************************************************************************
*/

// Using MPI?
#define L_BUILD_FOR_MPI				///< Enable MPI features in build

// Enable OMP support?
//#define L_ENABLE_OPENMP				///< Enable OpenMP features (experimental)

// Output Options
#define L_GRID_OUT_FREQ 500						///< How many timesteps before whole grid output
#define L_EXTRA_OUT_FREQ 1						///< Specific output frequency of body forces
#define L_OUTPUT_PRECISION 8					///< Precision of output (for text writers)
#define L_RESTART_OUT_FREQ (20*L_GRID_OUT_FREQ)	///< Frequency of write out of restart file
#define L_PROBE_OUT_FREQ 1000000				///< Write out frequency of probe output

// Types of output
//#define L_IO_LITE				///< ASCII dump on output
#define L_HDF5_OUTPUT				///< HDF5 dump on output
//#define L_LD_OUT				///< Write out lift and drag (all bodies)
//#define L_IO_FGA				///< Write the components of the macroscopic velocity in a .fga file. (To be used in Unreal Engine 4).
//#define L_PROBE_OUTPUT			///< Write out probe data

// Probe output options
#define L_PROBE_NUM_X 0						///< Number of probes in X direction
#define L_PROBE_NUM_Y 0						///< Number of probes in Y direction
#define L_PROBE_NUM_Z 0						///< Number of probes in Z direction
#define L_PROBE_MIN_X 0.5					///< Start position of probe array in X direction
#define L_PROBE_MIN_Y (0.4 + L_WALL_THICKNESS_BOTTOM)		///< Start position of probe array in Y direction
#define L_PROBE_MIN_Z 0.0					///< Start position of probe array in Z direction
#define L_PROBE_MAX_X 1.5					///< End position of probe array in X direction
#define L_PROBE_MAX_Y (1.6 + L_WALL_THICKNESS_BOTTOM)		///< End position of probe array in Y direction
#define L_PROBE_MAX_Z 0.0					///< End position of probe array in Z direction

// Forcing
#define L_GRAVITY_ON						///< Turn on gravity force
#define L_GRAVITY_FORCE 2.986e-4			///< Expression for the gravity force in dimensionless units
#define L_GRAVITY_DIRECTION eXDirection		///< Gravity direction (specify using enumeration)

// Initialisation
#define L_NO_FLOW					///< Initialise the domain with no flow
//#define L_INIT_VELOCITY_FROM_FILE			///< Read initial velocity from file
//#define L_RESTARTING					///< Initialise the GridObj with quantities read from a restart file

// LBM configuration
//#define L_USE_KBC_COLLISION				///< Use KBC collision operator instead of LBGK by default
//#define L_USE_BGKSMAG
#define L_CSMAG 0.3

/// Compute the time-averaged values of velocity, density and the velocity products.
//#define L_COMPUTE_TIME_AVERAGED_QUANTITIES


/*
*******************************************************************************
******************************** Time data ************************************
*******************************************************************************
*/

#define L_TOTAL_TIMESTEPS 2000				///< Number of time steps to run simulation for


/*
*******************************************************************************
**************************** Domain Dimensions ********************************
*******************************************************************************
*/

// MPI Data
#define L_MPI_XCORES 4		///< Number of MPI ranks to divide domain into in X direction
#define L_MPI_YCORES 2		///< Number of MPI ranks to divide domain into in Y direction
#define L_MPI_ZCORES 1		///< Number of MPI ranks to divide domain into in Z direction.

// Decomposition strategy
//#define L_MPI_SMART_DECOMPOSE		///< Use smart decomposition to improve load balancing
#define L_MPI_SD_MAX_ITER 1000		///< Max number of iterations to be used for smart decomposition algorithm

// Topology report
//#define L_MPI_TOPOLOGY_REPORT		///< Have the MPI Manager report on different combinations of X Y Z cores
#define L_MPI_TOP_XCORES 12			///< Max number of X MPI ranks to use for the topology report
#define L_MPI_TOP_YCORES 12			///< Max number of Y MPI ranks to use for the topology report
#define L_MPI_TOP_ZCORES 12			///< Max number of Z MPI ranks to use for the topology report

/*
*******************************************************************************
****************************** Physical Data **********************************
*******************************************************************************
*/

// Lattice properties
#define L_DIMS 3							///< Number of dimensions to the problem
#define L_RESOLUTION 10						///< Number of coarse lattice sites per unit length
#define L_TIMESTEP 0.1 / L_RESOLUTION		///< The timestep in non-dimensional units

// Non-dimensional domain dimensions
#define L_BX 5.0				///< End of domain in X (non-dimensional units)
#define L_BY 1.0				///< End of domain in Y (non-dimensional units)
#define L_BZ 1.0				///< End of domain in Z (non-dimensional units)

// Physical velocity
#define L_PHYSICAL_U 1.0		///< Reference velocity of the real fluid to model [m/s]

// Reference density	
#define L_PHYSICAL_RHO 1.225		///< Reference density in physical units


/*
*******************************************************************************
******************************** Fluid Data ***********************************
*******************************************************************************
*/

// Fluid data in lattice units
//#define L_USE_INLET_PROFILE		///< Use an inlet profile
//#define L_PARABOLIC_INLET		///< Use analytical parabolic inlet profile

// If not using an inlet profile, specify values or expressions here
#define L_UX0 1.0			///< Initial/inlet x-velocity
#define L_UY0 0.0			///< Initial/inlet y-velocity
#define L_UZ0 0.0			///< Initial/inlet z-velocity

#define L_RHOIN 1			///< Initial density. In lattice units. 
//#define L_NU 0			///< Dimensionless kinematic viscosity L_NU = 1/Re. Comment it to use L_RE instead.
#define L_RE 1000	///< Desired Reynolds number
//#define L_REYNOLDS_RAMP 1000	///< Defines over how many time steps to ramp the Reynolds number


/*
*******************************************************************************
****************************** Object Management ******************************
*******************************************************************************
*/

// General //
//#define L_GEOMETRY_FILE					///< If defined LUMA will read for geometry config file
//#define L_VTK_BODY_WRITE				///< Write out the bodies to a VTK file
//#define L_VTK_FEM_WRITE				///< Write out the FEM bodies to a VTK file

// IBM //
//#define L_IBM_ON				///< Turn on IBM
//#define L_UNIVERSAL_EPSILON_CALC		///< Do universal epsilon calculation (should be used if supports from different bodies overlap)

// FEM //
#define L_NB_ALPHA 0.25				///< Parameter for Newmark-Beta time integration (0.25 for 2nd order)
#define L_NB_DELTA 0.5				///< Parameter for Newmark-Beta time integration (0.5 for 2nd order)
#define L_RELAX 0.5				///< Under-relaxation for FSI coupling
//#define L_WRITE_TIP_POSITIONS			///< Turn on writing out filament tip positions (only works on flexible filaments)

/*
*******************************************************************************
********************************** Wall Data **********************************
*******************************************************************************
*/

// BC types (set to eFluid for periodic)
#define L_WALL_LEFT		eFluid			///< BC used on the left of the domain
#define L_WALL_RIGHT	eFluid			///< BC used on the right of the domain
#define L_WALL_BOTTOM	eSolid			///< BC used on the bottom of the domain
#define L_WALL_TOP		eSolid			///< BC used on the top of the domain
#define L_WALL_FRONT	eFluid			///< BC used on the front of the domain
#define L_WALL_BACK		eFluid			///< BC used on the bottom of the domain

// BC qualifiers
//#define L_REGULARISED_BOUNDARIES	///< Specify the velocity and pressure BCs to be regularised (Latt & Chopard)
//#define L_VELOCITY_RAMP 2			///< Defines time in dimensionless units over which to ramp up the inlet velocity
//#define L_PRESSURE_DELTA 0.0		///< Sets a desired pressure fluctuation away from L_RHOIN for a pressure boundary

// General
#define L_WALL_THICKNESS_BOTTOM (1.0 * L_COARSE_SITE_WIDTH)	///< Thickness of wall
#define L_WALL_THICKNESS_TOP (1.0 * L_COARSE_SITE_WIDTH)	///< Thickness of top wall
#define L_WALL_THICKNESS_LEFT (1.0 * L_COARSE_SITE_WIDTH)	///< Thickness of left wall
#define L_WALL_THICKNESS_RIGHT (1.0 * L_COARSE_SITE_WIDTH)	///< Thickness of right wall
#define L_WALL_THICKNESS_FRONT (1.0 * L_COARSE_SITE_WIDTH)	///< Thickness of front (3D) wall
#define L_WALL_THICKNESS_BACK (1.0 * L_COARSE_SITE_WIDTH)	///< Thickness of back (3D) wall


/*
*******************************************************************************
****************************** Multi-grid Data ********************************
*******************************************************************************
*/

#define L_NUM_LEVELS 1		///< Levels of refinement (0 = coarse grid only)
#define L_NUM_REGIONS 1		///< Number of refined regions (can be arbitrary if L_NUM_LEVELS = 0)
#define L_AUTO_SUBGRIDS		///< Activate auto sub-grid generation using the padding parameters below

// Auto-sub-grid configuration (if you want coincident edges then set to (-2.0 * dh))
#define L_PADDING_X_MIN 0.1		///< Padding between X start of each sub-grid and its child edge
#define L_PADDING_X_MAX 0.1		///< Padding between X end of each sub-grid and its child edge
#define L_PADDING_Y_MIN 0.1		///< Padding between Y start of each sub-grid and its child edge
#define L_PADDING_Y_MAX 0.1		///< Padding between Y end of each sub-grid and its child edge
#define L_PADDING_Z_MIN (-2.0 * dh)		///< Padding between Z start of each sub-grid and its child edge
#define L_PADDING_Z_MAX (2.0 * dh)		///< Padding between Z end of each sub-grid and its child edge

#if L_NUM_LEVELS != 0
// Position of each refined region

static double cRefStartX[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ .1 }
};
static double cRefEndX[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ 2.9 }
};
static double cRefStartY[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ .05 }
};
static double cRefEndY[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ .95 }
};
static double cRefStartZ[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ 0.05 }
};
static double cRefEndZ[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ 0.95 }
};

#endif


/*
*******************************************************************************
************************* Clean-up: NOT FOR EDITING ***************************
*******************************************************************************
*/

// Set probes
const static int cNumProbes[3] = { L_PROBE_NUM_X, L_PROBE_NUM_Y, L_PROBE_NUM_Z };
const static double cProbeLimsX[2] = { L_PROBE_MIN_X, L_PROBE_MAX_X };	///< Limits of X plane for array of probes
const static double cProbeLimsY[2] = { L_PROBE_MIN_Y, L_PROBE_MAX_Y };	///< Limits of Y plane for array of probes
const static double cProbeLimsZ[2] = { L_PROBE_MIN_Z, L_PROBE_MAX_Z };	///< Limits of Z plane for array of probes

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
#define L_BZ 0

#undef L_K
#define L_K 1

#undef L_MPI_ZCORES
#define L_MPI_ZCORES 1

// Set object limits for 2D
#undef L_BLOCK_MIN_Z
#define L_BLOCK_MIN_Z 0.0

#undef L_BLOCK_MAX_Z
#define L_BLOCK_MAX_Z 0.0

// Set z inlet velocity
#undef L_UZ0
#define L_UZ0 0.0

#endif

#if L_NUM_LEVELS == 0
// Set region info to default as no refinement
static double cRefStartX[1][1] = { 0.0 };
static double cRefEndX[1][1] = { 0.0 };
static double cRefStartY[1][1] = { 0.0 };
static double cRefEndY[1][1] = { 0.0 };
static double cRefStartZ[1][1] = { 0.0 };
static double cRefEndZ[1][1] = { 0.0 };

#undef L_NUM_REGIONS
#define L_NUM_REGIONS 1
#endif

#endif