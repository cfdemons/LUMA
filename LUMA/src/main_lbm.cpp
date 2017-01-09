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

/// \file main_lbm.cpp
/// \mainpage
///
/// --------------------------------------------------------------
///
/// ------ Lattice Boltzmann @ The University of Manchester ------
///
/// -------------------------- L-U-M-A ---------------------------
///
///  Copyright (C) 2015, 2016
///  E-mail contact: info@luma.manchester.ac.uk
///
/// This software is for academic use only and not available for
/// distribution without written consent.

#include "../inc/stdafx.h"			// Precompiled header
#include "../inc/GridObj.h"			// Grid class definition
#include "../inc/MpiManager.h"		// MPI manager class definition
#include "../inc/ObjectManager.h"	// Object manager class definition
#include "../inc/GridUtils.h"		// Grid utilities
#include "../inc/PCpts.h"			// Point cloud class

using namespace std;	// Use the standard namespace

// Static variable declarations
std::string GridUtils::path_str;
int MpiManager::MPI_coords[L_DIMS];

/// Entry point for the application
int main( int argc, char* argv[] )
{

	// Memory leak checking
#ifdef _DEBUG
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	// Set output to the terminal window
	_CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
	_CrtSetReportFile( _CRT_WARN, _CRTDBG_FILE_STDOUT );
	_CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
	_CrtSetReportFile( _CRT_ERROR, _CRTDBG_FILE_STDOUT );
	_CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
	_CrtSetReportFile( _CRT_ASSERT, _CRTDBG_FILE_STDOUT );
#endif


	/*
	****************************************************************************
	****************************** MPI INITIALISE ******************************
	****************************************************************************
	*/

#ifdef L_BUILD_FOR_MPI

	// Usual initialise
	MPI_Init( &argc, &argv );

#else

	// When not using MPI, set max ranks to 1 and current rank to 0.
	MpiManager::num_ranks = 1;
	MpiManager::my_rank = 0;

#endif

	// Reset the refined region z-limits if only 2D -- must be done before initialising the MPI manager
#if (L_DIMS != 3 && L_NUM_LEVELS)
	for (int i = 0; i < L_NUM_REGIONS; i++) {
		for (int l = 0; l < L_NUM_LEVELS; l++) {
			cRefStartZ[l][i] = 0;
			cRefEndZ[l][i] = 0;
		}
	}
#endif



	/*
	****************************************************************************
	**************************** GENERAL INITIALISE ****************************
	****************************************************************************
	*/

    // Timing variables
	clock_t t_start, secs; // Wall clock variables

	// Start clock to time initialisation
	t_start = clock();

	// Get the time and convert it to a serial stamp for the output directory creation
	time_t curr_time = time(NULL);	// Current system date/time
	struct tm* timeinfo = localtime(&curr_time);
	char timeout_char[80];
	std::strftime(timeout_char, 80, "./output_%Y-%m-%d_%H-%M-%S", timeinfo);
	std::string path_str(timeout_char);
	GridUtils::path_str = path_str;   // Set static path variable for output directory

#ifdef L_BUILD_FOR_MPI

	// Create MpiManager object
	MpiManager* mpim = MpiManager::getInstance();
	
	// Create stream for Mpi Logfile
	std::ofstream mpilog;
	MpiManager::logout = &mpilog;	// Assign pointer to logfile stream to MpiManager
	
	// Initialise the topology
	mpim->mpi_init();

	// Print out version number
	if (mpim->my_rank == 0) {
		std::cout << "Running LUMA -- Version " << LUMA_VERSION << std::endl;
#ifdef L_BUILD_FOR_MPI
		std::cout << "(Parallel Build: " << mpim->num_ranks << " Processes)" << std::endl;
#else
		std::cout << "(Serial Build)" << std::endl;
#endif
	}

	// Decompose the domain
	mpim->mpi_gridbuild();
#else

	// Create directory
	int result = GridUtils::createOutputDirectory(path_str);

	// Print out version number
	std::cout << "Running LUMA -- Version " << LUMA_VERSION << std::endl;

#endif
	
	// Create application log file
	std::ofstream logfile;
	logfile.open(GridUtils::path_str + "/log_rank" + to_string(MpiManager::my_rank) + ".out", std::ios::out);
	GridUtils::logfile = &logfile;	// Pass logfile reference to GridUtils class

	// TODO: Handle case when logfile doesn't open correctly
	if (!logfile.is_open()) {
		std::cout << "Logfile didn't open" << std::endl;
	}

	// Fix output format to screen
	cout.precision(L_OUTPUT_PRECISION);

	// Output start time
	*GridUtils::logfile << "LUMA -- Version " << LUMA_VERSION << std::endl;
	char* time_str = ctime(&curr_time);	// Format start time as string
    *GridUtils::logfile << "Simulation started at " << time_str;	// Write start time to log

#ifdef L_BUILD_FOR_MPI
	// Check that when using MPI at least 2 cores have been specified as have assumed so in implementation
	if (	L_MPI_XCORES < 2 || L_MPI_YCORES < 2
#if (L_DIMS == 3)
		|| L_MPI_ZCORES < 2
#endif
		) {
			std::cout << "Error: See Log File." << std::endl;
			*GridUtils::logfile << "When using MPI must use at least 2 cores in each direction. Exiting." << std::endl;
			MPI_Finalize();
			exit(LUMA_FAILED);
	}
#endif

	// Get time of MPI initialisation
#ifdef L_BUILD_FOR_MPI
	MPI_Barrier(mpim->world_comm);
	secs = clock() - t_start;
	double mpi_initialise_time = ((double)secs)/CLOCKS_PER_SEC*1000;
	*GridUtils::logfile << "MPI Topolgy initialised in "<< mpi_initialise_time << "ms." << std::endl;
#endif


	// Start clock again for next bit of initialisation
#ifdef L_BUILD_FOR_MPI
	MPI_Barrier(mpim->world_comm);
#endif
	t_start = clock();



	/*
	****************************************************************************
	**************************** LEVEL 0 INITIALISE ****************************
	****************************************************************************
	*/

	// Create the grid object (level = 0)
#ifdef L_BUILD_FOR_MPI
	// Call MPI constructor
	GridObj Grids(0, mpim->local_size, mpim->global_edge_ind, mpim->global_edge_pos);
#else
	// Call basic wrapper constructor
	GridObj Grids(0);
#endif

	// Log file information
	*GridUtils::logfile << "Grid size = " << L_N << "x" << L_M << "x" << L_K << endl;
#ifdef L_BUILD_FOR_MPI
	*GridUtils::logfile << "MPI size = " << L_MPI_XCORES << "x" << L_MPI_YCORES << "x" << L_MPI_ZCORES << endl;
	*GridUtils::logfile << "Coordinates on rank " << MpiManager::my_rank << " are (";
		for (size_t d = 0; d < L_DIMS; d++) {
			*GridUtils::logfile << "\t" << MpiManager::MPI_coords[d];
		}
		*GridUtils::logfile << "\t)" << std::endl;
#endif
	*GridUtils::logfile << "Number of time steps = " << std::to_string(L_TIMESTEPS) << endl;
	*GridUtils::logfile << "Physical grid spacing = " << std::to_string(Grids.dt) << endl;
	*GridUtils::logfile << "Lattice viscosity = " << std::to_string(Grids.nu) << endl;
	*GridUtils::logfile << "L0 relaxation time = " << std::to_string(1/Grids.omega) << endl;
	*GridUtils::logfile << "Lattice reference velocity " << std::to_string(L_UREF) << std::endl;
	// Reynolds Number
	*GridUtils::logfile << "Reynolds Number = " << std::to_string(L_RE) << endl;


	/*
	****************************************************************************
	************************** REFINEMENT INITIALISE ***************************
	****************************************************************************
	*/

	if (L_NUM_LEVELS != 0) {

		*GridUtils::logfile << "Initialising sub-grids..." << endl;

		// Loop over number of regions and add subgrids to Grids
		for (int reg = 0; reg < L_NUM_REGIONS; reg++) {

			// Try adding subgrids and let constructor initialise
			Grids.LBM_addSubGrid(reg);

		}

	}

	/*
	****************************************************************************
	************************ OBJECT MANAGER INITIALISE *************************
	****************************************************************************
	*/

	// Create Object Manager
	ObjectManager* objMan = ObjectManager::getInstance(&Grids);
	PCpts* _PCpts = NULL;
	*GridUtils::logfile << "Object Manager Created." << endl;

#ifdef L_IBM_ON

	*GridUtils::logfile << "Initialising IBM Objects..." << endl;

	// Build a body
	//		body_type == 1 is a rectangle/cuboid with rigid IBM,
	//		body_type == 2 is a circle/sphere with rigid IBM,
	//		body_type == 3 is a multi-body test case featuring both the above with rigid IBM
	//		body_type == 4 is a single inextensible flexible filament with Jacowire IBM
	//		body_type == 5 is an array of flexible filaments with Jacowire IBM
	//		body_type == 6 is a plate in 2D with rigid IBM
	//		body_type == 7 is the same as the previous case but with a Jacowire flexible flap added to the trailing edge
	//		body_type == 8 is a plate in 3D with rigid IBM
	//		body_type == 9 is the same as the previous case but with a rigid but moving filament array commanded by a single 2D Jacowire filament

#if defined L_INSERT_RECTANGLE_CUBOID
	objMan->ibm_buildBody(1);
	*GridUtils::logfile << "Case: Rectangle/Cuboid using IBM" << std::endl;

#elif defined L_INSERT_CIRCLE_SPHERE
	objMan->ibm_buildBody(2);
	*GridUtils::logfile << "Case: Circle/Sphere using IBM" << std::endl;

#elif defined L_INSERT_BOTH
	objMan->ibm_buildBody(3);
	*GridUtils::logfile << "Case: Rectangle/Cuboid + Circle/Sphere using IBM" << std::endl;

#elif defined L_INSERT_FILAMENT
	objMan->ibm_buildBody(4);
	*GridUtils::logfile << "Case: Single 2D filament using Jacowire IBM" << std::endl;

#elif defined L_INSERT_FILARRAY
	objMan->ibm_buildBody(5);
	*GridUtils::logfile << "Case: Array of filaments using Jacowire IBM" << std::endl;

#elif defined L_2D_RIGID_PLATE_IBM
	objMan->ibm_buildBody(6);
	*GridUtils::logfile << "Case: 2D rigid plate using IBM" << std::endl;

#elif defined L_2D_PLATE_WITH_FLAP
	objMan->ibm_buildBody(7);
	*GridUtils::logfile << "Case: 2D rigid plate using IBM with flexible flap" << std::endl;

#elif defined L_3D_RIGID_PLATE_IBM
	objMan->ibm_buildBody(8);
	*GridUtils::logfile << "Case: 3D rigid plate using IBM" << std::endl;

#elif defined L_3D_PLATE_WITH_FLAP
	objMan->ibm_buildBody(9);
	*GridUtils::logfile << "Case: 3D rigid plate using IBM with flexible 2D flap" << std::endl;

#endif

#ifdef L_IBB_FROM_FILE

	*GridUtils::logfile << "Initialising IB Body from File..." << endl;

	// Read in data from point cloud file
	_PCpts = new PCpts();
	objMan->io_readInCloud(_PCpts, eIBBCloud);
	delete _PCpts;
	*GridUtils::logfile << "Finished creating IBB Objects..." << endl;

#endif

#if !defined L_RESTARTING

	// Initialise the bodies (compute support etc.) using initial body positions and compute support from supplied grid
	objMan->ibm_initialise();

#endif

#endif // End IBM_ON

#ifdef L_BFL_ON


	*GridUtils::logfile << "Initialising BFL Objects from File..." << endl;

	// Read in input file to arrays
	_PCpts = new PCpts();
	objMan->io_readInCloud(_PCpts, eBFLCloud);
	delete _PCpts;
	*GridUtils::logfile << "Finished creating BFL Objects..." << endl;
	

#endif

#ifdef L_SOLID_FROM_FILE

	*GridUtils::logfile << "Initialising Solid Objects from File..." << endl;

	// Read in data from point cloud file
	_PCpts = new PCpts();
	objMan->io_readInCloud(_PCpts, eBBBCloud);
	delete _PCpts;
	*GridUtils::logfile << "Finished creating Solid Objects..." << endl;

#endif


	/*
	****************************************************************************
	************************* INITIALISE FROM RESTART **************************
	****************************************************************************
	*/

#ifdef L_RESTARTING

	////////////////////////
	// Restart File Input //
	////////////////////////

	// Read in
	Grids.io_restart(eRead);

	// Reinitialise IB bodies based on restart positions
#ifdef L_IBM_ON

	// Reinitialise the bodies (compute support etc.)
	ObjectManager::getInstance()->ibm_initialise();
	*GridUtils::logfile << "Reinitialising IB_bodies from restart data." << std::endl;

#endif

#endif


	/*
	****************************************************************************
	*************************** CLOSE INITIALISATION ***************************
	****************************************************************************
	*/

	// Get time of grid and object initialisation
#ifdef L_BUILD_FOR_MPI
	MPI_Barrier(mpim->world_comm);
#endif
	secs = clock() - t_start;
	double obj_initialise_time = ((double)secs)/CLOCKS_PER_SEC*1000;
	*GridUtils::logfile << "Grid & Object Initialisation completed in "<< obj_initialise_time << "ms." << std::endl;


	// Set the pointer to the hierarchy in the MpiManager
	MpiManager::Grids = &Grids;


#ifdef L_BUILD_FOR_MPI
	
	// Compute buffer sizes
	mpim->mpi_buffer_size();

	//  Build halo descriptors and sub-grid communicators
	mpim->mpi_buildCommunicators();

	// Compute load balance information
	mpim->mpi_updateLoadInfo();

#endif


	// Write out t = 0
#ifdef L_TEXTOUT
	*GridUtils::logfile << "Writing out to <Grids.out>..." << endl;
	Grids.io_textout("INITIALISATION");	// Do not change this tag!
#endif

#ifdef L_IO_LITE
	*GridUtils::logfile << "Writing out to IOLite file..." << endl;
	Grids.io_lite(Grids.t, "INITIALISATION");
#endif

#ifdef L_HDF5_OUTPUT
	*GridUtils::logfile << "Writing out to HDF5 file..." << endl;
	Grids.io_hdf5(Grids.t);
#endif


	*GridUtils::logfile << "Initialising LBM time-stepping..." << std::endl;


	/*
	****************************************************************************
	***************************** IB-LBM PROCEDURE *****************************
	****************************************************************************
	*/
	do {

		// Synchronise MPI processes before next time step starts
#ifdef L_BUILD_FOR_MPI
		MPI_Barrier(mpim->world_comm);
#endif

		if (MpiManager::my_rank == 0 && (Grids.t+1) % L_OUT_EVERY == 0)
			std::cout << "\n------ Time Step " << Grids.t + 1 << " of " << L_TIMESTEPS << " ------" << endl;


		///////////////////////
		// Launch LBM Kernel //
		///////////////////////

#ifdef L_USE_OPTIMISED_KERNEL
		Grids.LBM_multi_opt();		// Launch LBM kernel on top-level grid
#else
		Grids.LBM_multi();			// Main LBM kernel (same both with and without IBM)
#endif

		///////////////
		// Write Out //
		///////////////

		// Write out here
		if (Grids.t % L_OUT_EVERY == 0) {
#ifdef L_BUILD_FOR_MPI
			MPI_Barrier(mpim->world_comm);
#endif
#ifdef L_TEXTOUT
			*GridUtils::logfile << "Writing out to <Grids.out>..." << endl;
			Grids.io_textout("START OF TIMESTEP");
#endif
#ifdef L_IO_FGA
			*GridUtils::logfile << "Writing out to <.fga>..." << endl;
			Grids.io_fgaout();
#endif

#ifdef L_IO_LITE
			*GridUtils::logfile << "Writing out to IOLite file..." << endl;
			Grids.io_lite(Grids.t,"");
#endif

#ifdef L_HDF5_OUTPUT
			*GridUtils::logfile << "Writing out to HDF5 file..." << endl;
			Grids.io_hdf5(Grids.t);
#endif

#if (defined L_IBM_ON && defined L_VTK_BODY_WRITE)
			objMan->io_vtkIBBWriter(Grids.t);
#endif

#if (defined L_INSERT_FILAMENT || defined L_INSERT_FILARRAY || defined L_2D_RIGID_PLATE_IBM || \
	defined L_2D_PLATE_WITH_FLAP || defined L_3D_RIGID_PLATE_IBM || defined L_3D_PLATE_WITH_FLAP) \
	&& defined L_IBM_ON && defined L_IBBODY_TRACER
			*GridUtils::logfile << "Writing out flexible body position..." << endl;
			objMan->io_writeBodyPosition(Grids.t);
#endif
		}

		// Write out forces of objects
#ifdef L_LD_OUT
		if (Grids.t % L_OUT_EVERY_FORCES == 0) {

			*GridUtils::logfile << "Writing out object lift and drag" << endl;
			objMan->io_writeForceOnObject(Grids.t);

#ifdef L_IBM_ON
			*GridUtils::logfile << "Writing out flexible body lift and drag..." << endl;
			objMan->io_writeLiftDrag(Grids.t);
#endif
		}
#endif	

		// Probe output has different frequency
#ifdef L_PROBE_OUTPUT
		if (Grids.t % L_PROBE_OUT_FREQ == 0) {

			for (int n = 0; n < MpiManager::num_ranks; n++) {

				// Wait for rank accessing the file and only access if this rank's turn
#ifdef L_BUILD_FOR_MPI
				MPI_Barrier(mpim->world_comm);
				if (MpiManager::my_rank == n) 
#endif
				{

					*GridUtils::logfile << "Probe write out..." << endl;
					Grids.io_probeOutput();

				}

			}

		}
#endif


		/////////////////////////
		// Restart File Output //
		/////////////////////////
		if (Grids.t % L_RESTART_OUT_FREQ == 0) {

			// Write out
			Grids.io_restart(eWrite);
		}

	// Loop End
	} while (Grids.t < L_TIMESTEPS);


	/*
	****************************************************************************
	******************************* POST PROCESS *******************************
	****************************************************************************
	*/

#ifdef L_LOG_TIMINGS
	// TIMINGS FILE //
	/* Format is as follows:
	 * Mpi Init Time --- Obj Init Time --- Time Step Time L0 --- MPI Time L0 --- etc.
	 */

	std::ofstream timings;

	// Wait for rank accessing the file and only access if this rank's turn
	for (int n = 0; n < MpiManager::num_ranks; n++) {
		
#ifdef L_BUILD_FOR_MPI
		MPI_Barrier(mpim->world_comm);

		if (MpiManager::my_rank == n)
#endif
		{
			if (n == 0)	timings.open(GridUtils::path_str + "/timings.out",std::ios::out);
			else timings.open(GridUtils::path_str + "/timings.out",std::ios::app);
			GridObj* g = NULL;

			// Put in initialisation times
#ifdef L_BUILD_FOR_MPI
			timings << mpi_initialise_time;
#else
			timings << 0;
#endif
			timings << "\t" << obj_initialise_time;

			// Loop over expected grids
			for (int lev = 0; lev <= L_NUM_LEVELS; lev++) {
				for (int reg = 0; reg < L_NUM_REGIONS; reg++) {
					
					// Get the grid
					g = NULL;
					GridUtils::getGrid(MpiManager::Grids,lev,reg,g);

					// If grid does not exist the put in a zero
					if (g == NULL) {

						timings << "\t" << 0 << "\t" << 0;

					} else {

						// Add time step time on this grid then mpi overhead time
						timings << "\t" << g->timeav_timestep << "\t" << g->timeav_mpi_overhead;

					}

				}
			}

		}

		timings << std::endl;
		timings.close();
	}

	// END TIMINGS FILE //
#endif


	// Close log file
	curr_time = time(NULL);			// Current system date/time and string buffer
	time_str = ctime(&curr_time);	// Format as string
	*GridUtils::logfile << "Simulation completed at " << time_str << std::endl;		// Write end time to log file
	logfile.close();

	// Destroy ObjectManager
	ObjectManager::destroyInstance();

#ifdef L_BUILD_FOR_MPI
	// Close logfile
	MpiManager::logout->close();
	// Finalise MPI
	MPI_Finalize();
#endif

	return 0;
}
