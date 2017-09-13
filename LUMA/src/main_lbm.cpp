/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) The University of Manchester 2017
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * further distribution commericially or otherwise without written consent.
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
#include "../inc/GridManager.h"		// Grid manager class definition
#include "../inc/ObjectManager.h"	// Object manager class definition
#include "../inc/PCpts.h"			// Point cloud class

using namespace std;	// Use the standard namespace

// Static variable declarations
std::string GridUtils::path_str;

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

#endif

	// Reset the refined region z-limits if only 2D -- must be done before initialising the MPI manager
#if (L_DIMS != 3 && L_NUM_LEVELS)
	for (int i = 0; i < L_NUM_REGIONS; i++) {
		for (int l = 0; l < L_NUM_LEVELS; l++) {
			cRefStartZ[l][i] = 0.0;
			cRefEndZ[l][i] = 0.0;
		}
	}
#endif



	/*
	****************************************************************************
	**************************** GENERAL INITIALISE ****************************
	****************************************************************************
	*/

    // Timing variables
	clock_t t_start, secs;	// Wall clock variables
	double outer_loop_time = 0.0; 

	// Start clock to time initialisation
	t_start = clock();

	// Get the time and convert it to a serial stamp for the output directory creation
	time_t curr_time = time(NULL);	// Current system date/time
	struct tm* timeinfo = localtime(&curr_time);
	char timeout_char[80];
	std::strftime(timeout_char, 80, "./output_%Y-%m-%d_%H-%M-%S", timeinfo);
	std::string path_str(timeout_char);
	GridUtils::path_str = path_str;   // Set static path variable for output directory

	// Create version string
	std::string version_string("Running LUMA -- Version ");
	version_string += LUMA_VERSION;
	int rank = GridUtils::safeGetRank();

#ifdef L_BUILD_FOR_MPI
	// Create MpiManager object (and output directory)
	MpiManager* mpim = MpiManager::getInstance();

	// Add parallel build strap
	version_string += "\n(Parallel Build: " + std::to_string(mpim->num_ranks) + " Processes)";

	if (mpim->my_rank == 0)

#else

	// Add serial strap
	version_string += "\n(Serial Build)";

	// Create directory in serial mode
	GridUtils::createOutputDirectory(GridUtils::path_str);

#endif
	{
		// Print starting string
		std::cout << version_string << std::endl;

	}

	// Create application log file
	std::ofstream logfile;
	logfile.open(GridUtils::path_str + "/log_rank" + std::to_string(rank) + ".log", std::ios::out);
	GridUtils::logfile = &logfile;	// Pass logfile reference to GridUtils class

	// TODO: Handle case when logfile doesn't open correctly
	if (!logfile.is_open()) {
		std::cout << "Logfile didn't open" << std::endl;
	}

	// Fix output format to screen
	cout.precision(L_OUTPUT_PRECISION);

	// Output start time to application log
	*GridUtils::logfile << "LUMA -- Version " << LUMA_VERSION << std::endl;
	char* time_str = ctime(&curr_time);	// Format start time as string
    *GridUtils::logfile << "Simulation started at " << time_str;	// Write start time to log

	// Create the Grid Manager
	GridManager *gm = GridManager::getInstance();

#ifdef L_BUILD_FOR_MPI
	// Decompose the domain
	mpim->mpi_gridbuild(gm);
	
	// Get time of MPI initialisation
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

	// Create the first object in the hierarchy (level = 0)
	GridObj *const Grids = new GridObj(0);


	// Log file information
	L_INFO("L0 Grid size = " + std::to_string(L_N) + "x" + std::to_string(L_M) + "x" + std::to_string(L_K), GridUtils::logfile);
#ifdef L_BUILD_FOR_MPI
	L_INFO("MPI size = " + std::to_string(L_MPI_XCORES) + "x" + std::to_string(L_MPI_YCORES) + "x" + std::to_string(L_MPI_ZCORES), GridUtils::logfile);
	*GridUtils::logfile << "Coordinates on rank " << mpim->my_rank << " are (";
	for (size_t d = 0; d < L_DIMS; d++)
	{
		*GridUtils::logfile << "\t" << mpim->rank_coords[d];
	}
	*GridUtils::logfile << "\t)" << std::endl;
#endif
	L_INFO("Number of time steps to run = " + std::to_string(L_TOTAL_TIMESTEPS), GridUtils::logfile);
	L_INFO("Grid spacing = " + std::to_string(Grids->dh), GridUtils::logfile);
	L_INFO("Time step = " + std::to_string(Grids->dt), GridUtils::logfile);
	L_INFO("Lattice viscosity = " + std::to_string(Grids->nu), GridUtils::logfile);
	L_INFO("L0 relaxation time = " + std::to_string(1.0 / Grids->omega), GridUtils::logfile);
	L_INFO("Lattice reference velocity " + std::to_string(Grids->uref), GridUtils::logfile);
	L_INFO("Lattice Mach number " + std::to_string(Grids->uref / cs), GridUtils::logfile);
	// Reynolds Number
#ifdef L_NU
#if L_NU != 0
	L_INFO("Reynolds Number = " + std::to_string(1.0 / L_NU), GridUtils::logfile);
#endif
#else
	L_INFO("Reynolds Number = " + std::to_string(L_RE), GridUtils::logfile);
#endif


	/*
	****************************************************************************
	************************** REFINEMENT INITIALISE ***************************
	****************************************************************************
	*/

	if (L_NUM_LEVELS != 0) {

		*GridUtils::logfile << "Initialising sub-grids..." << endl;

#ifdef L_INIT_VELOCITY_FROM_FILE
		/*Loading the initial velocity field from a file is icompatible with subgrids 
		  because the subgrid initializer routine calls GridObject::initVelocity(), which will read
		  the velocity data from the file and try to put it in the subgrid (without interpolating). The number of data
		  points in the file might be different from the number of cells in the subgrid, so the initialization routine
		  will crash */
		L_ERROR("Loading the initial velocity field from a file is icompatible with subgrids. Exiting.", GridUtils::logfile);
#endif

		// Loop over number of regions and add subgrids to Grids
		for (int reg = 0; reg < L_NUM_REGIONS; reg++) {

			// Try adding subgrids and let constructor initialise
			Grids->LBM_addSubGrid(reg);

		}

	}

	// Set the rankGrids value which contains a list of how many subgrids each rank has access to
	Grids->LBM_calculateRankGrids();

	// Set the pointer to the hierarchy in the Grid Manager now all grids are built
	gm->setGridHierarchy(Grids);

	/*
	****************************************************************************
	************************ OBJECT MANAGER INITIALISE *************************
	****************************************************************************
	*/

	// Create Object Manager
	ObjectManager* objMan = ObjectManager::getInstance(Grids);
	*GridUtils::logfile << "Object Manager Created." << endl;

	// Read in the geometry config file
#ifdef L_GEOMETRY_FILE
	*GridUtils::logfile << "Reading geometry configuration file..." << endl;
	objMan->io_readInGeomConfig();
#endif

#if !defined L_RESTARTING

	/* Initialise the bodies (compute support etc.) using initial body positions
	 * and compute support from supplied grid. Only attempts to initialise IBM bodies
	 * in this way. */
	objMan->ibm_initialise();

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
	Grids->io_restart(eRead);

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

#ifdef L_BUILD_FOR_MPI
	
	// Compute buffer sizes
	mpim->mpi_buffer_size();
	
	//  Build writable data for all grids and sub-grid communicators
	mpim->mpi_buildCommunicators(gm);
	
	// Compute load balance information
	mpim->mpi_updateLoadInfo(gm);

#endif

	// Write out t = 0
#ifdef L_TEXTOUT
	*GridUtils::logfile << "Writing out to <Grids->out>..." << endl;
	Grids->io_textout("INITIALISATION");	// Do not change this tag!
#endif

#ifdef L_IO_LITE
	*GridUtils::logfile << "Writing out to IOLite file..." << endl;
	Grids->io_lite(Grids->t, "INITIALISATION");
#endif

#ifdef L_HDF5_OUTPUT
	*GridUtils::logfile << "Writing out to HDF5 file..." << endl;
	Grids->io_hdf5(Grids->t);
#endif

#ifdef L_VTK_BODY_WRITE
	*GridUtils::logfile << "Writing out Bodies to VTK file..." << endl;
	objMan->io_vtkBodyWriter(Grids->t);
#endif

<<<<<<< HEAD
#ifdef L_PROBE_OUTPUT

#ifdef L_BUILD_FOR_MPI
	for (int n = 0; n < mpim->num_ranks; n++)
#endif
	{
		// Wait for rank accessing the file and only access if this rank's turn
#ifdef L_BUILD_FOR_MPI
		MPI_Barrier(mpim->world_comm);
		if (mpim->my_rank == n)
#endif
		{
			L_INFO("Initial probe write out...", GridUtils::logfile);
			Grids->io_probeOutput();
		}
	}
=======
#ifdef L_VTK_FEM_WRITE
	*GridUtils::logfile << "Writing out FEM to VTK file..." << endl;
	objMan->io_vtkFEMWriter(Grids->t);
#endif

#ifdef L_WRITE_TIP_POSITIONS
	objMan->io_writeTipPositions(Grids->t);
>>>>>>> origin/MPI_IBM
#endif

#ifdef L_BUILD_FOR_MPI
	// Barrier before recording completion of initialisation
	MPI_Barrier(mpim->world_comm);
#endif
	
	*GridUtils::logfile << "Initialising LBM time-stepping..." << std::endl;

	if (rank == 0)
		std::cout << "Initialisation complete. Starting LBM time-stepping..." << std::endl;

	
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

#ifdef L_SHOW_TIME_TO_COMPLETE
		// Start clock for timing outer loop
		t_start = clock();
#endif
		if ((Grids->t + 1) % L_OUT_EVERY == 0 && rank == 0)
			std::cout << "\rTime Step " << Grids->t + 1 << " of " << L_TOTAL_TIMESTEPS << " ------>\n" << std::flush;


		///////////////////////
		// Launch LBM Kernel //
		///////////////////////

		Grids->LBM_multi_opt();		// Launch LBM kernel on top-level grid


		///////////////
		// Write Out //
		///////////////

		// Write out here
		if (Grids->t % L_OUT_EVERY == 0)
		{
#ifdef L_BUILD_FOR_MPI
			MPI_Barrier(mpim->world_comm);
#endif
			// Write out the time an outer loop is taking to the log file
			L_INFO("Outer loop taking " + std::to_string(outer_loop_time) + 
				"ms. Approximate MLUPS for active sites only = " + 
				std::to_string(gm->activeCellOps / (outer_loop_time * 1000)),
				GridUtils::logfile);

#ifdef L_TEXTOUT
			*GridUtils::logfile << "Writing out to <Grids->out>..." << endl;
			Grids->io_textout("START OF TIMESTEP");
#endif
#ifdef L_IO_FGA
			*GridUtils::logfile << "Writing out to <.fga>..." << endl;
			Grids->io_fgaout();
#endif

#ifdef L_IO_LITE
			*GridUtils::logfile << "Writing out to IOLite file..." << endl;
			Grids->io_lite(Grids->t,"");
#endif

#ifdef L_HDF5_OUTPUT
			*GridUtils::logfile << "Writing out to HDF5 file..." << endl;
			Grids->io_hdf5(Grids->t);
#endif

#ifdef L_VTK_BODY_WRITE
			*GridUtils::logfile << "Writing out Bodies to VTK file..." << endl;
			objMan->io_vtkBodyWriter(Grids->t);
#endif

#ifdef L_VTK_FEM_WRITE
			*GridUtils::logfile << "Writing out FEM to VTK file..." << endl;
			objMan->io_vtkFEMWriter(Grids->t);
#endif

#ifdef L_WRITE_TIP_POSITIONS
			objMan->io_writeTipPositions(Grids->t);
#endif

#if (defined L_IBM_ON && defined L_IBBODY_TRACER)
			*GridUtils::logfile << "Writing out flexible body position..." << endl;
			objMan->io_writeBodyPosition(Grids->t);
#endif

		}

		// Completion time
#ifdef L_SHOW_TIME_TO_COMPLETE
		if (rank == 0 && (Grids->t % L_OUT_EVERY == 0 || Grids->t < 10))
		{
			int hms[3];
			GridUnits::secs2hms((L_TOTAL_TIMESTEPS - Grids->t) * outer_loop_time / 1000, &hms[0]);
			if (Grids->t % L_OUT_EVERY != 0) std::cout << "\r";
			std::cout << " Time to complete approx. " << hms[0] << " [h] " << hms[1] << " [m] " << hms[2] << " [s]     " << std::flush;
		}
#endif

		// Write out forces of objects
#if (defined L_LD_OUT && defined L_GEOMETRY_FILE)
		if (Grids->t % L_OUT_EVERY_FORCES == 0) {

			*GridUtils::logfile << "Writing out object lift and drag" << endl;
			objMan->io_writeForcesOnObjects(Grids->t);

#ifdef L_IBM_ON
			*GridUtils::logfile << "Writing out flexible body lift and drag..." << endl;
			objMan->io_writeLiftDrag();
#endif
		}
#endif	

		// Probe output has different frequency
#ifdef L_PROBE_OUTPUT
		if (Grids->t % L_PROBE_OUT_FREQ == 0) {

#ifdef L_BUILD_FOR_MPI
			for (int n = 0; n < mpim->num_ranks; n++)
#endif
			{

				// Wait for rank accessing the file and only access if this rank's turn
#ifdef L_BUILD_FOR_MPI
				MPI_Barrier(mpim->world_comm);
				if (mpim->my_rank == n) 
#endif
				{

					*GridUtils::logfile << "Probe write out..." << endl;
					Grids->io_probeOutput();

				}

			}

		}
#endif


		/////////////////////////
		// Restart File Output //
		/////////////////////////
		if (Grids->t % L_RESTART_OUT_FREQ == 0) {

			// Write out
			Grids->io_restart(eWrite);
		}


#ifdef L_SHOW_TIME_TO_COMPLETE
		// Update outer loop time (inc. effects of writing out for accuracy)
		outer_loop_time *= Grids->t - 1;
		outer_loop_time += (static_cast<double>(clock() - t_start)) / CLOCKS_PER_SEC * 1000;
		outer_loop_time /= Grids->t;
#endif

	// Loop End
	} while (Grids->t < L_TOTAL_TIMESTEPS);


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
#ifdef L_BUILD_FOR_MPI
	int num_ranks = mpim->num_ranks;
#else
	int num_ranks = 1;
#endif

	for (int n = 0; n < num_ranks; n++)
	{
		
#ifdef L_BUILD_FOR_MPI
		MPI_Barrier(mpim->world_comm);
#endif
		if (rank == n)
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
					GridUtils::getGrid(Grids,lev,reg,g);

					// If grid does not exist the put in a zero
					if (g == NULL) {

						timings << "\t" << 0 << "\t" << 0;

					} else {

						// Add time step time on this grid then mpi overhead time
						timings << "\t" << g->timeav_timestep << "\t" << 
#ifdef L_BUILD_FOR_MPI
							g->timeav_mpi_overhead;
#else
							0;
#endif	
							

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

	// Destroy singletons
	ObjectManager::destroyInstance();
	MpiManager::destroyInstance();
	GridManager::destroyInstance();

	// Destroy hierarchy
	delete Grids;

#ifdef L_BUILD_FOR_MPI
	// Finalise MPI
	MPI_Finalize();
#endif

	return 0;
}
