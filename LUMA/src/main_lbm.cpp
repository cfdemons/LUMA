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

// LUMA.cpp : Defines the entry point for the console application.

#include "../inc/stdafx.h"			// Precompiled header
#include "../inc/definitions.h"		// Definitions file
#include "../inc/globalvars.h"		// Global variable references
#include "../inc/GridObj.h"			// Grid class definition
#include "../inc/MpiManager.h"		// MPI manager class definition
#include "../inc/ObjectManager.h"	// Object manager class definition

using namespace std;	// Use the standard namespace

// Static variable declarations
std::string GridUtils::path_str;
int MpiManager::MPI_coords[dims];

// Entry point
int main( int argc, char* argv[] )
{

	// Memeory leak checking
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

#ifdef BUILD_FOR_MPI

	// Usual initialise
	MPI_Init( &argc, &argv );

#else

	// When not using MPI, set max ranks to 1 and current rank to 0.
	MpiManager::num_ranks = 1;
	MpiManager::my_rank = 0;

#endif

	// Reset the refined region z-limits if only 2D -- must be done before initialising the MPI manager
#if (dims != 3 && NumLev)
	for (int i = 0; i < NumReg; i++) {
		for (int l = 0; l < NumLev; l++) {
			RefZstart[l][i] = 0;
			RefZend[l][i] = 0;
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

#ifdef BUILD_FOR_MPI

	// Create MpiManager object
	MpiManager* mpim = MpiManager::getInstance();
	
	// Create stream for Mpi Logfile
	std::ofstream mpilog;
	MpiManager::logout = &mpilog;	// Assign pointer to logfile stream to MpiManager
	
	// Initialise the topology
	mpim->mpi_init();

	// Print out version number
	if (mpim->my_rank == 0) std::cout << "Running LUMA -- Version " << LUMA_VERSION << std::endl;

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
	cout.precision(output_precision);

	// Output start time
	*GridUtils::logfile << "LUMA -- Version " << LUMA_VERSION << std::endl;
	char* time_str = ctime(&curr_time);	// Format start time as string
    *GridUtils::logfile << "Simulation started at " << time_str;	// Write start time to log

#ifdef BUILD_FOR_MPI
	// Check that when using MPI at least 2 cores have been specified as have assumed so in implementation
	if (	Xcores < 2 || Ycores < 2
#if (dims == 3)
		|| Zcores < 2
#endif
		) {
			std::cout << "Error: See Log File." << std::endl;
			*GridUtils::logfile << "When using MPI must use at least 2 cores in each direction. Exiting." << std::endl;
			MPI_Finalize();
			exit(LUMA_FAILED);
	}
#endif

	// Get time of MPI initialisation
#ifdef BUILD_FOR_MPI
	MPI_Barrier(mpim->my_comm);
	secs = clock() - t_start;
	double mpi_initialise_time = ((double)secs)/CLOCKS_PER_SEC*1000;
	*GridUtils::logfile << "MPI Topolgy initialised in "<< mpi_initialise_time << "ms." << std::endl;
#endif


	// Start clock again for next bit of initialisation
#ifdef BUILD_FOR_MPI
	MPI_Barrier(mpim->my_comm);
#endif
	t_start = clock();



	/*
	****************************************************************************
	**************************** LEVEL 0 INITIALISE ****************************
	****************************************************************************
	*/

	// Create the grid object (level = 0)
#ifdef BUILD_FOR_MPI
	// Call MPI constructor
	GridObj Grids(0, mpim->local_size, mpim->global_edge_ind, mpim->global_edge_pos);
#else
	// Call basic wrapper constructor
	GridObj Grids(0);
#endif

	// Log file information
	*GridUtils::logfile << "Grid size = " << N << "x" << M << "x" << K << endl;
#ifdef BUILD_FOR_MPI
	*GridUtils::logfile << "MPI size = " << Xcores << "x" << Ycores << "x" << Zcores << endl;
	*GridUtils::logfile << "Coordinates on rank " << MpiManager::my_rank << " are (";
		for (size_t d = 0; d < dims; d++) {
			*GridUtils::logfile << "\t" << MpiManager::MPI_coords[d];
		}
		*GridUtils::logfile << "\t)" << std::endl;
#endif
	*GridUtils::logfile << "Number of time steps = " << std::to_string(T) << endl;
	*GridUtils::logfile << "Physical grid spacing = " << std::to_string(Grids.dt) << endl;
	*GridUtils::logfile << "Lattice viscosity = " << std::to_string(Grids.nu) << endl;
	*GridUtils::logfile << "L0 relaxation time = " << std::to_string(1/Grids.omega) << endl;
	*GridUtils::logfile << "Lattice reference velocity " << std::to_string(u_ref) << std::endl;
	// Reynolds Number
	*GridUtils::logfile << "Reynolds Number = " << std::to_string(Re) << endl;


	/*
	****************************************************************************
	************************** REFINEMENT INITIALISE ***************************
	****************************************************************************
	*/

	if (NumLev != 0) {

		*GridUtils::logfile << "Initialising sub-grids..." << endl;

		// Loop over number of regions and add subgrids to Grids
		for (int reg = 0; reg < NumReg; reg++) {

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
	*GridUtils::logfile << "Object Manager Created." << endl;

#ifdef IBM_ON

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

#if defined INSERT_RECTANGLE_CUBOID
	objMan->ibm_build_body(1);
	*GridUtils::logfile << "Case: Rectangle/Cuboid using IBM" << std::endl;

#elif defined INSERT_CIRCLE_SPHERE
	objMan->ibm_build_body(2);
	*GridUtils::logfile << "Case: Circle/Sphere using IBM" << std::endl;

#elif defined INSERT_BOTH
	objMan->ibm_build_body(3);
	*GridUtils::logfile << "Case: Rectangle/Cuboid + Circle/Sphere using IBM" << std::endl;

#elif defined INSERT_FILAMENT
	objMan->ibm_build_body(4);
	*GridUtils::logfile << "Case: Single 2D filament using Jacowire IBM" << std::endl;

#elif defined INSERT_FILARRAY
	objMan->ibm_build_body(5);
	*GridUtils::logfile << "Case: Array of filaments using Jacowire IBM" << std::endl;

#elif defined _2D_RIGID_PLATE_IBM
	objMan->ibm_build_body(6);
	*GridUtils::logfile << "Case: 2D rigid plate using IBM" << std::endl;

#elif defined _2D_PLATE_WITH_FLAP
	objMan->ibm_build_body(7);
	*GridUtils::logfile << "Case: 2D rigid plate using IBM with flexible flap" << std::endl;

#elif defined _3D_RIGID_PLATE_IBM
	objMan->ibm_build_body(8);
	*GridUtils::logfile << "Case: 3D rigid plate using IBM" << std::endl;

#elif defined _3D_PLATE_WITH_FLAP
	objMan->ibm_build_body(9);
	*GridUtils::logfile << "Case: 3D rigid plate using IBM with flexible 2D flap" << std::endl;

#endif

#if !defined RESTARTING

	// Initialise the bodies (compute support etc.) using initial body positions and compute support from supplied grid
	objMan->ibm_initialise(Grids);
	*GridUtils::logfile << "Number of markers requested = " << num_markers << std::endl;

#endif

#endif

#ifdef BFL_ON


	*GridUtils::logfile << "Initialising BFL Objects..." << endl;

	// Read in input file to arrays
	PCpts* _PCpts = new PCpts();
	objMan->readInPCData(_PCpts);

	// Call BFL body builder if there are points on this rank
	if (!_PCpts->x.empty())	objMan->bfl_build_body(_PCpts);

	*GridUtils::logfile << "Finished creating BFL Objects..." << endl;
	delete _PCpts;
	

#endif

#ifdef SOLID_FROM_FILE

	*GridUtils::logfile << "Initialising Solid Objects from File..." << endl;

	// Read in data from point cloud file
	PCpts* _PCpts = new PCpts();
	objMan->readInPointData(_PCpts);
	delete _PCpts;

#endif


	/*
	****************************************************************************
	************************* INITIALISE FROM RESTART **************************
	****************************************************************************
	*/

#ifdef RESTARTING

	////////////////////////
	// Restart File Input //
	////////////////////////

	for  (int n = 0; n < MpiManager::num_ranks; n++) {

		// Wait for rank accessing the file and only access if this rank's turn
#ifdef BUILD_FOR_MPI
		MPI_Barrier(mpim->my_comm);

		if (MpiManager::my_rank == n)
#endif
		{

			// Read in
			Grids.io_restart(false);

		}

	}

	// Re-initialise IB bodies based on restart positions
#ifdef IBM_ON

	// Re-initialise the bodies (compute support etc.)
	Grids.ibm_initialise();
	*GridUtils::logfile << "Reinitialising IB_bodies from restart data." << std::endl;

#endif

#endif


	/*
	****************************************************************************
	*************************** CLOSE INITIALISATION ***************************
	****************************************************************************
	*/

	// Write out t = 0
#ifdef TEXTOUT
	*GridUtils::logfile << "Writing out to <Grids.out>..." << endl;
	Grids.io_textout("INITIALISATION");	// Do not change this tag!
#endif

#ifdef VTK_WRITER
	*GridUtils::logfile << "Writing out to VTK file..." << endl;
	Grids.io_vtkwriter(0.0);
#ifdef IBM_ON
    objMan->io_vtk_IBwriter(0.0);
#endif
#endif

#ifdef TECPLOT
		for (int n = 0; n < MpiManager::num_ranks; n++) {
			// Wait for rank accessing the file and only access if this rank's turn
#ifdef BUILD_FOR_MPI
			MPI_Barrier(mpim->my_comm);

			if (MpiManager::my_rank == n)
#endif
			{
				*GridUtils::logfile << "Writing out to TecPlot file" << endl;
				Grids.io_tecplot(Grids.t);
			}
		}

#endif

#ifdef IO_LITE
		*GridUtils::logfile << "Writing out to IOLite file" << endl;
		Grids.io_lite(Grids.t);
#endif

	// Get time of grid and object initialisation
#ifdef BUILD_FOR_MPI
	MPI_Barrier(mpim->my_comm);
#endif
	secs = clock() - t_start;
	double obj_initialise_time = ((double)secs)/CLOCKS_PER_SEC*1000;
	*GridUtils::logfile << "Grid & Object Initialisation completed in "<< obj_initialise_time << "ms." << std::endl;


	// Set the pointer to the hierarchy in the MpiManager
	MpiManager::Grids = &Grids;

#ifdef BUILD_FOR_MPI
	// Compute buffer sizes
	MPI_Barrier(mpim->my_comm);
	t_start = clock();
	
	mpim->mpi_buffer_size();	// Call buffer sizing routine
	
	secs = clock() - t_start;
	*GridUtils::logfile << "Preallocating MPI buffers completed in "<< ((double)secs)/CLOCKS_PER_SEC*1000 << "ms." << std::endl;
#endif

	*GridUtils::logfile << "Initialising LBM time-stepping..." << std::endl;


	/*
	****************************************************************************
	***************************** IB-LBM PROCEDURE *****************************
	****************************************************************************
	*/
	do {

		// Synchronise MPI processes before next time step starts
#ifdef BUILD_FOR_MPI
		MPI_Barrier(mpim->my_comm);
#endif

		cout << "\n------ Time Step " << Grids.t+1 << " of " << T << " ------" << endl;


		///////////////////////
		// Launch LBM Kernel //
		///////////////////////
#ifdef IBM_ON
		Grids.LBM_multi(true);	// IBM requires predictor-corrector calls
#else
		Grids.LBM_multi(false);	// Just called once as no IBM
#endif



		///////////////
		// Write Out //
		///////////////

		// Write out here
		if (Grids.t % out_every == 0) {
#ifdef BUILD_FOR_MPI
			MPI_Barrier(mpim->my_comm);
#endif
#ifdef TEXTOUT
			*GridUtils::logfile << "Writing out to <Grids.out>" << endl;
			Grids.io_textout("START OF TIMESTEP");
#endif

#ifdef VTK_WRITER
			*GridUtils::logfile << "Writing out to VTK file" << endl;
			Grids.io_vtkwriter(Grids.t);
#ifdef IBM_ON
            objMan->io_vtk_IBwriter(Grids.t);
#endif
#endif

#ifdef TECPLOT
			for (int n = 0; n < MpiManager::num_ranks; n++) {
				// Wait for rank accessing the file and only access if this rank's turn
#ifdef BUILD_FOR_MPI
				MPI_Barrier(mpim->my_comm);

				if (MpiManager::my_rank == n)
#endif
				{
					*GridUtils::logfile << "Writing out to TecPlot file" << endl;
					Grids.io_tecplot(Grids.t);
				}
			}
#endif

#ifdef IO_LITE
		*GridUtils::logfile << "Writing out to IOLite file" << endl;
		Grids.io_lite(Grids.t);
#endif

#if (defined INSERT_FILAMENT || defined INSERT_FILARRAY || defined _2D_RIGID_PLATE_IBM || \
	defined _2D_PLATE_WITH_FLAP || defined _3D_RIGID_PLATE_IBM || defined _3D_PLATE_WITH_FLAP) \
	&& defined IBM_ON && defined IBBODY_TRACER
			*GridUtils::logfile << "Writing out flexible body position" << endl;
			objMan->io_write_body_pos(Grids.t);
#endif
#if defined LD_OUT && defined IBM_ON
			*GridUtils::logfile << "Writing out flexible body lift and drag" << endl;
			objMan->io_write_lift_drag(Grids.t);
#endif
		}
		
		// Write out lift and drag coefficients of solid object
#ifdef COMPUTE_LIFT_AND_DRAG
		if (Grids.t % out_every_forces == 0) {
			*GridUtils::logfile << "Writing out object lift and drag" << endl;
			Grids.io_writeForceOnObject(Grids.t);
		}
#endif		

		// Probe output has different frequency
#ifdef PROBE_OUTPUT
		if (Grids.t % out_every_probe == 0) {

			for (int n = 0; n < MpiManager::num_ranks; n++) {

				// Wait for rank accessing the file and only access if this rank's turn
#ifdef BUILD_FOR_MPI
				MPI_Barrier(mpim->my_comm);
				if (MpiManager::my_rank == n) 
#endif
				{

					*GridUtils::logfile << "Probe write out" << endl;
					Grids.io_probeOutput();

				}

			}

		}
#endif


		/////////////////////////
		// Restart File Output //
		/////////////////////////
		if (Grids.t % restart_out_every == 0) {

			for (int n = 0; n < MpiManager::num_ranks; n++) {

				// Wait for turn and access one rank at a time
#ifdef BUILD_FOR_MPI
				MPI_Barrier(mpim->my_comm);

				if (MpiManager::my_rank == n)
#endif
				{

				// Write out
				Grids.io_restart(true);

				}

			}

		}

	// Loop End
	} while (Grids.t < T);


	/*
	****************************************************************************
	******************************* POST PROCESS *******************************
	****************************************************************************
	*/

#ifdef LOG_TIMINGS
	// TIMINGS FILE //
	/** Format is as follows:
	 * Mpi Init Time --- Obj Init Time --- Time Step Time L0 --- MPI Time L0 --- etc.
	 **/

	std::ofstream timings;

	// Wait for rank accessing the file and only access if this rank's turn
	for (int n = 0; n < MpiManager::num_ranks; n++) {
		
#ifdef BUILD_FOR_MPI
		MPI_Barrier(mpim->my_comm);

		if (MpiManager::my_rank == n)
#endif
		{
			if (n == 0)	timings.open(GridUtils::path_str + "/timings.out",std::ios::out);
			else timings.open(GridUtils::path_str + "/timings.out",std::ios::app);
			GridObj* g = NULL;

			// Put in initialisation times
#ifdef BUILD_FOR_MPI
			timings << mpi_initialise_time;
#else
			timings << 0;
#endif
			timings << "\t" << obj_initialise_time;

			// Loop over expected grids
			for (int lev = 0; lev <= NumLev; lev++) {
				for (int reg = 0; reg < NumReg; reg++) {
					
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

#ifdef BUILD_FOR_MPI
	// Close logfile
	MpiManager::logout->close();
	// Finalise MPI
	MPI_Finalize();
	// Destroy MpiManager
	MpiManager::destroyInstance();
#endif

	return 0;
}
