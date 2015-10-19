// LatBo.cpp : Defines the entry point for the console application.

/*
	**************************************************************************************************************
	**************************************************************************************************************
	**																											**
	**												LatBo MAIN													**
	**																											**
	**************************************************************************************************************
	**************************************************************************************************************
*/

#include "../inc/stdafx.h"
#include "../inc/definitions.h"	// Definitions file
#include "../inc/globalvars.h"		// Global variable references
#include "../inc/GridObj.h"		// Grid class
#include "../inc/MPI_manager.h"	// MPI manager class

using namespace std;	// Use the standard namespace

// Entry point
int main( int argc, char* argv[] )
{
	/*
	***************************************************************************************************************
	*********************************************** MPI INITIALISE ************************************************
	***************************************************************************************************************
	*/
#ifdef BUILD_FOR_MPI
	
	// Filename for verbose write out
#ifdef MPI_VERBOSE
	string filename;
#endif
	
	// Usual initialise
	MPI_Init( &argc, &argv );
	
#endif

	// Reset the refined region z-limits if only 2D -- must be done before initialising the MPI manager
#if (dims != 3)
	for (int i = 0; i < NumReg; i++) {
		RefZstart[i] = 0;
		RefZend[i] = 0;
	}
#endif

#ifdef BUILD_FOR_MPI
	// Create MPI_manager object
	MPI_manager mpim;

	// Initialise the topology
	mpim.mpi_init();

	// Decompose the domain
	mpim.mpi_gridbuild();
#endif

	/*
	***************************************************************************************************************
	********************************************* GENERAL INITIALISE **********************************************
	***************************************************************************************************************
	*/
	
	// Create a log file
	std::ofstream logfile;
	string fNameRank;
#ifdef BUILD_FOR_MPI
	fNameRank = to_string(mpim.my_rank);
#else
	fNameRank = to_string(0);
#endif
	logfile.open("./output/log_rank" + fNameRank + ".out", std::ios::out);

	// Fix output format
	cout.precision(6);

	// Timing variables
	clock_t t_start, t_end, secs; // Wall clock variables
	double timeav_mpi_overhead = 0.0, timeav_timestep = 0.0;	// Variables for measuring performance

	// Output start time
	time_t curr_time = time(NULL);	// Current system date/time
	char* time_str = ctime(&curr_time);	// Format as string
    logfile << "Simulation started at " << time_str;	// Write start time to log

#ifdef BUILD_FOR_MPI
	// Check that when using MPI at least 2 cores have been specified as have assumed so in implementation
	if (	Xcores < 2 || Ycores < 2
#if (dims == 3)
		|| Zcores < 2
#endif
		) {
			std::cout << "Error: See Log File." << std::endl;
			logfile << "When using MPI must use at least 2 cores in each direction. Exiting." << std::endl;
			MPI_Finalize();
			exit(EXIT_FAILURE);
	}
#endif



	/* ***************************************************************************************************************
	*********************************************** LEVEL 0 INITIALISE ***********************************************
	*************************************************************************************************************** */

	// Create the grid object (level = 0)
#ifdef BUILD_FOR_MPI
	// Call MPI constructor
	GridObj Grids(0, mpim.my_rank, mpim.local_size, mpim.global_edge_ind, mpim.global_edge_pos, mpim.MPI_coords, &logfile);
#else
	// Call basic wrapper constructor
	GridObj Grids(0, &logfile);
#endif	
	
	// Log file information
	logfile << "Grid size = " << N << "x" << M << "x" << K << endl;
#ifdef BUILD_FOR_MPI
	logfile << "MPI size = " << Xcores << "x" << Ycores << "x" << Zcores << endl;
	logfile << "Coordinates on rank " << Grids.my_rank << " are (";
		for (size_t d = 0; d < dims; d++) {
			logfile << "\t" << mpim.MPI_coords[d];
		}
		logfile << "\t)" << std::endl;
#endif
	logfile << "Number of time steps = " << T << endl;
	logfile << "Physical grid spacing = " << Grids.dt << endl;
	logfile << "Lattice viscosity = " << Grids.nu << endl;
	logfile << "L0 relaxation time = " << (1/Grids.omega) << endl;
	logfile << "Lattice inlet velocity " << Grids.gUtils.vecnorm(u_0x,u_0y,u_0z) << std::endl;
	// Reynolds Number
	logfile << "Reynolds Number = " << Re << endl;


	/* ***************************************************************************************************************
	**************************************** REFINED LEVELS INITIALISE ***********************************************
	*************************************************************************************************************** */

	if (NumLev != 0) {

		// Loop over number of regions and add subgrids to Grids
		for (int reg = 0; reg < NumReg; reg++) {

			// Try adding subgrids and let constructor initialise
			Grids.LBM_addSubGrid(reg);

		}

	}


	/* ***************************************************************************************************************
	********************************************* IBM INITIALISE *****************************************************
	*************************************************************************************************************** */

#ifdef IBM_ON
	
	logfile << "Initialising IBM..." << endl;

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
	Grids.ibm_build_body(1);
	logfile << "Case: Rectangle/Cuboid using IBM" << std::endl;

#elif defined INSERT_CIRCLE_SPHERE
	Grids.ibm_build_body(2);
	logfile << "Case: Circle/Sphere using IBM" << std::endl;

#elif defined INSERT_BOTH
	Grids.ibm_build_body(3);
	logfile << "Case: Rectangle/Cuboid + Circle/Sphere using IBM" << std::endl;

#elif defined INSERT_FILAMENT
	Grids.ibm_build_body(4);
	logfile << "Case: Single 2D filament using Jacowire IBM" << std::endl;

#elif defined INSERT_FILARRAY
	Grids.ibm_build_body(5);
	logfile << "Case: Array of filaments using Jacowire IBM" << std::endl;

#elif defined _2D_RIGID_PLATE_IBM
	Grids.ibm_build_body(6);
	logfile << "Case: 2D rigid plate using IBM" << std::endl;

#elif defined _2D_PLATE_WITH_FLAP
	Grids.ibm_build_body(7);
	logfile << "Case: 2D rigid plate using IBM with flexible flap" << std::endl;

#elif defined _3D_RIGID_PLATE_IBM
	Grids.ibm_build_body(8);
	logfile << "Case: 3D rigid plate using IBM" << std::endl;

#elif defined _3D_PLATE_WITH_FLAP
	Grids.ibm_build_body(9);
	logfile << "Case: 3D rigid plate using IBM with flexible 2D flap" << std::endl;

#endif

#if !defined RESTARTING

	// Initialise the bodies (compute support etc.) using initial body positions
	Grids.ibm_initialise();
	logfile << "Number of markers requested = " << num_markers << std::endl;

#endif

#endif


	/* ***************************************************************************************************************
	****************************************** READ IN RESTART DATA **************************************************
	*************************************************************************************************************** */

#ifdef RESTARTING

	////////////////////////
	// Restart File Input //
	////////////////////////

	// Loop over ranks and read in information one at a time
	unsigned int max_ranks;
#ifdef BUILD_FOR_MPI
	max_ranks = mpim.num_ranks;
#else
	max_ranks = 1;
#endif

	for (unsigned int n = 0; n < max_ranks; n++) {

		// Wait for rank accessing the file and only access if this rank's turn
#ifdef BUILD_FOR_MPI
		MPI_Barrier(mpim.my_comm);

		if (mpim.my_rank == n)
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
	logfile << "Reinitialising IB_bodies from restart data." << std::endl;

#endif

#endif


	/* ***************************************************************************************************************
	******************************************* CLOSE INITIALISATION *************************************************
	*************************************************************************************************************** */
	
	// Write out t = 0
#ifdef TEXTOUT
	logfile << "Writing out to <Grids.out>..." << endl;
	Grids.io_textout("INITIALISATION");	// Do not change this tag!
#endif
#ifdef VTK_WRITER
	logfile << "Writing out to VTK file..." << endl;
	Grids.vtk_writer(0.0);
#endif
	
	logfile << "Initialisation Complete." << endl << "Initialising LBM time-stepping..." << endl;


	/* ***************************************************************************************************************
	********************************************** LBM PROCEDURE *****************************************************
	*************************************************************************************************************** */
	do {

		cout << "\n------ Time Step " << Grids.t+1 << " of " << T << " ------" << endl;
		
		// Start the clock
		t_start = clock();

		///////////////////////
		// Launch LBM Kernel //
		///////////////////////
#ifdef IBM_ON
		Grids.LBM_multi(true);	// IBM requires predictor-corrector calls
#else
		Grids.LBM_multi(false);	// Just called once as no IBM
#endif

		// Print Time of loop
		t_end = clock();
		secs = t_end - t_start;
		printf("Last time step took %f second(s)\n", ((double)secs)/CLOCKS_PER_SEC);

		// Update average timestep time
		timeav_timestep *= (Grids.t-1);
		timeav_timestep += ((double)secs)/CLOCKS_PER_SEC;
		timeav_timestep /= Grids.t;

		
		///////////////
		// Write Out //
		///////////////
		if (Grids.t % out_every == 0) {
#ifdef BUILD_FOR_MPI
			MPI_Barrier(mpim.my_comm);
#endif
#ifdef TEXTOUT
			logfile << "Writing out to <Grids.out>" << endl;
			Grids.io_textout("START OF TIMESTEP");
#endif
#ifdef VTK_WRITER
			logfile << "Writing out to VTK file" << endl;
			Grids.vtk_writer(Grids.t);
#endif
#if (defined INSERT_FILAMENT || defined INSERT_FILARRAY || defined _2D_RIGID_PLATE_IBM || \
	defined _2D_PLATE_WITH_FLAP || defined _3D_RIGID_PLATE_IBM || defined _3D_PLATE_WITH_FLAP) \
	&& defined IBM_ON && defined IBBODY_TRACER
			logfile << "Writing out flexible body position" << endl;
			Grids.io_write_body_pos();
#endif
#if defined LD_OUT && defined IBM_ON
			logfile << "Writing out flexible body lift and drag" << endl;
			Grids.io_write_lift_drag();
#endif
			// Performance data
			logfile << "Time stepping taking an average of " << timeav_timestep*1000 << "ms" << std::endl;

		}


		/////////////////////////
		// Restart File Output //
		/////////////////////////
		if (Grids.t % restart_out_every == 0) {

			// Loop over ranks and write out information sequentially
			int max_ranks;
#ifdef BUILD_FOR_MPI
			max_ranks = mpim.num_ranks;
#else
			max_ranks = 1;
#endif

			for (int n = 0; n < max_ranks; n++) {

				// Wait for turn and access one rank at a time
#ifdef BUILD_FOR_MPI
				MPI_Barrier(mpim.my_comm);

				if (mpim.my_rank == n)  
#endif
				{
										
				// Write out
				Grids.io_restart(true);

				}

			}

		}


		///////////////////////
		// MPI Communication //
		///////////////////////
#ifdef BUILD_FOR_MPI
		
		// Start the clock
		t_start = clock();

		// Loop over directions in Cartesian topology
		for (int dir = 0; dir < MPI_dir; dir++) {

			////////////////////////
			// Buffer Information //
			////////////////////////

			MPI_Barrier(mpim.my_comm);
			// Pass direction and Grids by reference
			mpim.mpi_buffer_pack( dir, Grids );

			//////////////////////
			// Send Information //
			//////////////////////

			// Find opposite direction (neighbour it receives from)
			unsigned int opp_dir;
			// Opposite of even directions is +1, odd is -1 based on MPI_cartlab
			if ( (dir + 2) % 2 == 0) {
				opp_dir = dir + 1;
			} else {
				opp_dir = dir - 1;
			}

			MPI_Barrier(mpim.my_comm);

			// Send info to neighbour while receiving from other neighbour
			MPI_Sendrecv_replace( &mpim.f_buffer.front(), mpim.f_buffer.size(), MPI_DOUBLE, mpim.neighbour_rank[dir], dir, 
				mpim.neighbour_rank[opp_dir], dir, mpim.my_comm, &mpim.stat);

			

#ifdef MPI_VERBOSE
			// Write out buffer
			MPI_Barrier(mpim.my_comm);
			std::ofstream logout( "./output/mpiLog_Rank_" + std::to_string(mpim.my_rank) + ".out", std::ios::out | std::ios::app );
			logout << "Direction " << dir << "; Sending to " << mpim.neighbour_rank[dir] << "; Receiving from " << mpim.neighbour_rank[opp_dir] << std::endl;
			filename = "./output/mpiBuffer_Rank" + std::to_string(mpim.my_rank) + "_Dir" + std::to_string(dir) + ".out";
			mpim.writeout_buf(filename);
			logout.close();
#endif

			//////////////////////////////
			// Copy from buffer to grid //
			//////////////////////////////

			MPI_Barrier(mpim.my_comm);
			// Pass direction and Grids by reference
			mpim.mpi_buffer_unpack( dir, Grids );


		}

		// Print Time of MPI comms
		t_end = clock();
		secs = t_end - t_start;
		printf("MPI overhead took %f second(s)\n", ((double)secs)/CLOCKS_PER_SEC);

		// Update average MPI overhead time
		timeav_mpi_overhead *= (Grids.t-1);
		timeav_mpi_overhead += ((double)secs)/CLOCKS_PER_SEC;
		timeav_mpi_overhead /= Grids.t;

#ifdef TEXTOUT
	if (Grids.t % out_every == 0) {
		MPI_Barrier(mpim.my_comm);
		logfile << "Writing out to <Grids.out>" << endl;
		Grids.io_textout("POST MPI COMMS");
	}
#endif

	if (Grids.t % out_every == 0) {
		// Performance Data
		logfile << "MPI overhead taking an average of " << timeav_mpi_overhead*1000 << "ms" << std::endl;
	}


#endif

	// Loop End
	} while (Grids.t < T);


	/* ***************************************************************************************************************
	*********************************************** POST PROCESS *****************************************************
	*************************************************************************************************************** */


	// Close log file
	curr_time = time(NULL);			// Current system date/time and string buffer
	time_str = ctime(&curr_time);	// Format as string
	logfile << "Simulation completed at " << time_str << std::endl;		// Write end time to log file
	logfile.close();

#ifdef BUILD_FOR_MPI
	// Finalise MPI
	MPI_Finalize();
#endif

	return 0;
}
