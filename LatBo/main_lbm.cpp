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

#include "stdafx.h"
#include "definitions.h"	// Definitions file
#include "globalvars.h"		// Global variable references
#include "GridObj.h"		// Grid class
#include "MPI_manager.h"	// MPI manager class

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
	logfile.open("./Output/log_rank" + fNameRank + ".out", std::ios::out);

	// Fix output format
	cout.precision(4);

	// Timing variables
	clock_t t_start, t_end, secs; // Wall clock variables
	double tval = 0;		// Actual value of physical time (for multi-grid do not necessarily have unit time step)

	// Output start time
	time_t curr_time = time(NULL);	// Current system date/time
	char* time_str = ctime(&curr_time);	// Format as string
    logfile << "Simulation started at " << time_str;	// Write start time to log



	/* ***************************************************************************************************************
	*********************************************** LEVEL 0 INITIALISE ***********************************************
	*************************************************************************************************************** */

	// Create the grid object (level = 0)
#ifdef BUILD_FOR_MPI
	// Call MPI constructor
	GridObj Grids(0, mpim.my_rank, mpim.local_size, mpim.global_edge_ind, mpim.global_edge_pos, &logfile);
#else
	// Call basic wrapper constructor
	GridObj Grids(0, &logfile);
#endif

	// Number of loops based on L0 time step
	const int totalloops = (int)(T/Grids.dt);		// Total number of loops to be performed (computed)#
	
	
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
	logfile << "Number of time steps = " << totalloops << endl;
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

			/* Check to make sure the refined region lies completely on this rank's L0 grid is performed
			 * in the init_refined_lab() routine. If the code has got this far without exiting then
			 * we are OK to add the subgrids.
			 */

			// Add the subgrid and let constructor initialise
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

	// Initialise the bodies (compute support etc.)
	Grids.ibm_initialise();

	logfile << "Number of markers requested = " << num_markers << std::endl;

#endif


	/* ***************************************************************************************************************
	******************************************* CLOSE INITIALISATION *************************************************
	*************************************************************************************************************** */
	
	// Write out t = 0
#ifdef TEXTOUT
	logfile << "Writing Output to <Grids.out>..." << endl;
	Grids.io_textout();
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

		cout << "\n------ Time Step " << Grids.t+1 << " of " << totalloops << " ------" << endl;
		
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

		// Increment physical time
		tval += Grids.dt;

		
		///////////////
		// Write Out //
		///////////////
		if (Grids.t % out_every == 0) {
#ifdef BUILD_FOR_MPI
			MPI_Barrier(mpim.my_comm);
#endif
#ifdef TEXTOUT
			logfile << "Writing out to <Grids.out>" << endl;
			Grids.io_textout();
#endif
#ifdef VTK_WRITER
			logfile << "Writing out to VTK file" << endl;
			Grids.vtk_writer(tval);
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

		}

#ifdef BUILD_FOR_MPI
		///////////////////////
		// MPI Communication //
		///////////////////////
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
			MPI_Sendrecv_replace( &mpim.f_buffer.front(), mpim.f_buffer.size(), MPI_DOUBLE, mpim.neighbour_rank[dir], 999, 
				mpim.neighbour_rank[opp_dir], 999, mpim.my_comm, &mpim.stat);

#ifdef MPI_VERBOSE
			// Write out buffer
			MPI_Barrier(mpim.my_comm);
			filename = "./Output/Buffer_Rank" + std::to_string(mpim.my_rank) + "_Dir" + std::to_string(dir) + ".out";
			mpim.writeout_buf(filename);
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

#ifdef TEXTOUT
		MPI_Barrier(mpim.my_comm);
		logfile << "Writing out (post-MPI) to <Grids.out>" << endl;
		Grids.io_textout();
#endif

#endif

	// Loop End
	} while (tval < T);


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
