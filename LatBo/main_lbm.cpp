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

using namespace std;	// Use the standard namespace

// Entry point (Compiling in VS on Windows)
int _tmain( )
{
	/*
	***************************************************************************************************************
	********************************************* GENERAL INITIALISE **********************************************
	***************************************************************************************************************
	*/

#if (dims != 3)
	// Reset the refined region z-limits if only 2D
	for (int i = 0; i < NumReg; i++) {
		RefZstart[i] = 0;
		RefZend[i] = 0;
	}
#endif

	// Create a log file
	std::ofstream logfile;
	logfile.open("./Output/log.out", std::ios::out);

	// Fix output format
	cout.precision(4);

	// Timing variables
	clock_t t_start, t_end; // Wall clock variables
	int t = 0;				// Time step counter, initially zero
	double tval = 0;		// Actual value of physical time (for multi-grid do not necessarily have unit time step)
	int fileNum = 0;		// Output file number (1 per timestep)



	/* ***************************************************************************************************************
	*********************************************** LEVEL 0 INITIALISE ***********************************************
	*************************************************************************************************************** */

	// Create grid level 0 and let constructor initialise
	GridObj Grids(0);

	// Number of loops based on L0 time step
	const int totalloops = (int)(T/Grids.dt);		// Total number of loops to be performed (computed)#
	logfile << "Number of time steps = " << totalloops << endl;
	logfile << "Lattice viscosity = " << Grids.nu << endl;
	logfile << "L0 relaxation time = " << (1/Grids.omega) << endl;

	/* ***************************************************************************************************************
	**************************************** REFINED LEVELS INITIALISE ***********************************************
	*************************************************************************************************************** */

	if (NumLev != 0) {

		// Loop over number of regions and add subgrids to Grids
		for (int reg = 0; reg < NumReg; reg++) {

			// Add the subgrid and let constructor initialise
			Grids.LBM_addSubGrid(reg);
		}

	}


	/* ***************************************************************************************************************
	********************************************* IBM INITIALISE *****************************************************
	*************************************************************************************************************** */

#ifdef IBM_ON
	
	logfile << "Initialising IBM..." << endl;

	// Build a body -- type == 1 is a rectangle/cuboid, type == 2 is a circle/sphere
#if defined INSERT_RECTANGLE_CUBOID
	Grids.build_body(1);
#elif defined INSERT_CIRCLE_SPHERE
	Grids.build_body(2);
#endif

	// Initialise the body (compute support etc.)
	Grids.ibm_initialise();

#endif


	/* ***************************************************************************************************************
	******************************************* CLOSE INITIALISATION *************************************************
	*************************************************************************************************************** */

	// Write out t = 0
#ifdef ENSIGHTGOLD
	logfile << "Writing out to EnSight file..." << endl;
	Grids.genVec(fileNum);
	Grids.genScal(fileNum);
#endif
#ifdef TEXTOUT
	logfile << "Writing Output to <Grids.out>..." << endl;
	Grids.LBM_textout(t);
#endif
#ifdef VTK_WRITER
	logfile << "Writing out to VTK file..." << endl;
	Grids.vtk_writer(t, 0.0);
#endif

	logfile << "Initialisation Complete." << endl << "Initialising LBM time-stepping..." << endl;

	// Reynolds Number
	logfile << "Reynolds Number = " << Re << endl;


	/* ***************************************************************************************************************
	********************************************** LBM PROCEDURE *****************************************************
	*************************************************************************************************************** */

	// LBM
	do {

		cout << "\n------ Time Step " << t+1 << " of " << totalloops << " ------" << endl;

		// Start the clock
		t_start = clock();

		// Call LBM kernel from L0
#ifdef IBM_ON
		Grids.LBM_multi(true);	// IBM requires predictor-corrector calls
#else
		Grids.LBM_multi(false);	// Just called once as no IBM
#endif

		// Print Time of loop
		t_end = clock();
		clock_t secs = t_end - t_start;
		printf("Last time step took %f second(s)\n", ((double)secs)/CLOCKS_PER_SEC);

		// Increment counters
		t++;
		tval += Grids.dt;


		// Write out
		if (t % out_every == 0) {
#ifdef ENSIGHTGOLD
		logfile << "Writing out to EnSight file" << endl;
		fileNum++;
		Grids.genVec(fileNum);
		Grids.genScal(fileNum);
#endif
#ifdef TEXTOUT
		logfile << "Writing out to <Grids.out>..." << endl;
		Grids.LBM_textout(t);
#endif
#ifdef VTK_WRITER
	logfile << "Writing out to VTK file" << endl;
	Grids.vtk_writer(t, tval);
#endif
		}

	} while (tval < T);


	/* ***************************************************************************************************************
	*********************************************** POST PROCESS *****************************************************
	*************************************************************************************************************** */

	// End of loop write out
#ifdef ENSIGHTGOLD
	Grids.genCase(totalloops);
	Grids.genGeo();
#endif

	// Close log file
	logfile.close();

	return 0;
}

// Entry point (compiling in Code::Blocks on Linux)
#ifndef _WIN32
int main()
{

    _tmain();

    return 0;
}
#endif
