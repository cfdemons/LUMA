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
	const int totalloops = (int)(T/Grids.dt);		// Total number of loops to be performed (computed)


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
	******************************************* CLOSE INITIALISATION *************************************************
	*************************************************************************************************************** */

	// Write out t = 0
#ifdef ENSIGHTGOLD
	std::cout << "Writing out to EnSight file" << endl;
	Grids.genVec(fileNum);
	Grids.genScal(fileNum);
#endif
#ifdef TEXTOUT
	cout << "Writing Output to <Grids.out>..." << endl;
	Grids.lbm_write3(t);
#endif
#ifdef VTK_WRITER
	std::cout << "Writing out to VTK file" << endl;
	Grids.vtk_writer(t, 0.0);
#endif

	cout << "Initialisation Complete." << endl << "Initialising LBM time-stepping..." << endl;

	// Reynolds Number
	cout << "Reynolds Number = " << Grids.Re << endl;



	/* ***************************************************************************************************************
	********************************************** LBM PROCEDURE *****************************************************
	*************************************************************************************************************** */

	// LBM
	do {

		cout << "\n///////////////// Time Step " << t+1 << " /////////////////";

		// Start the clock
		t_start = clock();

		// Call LBM kernel from L0
		Grids.LBM_multi();

		// Print Time of loop
		t_end = clock();
		clock_t secs = t_end - t_start;
		printf("\nLast time step took %f second(s)\n", ((double)secs)/CLOCKS_PER_SEC);

		// Increment counters
		t++;
		tval += Grids.dt;


		// Write out
#ifdef ENSIGHTGOLD
		std::cout << "Writing out to EnSight file" << endl;
		fileNum++;
		Grids.genVec(fileNum);
		Grids.genScal(fileNum);
#endif
#ifdef TEXTOUT
		cout << "Writing out to <Grids.out>..." << endl;
		Grids.lbm_write3(t);
#endif
#ifdef VTK_WRITER
	std::cout << "Writing out to VTK file" << endl;
	Grids.vtk_writer(t, tval);
#endif

	} while (tval < T);


	/* ***************************************************************************************************************
	*********************************************** POST PROCESS *****************************************************
	*************************************************************************************************************** */

	// End of loop write out
#ifdef ENSIGHTGOLD
	Grids.genCase(totalloops);
	Grids.genGeo();
#endif

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
