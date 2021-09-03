/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
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

#ifndef PLEADAPTER_H
#define PLEADAPTER_H

//#ifdef _MSC_VER
//#define _NOEXCEPT noexcept
//#else
//#define NOEXCEPT noexcept
//#endif

#include "GridManager.h"
#include "MpiManager.h"
//#include "PLEInterface.h"
#include "../inc/InOutData.h"

// Other header files
// YAML reader - Used to read the adapter's configuration file.
//#include "yaml-cpp/yaml.h"

// PLE files
//#include "ple_defs.h"
#include "ple_coupling.h"
#include "ple_locator.h"

class PLEAdapter
{

private:

    //- Structure of the configuration of each coupling interface.
    //  Every interface needs to know the coupling mesh, and the kinds
    //  of data that are exchanged.
    struct InterfaceConfig
    {
        std::string meshName;
		std::string coordFileName;
        std::vector<std::string> writeData;
        std::vector<std::string> readData;
		std::vector<int> dimensions;
		std::vector<double> position;
		std::vector<double> offset;
    };

	//- Tolerance in percentage for the location of points
	double tolerance_;

    //- Configuration interfaces
    std::vector<struct InterfaceConfig> interfacesConfig_;

	//- Data to exchange with GASCANS
	InOutData* exchangeData_;

    // Configuration parameters used in the Adapter

        //- Name of the configuration file
	    std::string adapterConfigFileName_; 
	
	    //- Remember if there were errors in the read() method
       // bool errorsInConfigure = false;

        //- PLE participant name
        std::string participantName_;

        //- preCICE configuration file name
        //std::string preciceConfigFilename_;

		//- Number of coupled apps
		int numApps_;

		//- PLE app ID for LUMA
		int appID_;

		//- PLE app ID for CS
		int CSAppID_;

    //- Interfaces
    //std::vector<PLEInterface *> interfaces_;

	//- PLE locators
	std::vector<ple_locator_t *> locators_;

	//- PLE locators functions
	static ple_lnum_t meshExtents(const void *mesh, ple_lnum_t n_max_extents, double tolerance, double extents[]);
	static void pointInMesh(const void *mesh,
		float tolerance_base,
		float tolerance_fraction,
		ple_lnum_t n_points,
		const ple_coord_t point_coords[],
		const int point_tag[],
		ple_lnum_t location[],
		float distance[]);


	//- Data. Each position in the vector corresponds to a locator in the locators vector. 
	std::vector<std::vector<double>> coordinates_;  // Coordinates of the data. Format (x0,y0,z0,x1,y1,z1...,xn,yn,zn)
													// They are in the local LUMA coordinate system. LUMA_coord = World_coord - offset_
	std::vector<std::map<std::string, std::vector<double>>> vectorData_;
	std::vector<std::map<std::string, std::vector<double>>> scalarData_;
	std::vector<std::vector<double>> offset_; // Offset between the start of LUMA mesh and the start of the world coordinate system. 
	                                          // The world coordinate system is the same for LUMA and Code Saturne. 

	void exchangeMessage(const std::string &messageToSend, std::string *messageReceived);


    //- PLE data
    ple_coupling_mpi_set_t * pleSets_ = NULL;
	
	//- Sync flag
	int pleCouplingFlag_ = 0;  // Syncrhonised



	// I think that to be able to use any of this I'll have to make this class a friend of GridObj or GridManager
	GridManager * LUMAGrid_;

	// PLEAdapter also needs access to LUMA's MPI configuration
	MpiManager* lumaMpi_;

	// MPI intercommunicator with Code Saturne
	MPI_Comm mpiCS_;

	// LUMA root MPI rank
	//int lumaRootRank_;

	// Number of MPI ranks for LUMA
	//int lumaNumRanks_;

	// CS root MPI rank
	int CSRootRank_;

	// Number of MPI ranks for CS
	int CSNumRanks_;

	// LUMA-CS coupling tag
	const int  csLumaCouplingTag_ = 'C' + 'S' + '_' + 'C' + 'O' + 'U' + 'P' + 'L' + 'I' + 'N' + 'G';


    //- Solver interface initialized
    bool PLEInitialized_ = false;

    // Timesteps

        //- Coupling time step
        double couplingTimeStep_;

        //- Timestep used by the solver
        //T timestepSolver_;

    // Configuration

        //- Read the adapter's configuration file
        bool configFileRead();

		//- Create and configure PLE locator
		void addPLELocator(int i);

        //- Check the adapter's configuration file
       // bool configFileCheck();

    // Methods communicating with preCICE

        //- Initialize preCICE and exchange the first data
        //void initialize();

        //- Advance preCICE
        //void advance();

        //- Read the coupling data at each interface
        //void readCouplingData();

        //- Write the coupling data at each interface
        //void writeCouplingData();

        //- Adjust the timestep of the solver according to preCICE
        // void adjustSolverTimeStep();

        //- Determine if the coupling timestep has been completed
        //bool isCouplingTimeWindowComplete();

    //- Destroy the PLE interface and delete the allocated
    //  memory in a proper way. Called by the destructor.
    void teardown();


	//- Data handling
	void setCoordinates(int Nx, int Ny, int Nz, double x, double y, double z, double dx, int i);
	void initialiseVectorData(std::string name, double initValue, int i);
	void initialiseScalarData(std::string name, double initValue, int i);

public:

        //- Constructor
		PLEAdapter();

	    //- Initialise
	    void init(std::string adapterConfigFileName, double timestepSolver, GridManager* lumaGrid, MpiManager* lumaMpi);

        //- Setup the adapter's configuration
		bool configure();

		//- Set synchronisation flag
		void setSyncFlag(int flag);

        //- Reads data from the LUMA grid and sends it to Code_Saturne
        void sendData();

		//- Reads data received from PLE (Code_Saturne) and incorporates it to the LUMA grid. 
		void receiveData();

		//- Finalize and destroy preCICE
		void finalize();

		//- Syncronise the coupled solvers
		bool synchronise(int flags);

        //- Destructor
        ~PLEAdapter();

};


#endif
