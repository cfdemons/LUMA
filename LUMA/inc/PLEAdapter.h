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
#include "PLEInterface.h"
//#include "InOutData.h"

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
    };

    //- Configuration interfaces
    std::vector<struct InterfaceConfig> interfacesConfig_;

	//- Data to exchange with GASCANS
	InOutData<T>* exchangeData_;

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

		//- PLE app ID
		int appID_;

		//- PLE app ID for CS
		int CSAppID_;

    //- Interfaces
    //std::vector<PLEInterface *> interfaces_;

	//- PLE locators
	std::vector<ple_locator_t *> locators_;

    //- PLE data
    ple_coupling_mpi_set_t * pleSets_ = NULL;
	
	//- Sync flag
	int pleCouplingFlag_ = 0;

	//- PLE locator arrays (one locator for each coupled mesh with PLE). 
	//std::vector<ple_locator_t*> locators_; 

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
		void addPLELocator();

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

public:

        //- Constructor
	    void init(std::string adapterConfigFileName, double timestepSolver, GridManager* lumaGrid, MpiManager* lumaMpi);

        //- Setup the adapter's configuration
		bool configure();

        //- Called at each time step. Returns the data to exchange with LUMIS. 
        void execute();

		//- Finalize and destroy preCICE
		void finalize();

		//- Syncronise the coupled solvers
		bool synchronise(int flags);

        //- Destructor
        ~PLEAdapter();

};


#endif
