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

#ifndef PLEINTERFACE_H
#define PLEINTERFACE_H

//#include <string>
//#include <vector>

// PLE files
//#include "ple_defs.h"
#include "ple_coupling.h"
#include "ple_locator.h"

class PLEInterface
{
private:

	//- Associated PLE locator
	ple_locator_t *locator; 

	//- Number of data points
	int numDataLocations_ = 0;

	//- Data to read from the InOutData
	std::vector<std::string> readData_;

	//- Data to write to the InOutData
	std::vector<std::string> writeData_;

	//- Configures PLE locator with the coordinates in mesh. 
	void configureLocator(InOutRepo<T>& mesh);

    //- preCICE solver interface
    //precice::SolverInterface & precice_;

    //- Mesh name used in the preCICE configuration
    //std::string meshName_;

    //- Mesh ID assigned by preCICE to the interface
    //int meshID_;

    //- Names of the file to read the mesh from
    //std::string coordFileName_;

    //- Vertex IDs assigned by preCICE
    //int * vertexIDs_;

	//- Read data IDs assigned by preCICE
	//std::vector<int> readDataIDs_;

	//- Write data IDs assigned by preCICE
	//std::vector<int> writeDataIDs_;

public:

    //- Constructor
	PLEInterface 
	(
		precice::SolverInterface &precice,
		std::string meshName,
		std::vector<std::string>& readData,
		std::vector<std::string>& writeData,
		InOutRepo<T>& mesh
	);

    //  data from the buffer 
    void readCouplingData(InOutRepo<T> &data);

    //  data and write them into the buffer
    void writeCouplingData(InOutRepo<T> &data);

    //- Destructor
    ~PLEInterface();

};


#endif
