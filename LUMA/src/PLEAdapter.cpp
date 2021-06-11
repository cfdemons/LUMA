
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

#include <iostream>
#include "../inc/PLEAdapter.h"

bool PLEAdapter::PLESynchronise()
{
	//Synchronise the coupled solvers by calling ple_coupling_mpi_set_synchronise

	return false;
}

void PLEAdapter::init(std::string adapterConfigFileName, double timestepSolver, GridManager* lumaGrid)
{
	adapterConfigFileName_ = adapterConfigFileName;
	LUMAGrid_ = lumaGrid;

    return;
}

bool PLEAdapter::configFileRead()
{
	//Read the configuration file //

	// Hardcoded for the moment
	participantName_ = "../LUMA";

	/*std::cout << "Reading the adapter's YAML configuration file " << adapterConfigFileName_ << "..." << std::endl;

    if (!configFileCheck()) return false;

    // Load the YAML file
    YAML::Node adapterConfig_ = YAML::LoadFile(adapterConfigFileName_);

    // Read the preCICE participant name
    participantName_ = adapterConfig_["participant"].as<std::string>();
    std::cout << "  participant : " << participantName_ << std::endl;

    // Read the preCICE configuration file name
    preciceConfigFilename_ = adapterConfig_["precice-config-file"].as<std::string>();
	std::cout << "  precice-config-file : " << preciceConfigFilename_ << std::endl;

    YAML::Node adapterConfigInterfaces = adapterConfig_["interfaces"];
    std::cout << "  interfaces : " << std::endl;
    for (int i = 0; i < adapterConfigInterfaces.size(); i++)
    {
        struct InterfaceConfig interfaceConfig;
        interfaceConfig.meshName = adapterConfigInterfaces[i]["mesh"].as<std::string>();
        std::cout << "  - mesh      : " << interfaceConfig.meshName << std::endl;

		interfaceConfig.coordFileName = adapterConfigInterfaces[i]["coordFileName"].as<std::string>();
		std::cout << "  - coordinate file Name      : " << interfaceConfig.coordFileName << std::endl;

        if (adapterConfigInterfaces[i]["write-data"])
        {
            std::cout << "    write-data : " << std::endl;
            if (adapterConfigInterfaces[i]["write-data"].size() > 0)
            {
                for (int j = 0; j < adapterConfigInterfaces[i]["write-data"].size(); j++)
                {
                    interfaceConfig.writeData.push_back(adapterConfigInterfaces[i]["write-data"][j].as<std::string>());
                    std::cout << "      " << interfaceConfig.writeData[j] << std::endl;
                }
            }
            else
            {
                interfaceConfig.writeData.push_back(adapterConfigInterfaces[i]["write-data"].as<std::string>());
                std::cout << "      " << interfaceConfig.writeData[i] << std::endl;
            }
        }

        if (adapterConfigInterfaces[i]["read-data"])
        {
            std::cout << "    read-data : " << std::endl;
            if (adapterConfigInterfaces[i]["read-data"].size() > 0)
            {
                for (int j = 0; j < adapterConfigInterfaces[i]["read-data"].size(); j++)
                {
                    interfaceConfig.readData.push_back(adapterConfigInterfaces[i]["read-data"][j].as<std::string>());
                    std::cout << "      " << interfaceConfig.readData[j] << std::endl;
                }
            }
            else
            {
                interfaceConfig.readData.push_back(adapterConfigInterfaces[i]["read-data"].as<std::string>());
                std::cout << "      " << interfaceConfig.readData[i] << std::endl;
            }
        }

        interfacesConfig_.push_back(interfaceConfig);
    }

    // Set the subcyclingAllowed_ switch
    if (adapterConfig_["subcycling"])
    {
        subcyclingAllowed_ = adapterConfig_["subcycling"].as<bool>();
    }
    std::cout << "    subcycling : " << subcyclingAllowed_ << std::endl;*/

    return true;
}

bool PLEAdapter::configure()
{
    //Create a PLE locator for each PLE interface in the configuration file. 
	// ple_locator_set_mesh. Do I need PLEInterface? Or just a vector of ple_locators? 
	// I'll start with just a vector of PLE locators

	// Read PLE adapter configuration file. 
	configFileRead(); 

	// Configure MPI and PLE
	int app_num = ple_coupling_mpi_name_to_id(MPI_COMM_WORLD, participantName_.c_str());
	std::cout << app_num << std::endl;

	



  return true;
}


void PLEAdapter::execute()
{
    // The solver has already solved the equations for this timestep.
    // Now call the adapter's methods to perform the coupling.

    // Perform the exchange of information
	// ple_locator_get_n_dist_points()
	// ple_locator_get_dist_locations() 
	// .... (see power point and recorded video). 


}

// I don't really know if I need this. Because I don't think PLE has a finalise function
void PLEAdapter::finalize()
{
    if (NULL != pleSets_ )
    {
	
        // Finalize the preCICE solver interface
       // precice_->finalize();

        PLEInitialized_ = false;

        // Delete the solver interface and all the related data
        teardown();
    }
    else
    {
		std::cout << "Could not finalize PLE." << std::endl;
    }

    return;
}

void PLEAdapter::teardown()
{
   
	ple_coupling_mpi_set_t * pleSets_ = NULL;

	//- PLE locator arrays (one locator for each coupled mesh with PLE). 
	std::vector<ple_locator_t*> locators_;
	std::cout << "Finalizing the PLE interface..." << std::endl;


    if (NULL != pleSets_)
    {
        std::cout << "Destroying the PLE sets..." << std::endl;
        delete pleSets_;
        pleSets_ = NULL;
    }

    // Delete the PLE locators
    if (locators_.size() > 0)
    {
        std::cout << "Deleting PLE locators..." << std::endl;
        for (int i = 0; i < locators_.size(); i++)
        {
			locators_.at(i) = ple_locator_destroy(locators_.at(i));
            delete locators_.at(i);
        }
        locators_.clear();
    }
    
    return;
}

PLEAdapter::~PLEAdapter()
{
    teardown();

    return;
}

