
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

bool PLEAdapter::synchronise(int flags)
{

	//Synchronise the coupled solvers by calling ple_coupling_mpi_set_synchronise

	int i;

	int sync_flags = 0;
	int stop_mask = pleCouplingFlag_ & PLE_COUPLING_STOP;
	int leader_id = -1;
	double ts_min = -1.;

	int n_apps = ple_coupling_mpi_set_n_apps(pleSets_);
	int app_id = ple_coupling_mpi_set_get_app_id(pleSets_);

	int reset_flags[] = { PLE_COUPLING_NEW_ITERATION,
						 PLE_COUPLING_REDO_ITERATION };

	const int *app_status = NULL;
	const double *app_ts = NULL;

	ple_coupling_mpi_set_info_t ai;

	/* Set synchronization flag */

	app_status = ple_coupling_mpi_set_get_status(pleSets_);

	sync_flags = app_status[app_id];
	for (i = 0; i < 2; i++) {
		if (sync_flags & reset_flags[i])
			sync_flags -= reset_flags[i];
	}
	sync_flags = sync_flags | flags | stop_mask;

	if (LUMAGrid_->Grids->t >= L_TOTAL_TIMESTEPS)
		sync_flags = sync_flags | PLE_COUPLING_STOP;
	else {
		sync_flags = sync_flags | PLE_COUPLING_NEW_ITERATION;
		if (LUMAGrid_->Grids->t == L_TOTAL_TIMESTEPS - 1)
			sync_flags = sync_flags | PLE_COUPLING_LAST;
	}

	// LUMA is not planned to have redo iteration capability. 
	/*if (flags & PLE_COUPLING_REDO_ITERATION) {
		if (sync_flags & PLE_COUPLING_NEW_ITERATION)
			sync_flags -= PLE_COUPLING_NEW_ITERATION;
		if (sync_flags & PLE_COUPLING_STOP)
			sync_flags -= PLE_COUPLING_STOP;
	}*/

	/* Synchronize applications */

	printf("LUMA: Hi before ple synchronise. \n");

	ple_coupling_mpi_set_synchronize(pleSets_, sync_flags, LUMAGrid_->Grids->dt);

	printf("CS: Hi after ple synchronise. \n");

	app_status = ple_coupling_mpi_set_get_status(pleSets_);
	
	app_ts = ple_coupling_mpi_set_get_timestep(pleSets_);

	/* Check if we should use the smallest time step */
	//if (app_status[app_id] & PLE_COUPLING_TS_MIN)
	//	ts_min = LUMAGrid_->Grids->dt;

	/* Loop on applications */

	for (i = 0; i < n_apps; i++) {

		if (app_status[i] & PLE_COUPLING_NO_SYNC)
			continue;

		/* Handle leader or minimum time step update */

		/*if (app_status[i] & PLE_COUPLING_TS_LEADER) {
			if (leader_id > -1) {
				ple_coupling_mpi_set_info_t ai_prev = ple_coupling_mpi_set_get_info(pleSets_, ts_min);
				ai = ple_coupling_mpi_set_get_info(pleSets_, i);
				
				L_INFO("Application " + ai.app_name + " ("+ ai.app_type + ") tried to set the group time step," +
					   " application" + ai_prev.app_name + " (" + ai_prev.app_type + ") has already done so.", GridUtils::logfile);
			}
			else {
				leader_id = i;
				couplingTimeStep_ = app_ts[i] / _cs_coupling_ts_multiplier;
			}
		}
		else if (app_status[i] & PLE_COUPLING_TS_MIN) {
			if (ts_min > 0)
				ts_min = CS_MIN(ts_min, app_ts[i]);
		}*/

		/* Handle time stepping behavior */

		if (app_status[i] & PLE_COUPLING_STOP)
		{
			if (L_TOTAL_TIMESTEPS > LUMAGrid_->Grids->dt)
			{
				ai = ple_coupling_mpi_set_get_info(pleSets_ , i);
				L_INFO("Application " + std::string(ai.app_name) + " (" + std::string(ai.app_type) + ") requested calculation stop.", GridUtils::logfile);
				
				//TODO: This sets the current time step to the maximum number of time steps. The simulation should stop. 
				LUMAGrid_->Grids->t = L_TOTAL_TIMESTEPS;
			}
		}
		/*else if (app_status[i] & PLE_COUPLING_REDO_ITERATION) {
			ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
			bft_error
			(__FILE__, __LINE__, 0,
				_("\nApplication \"%s\" (%s) requested restarting iteration,\n"
					"but this is not currently handled."),
				ai.app_name, ai.app_type);
		}*/
		else if (!(app_status[i] & PLE_COUPLING_NEW_ITERATION))
		{
			ai = ple_coupling_mpi_set_get_info(pleSets_, i);
			L_ERROR("Application " + std::string(ai.app_name) + " (" + std::string(ai.app_type) + ") synchronized with status flag " +
				     std::to_string(app_status[i]) + " which does not specify a known behavior.", GridUtils::logfile);
		}

		if (app_status[i] & PLE_COUPLING_LAST) {
			if (L_TOTAL_TIMESTEPS > LUMAGrid_->Grids->t + 1)
			{
				ai = ple_coupling_mpi_set_get_info(pleSets_, i);
				L_INFO("Application " + std::string(ai.app_name) + " (" + std::string(ai.app_type) + ") requested last iteration.", GridUtils::logfile);
				LUMAGrid_->Grids->t = L_TOTAL_TIMESTEPS -1;
			}
		}

	} /* end of loop on applications */

	//if (ts_min > 0)
	//	*ts = ts_min / _cs_coupling_ts_multiplier;

	return true;
}

void PLEAdapter::init(std::string adapterConfigFileName, double timestepSolver, GridManager* lumaGrid, MpiManager* lumaMpi)
{
	adapterConfigFileName_ = adapterConfigFileName;
	LUMAGrid_ = lumaGrid;
	lumaMpi_ = lumaMpi;

    return;
}

bool PLEAdapter::configFileRead()
{
	//Read the configuration file //

	// Hardcoded for the moment
	participantName_ = "LUMA";

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

void PLEAdapter::addPLELocator()
{
	//-  Create PLE locator
	locators_.push_back(ple_locator_create(mpiCS_, CSNumRanks_, CSRootRank_));

	//- Call PLE locator set mesh

}

bool PLEAdapter::configure()
{
    //Create a PLE locator for each PLE interface in the configuration file. 
	// ple_locator_set_mesh. Do I need PLEInterface? Or just a vector of ple_locators? 
	// I'll start with just a vector of PLE locators

	// Read PLE adapter configuration file. 
	configFileRead(); 

	// ** DEBUG!!!! This data should be read from the Config file but I hard code it here just to test. **//
	std::string app_type = "LUMA";
	std::string app_name = "LEFT"; // Name of the domain. It is in the mpi command line, or I can also write it in the config file?
								   // I think it is better to read from mpi command line so that there is only one config file. I don't like it but it is 
	                               // what CS does.
	InterfaceConfig CS_inlet_position;
	CS_inlet_position.dimensions.push_back(1); 
	CS_inlet_position.dimensions.push_back(20);
	CS_inlet_position.dimensions.push_back(10);

	CS_inlet_position.meshName = "CS_inlet";
	CS_inlet_position.position.push_back(1.5);
	CS_inlet_position.position.push_back(0.0);
	CS_inlet_position.position.push_back(0.0);

	CS_inlet_position.readData.push_back(" ");
	CS_inlet_position.writeData.push_back("v");

	interfacesConfig_.push_back(CS_inlet_position);


	pleSets_ = ple_coupling_mpi_set_create(pleCouplingFlag_, "LUMA", "LEFT", MPI_COMM_WORLD, lumaMpi_->ple_comm);  //or lumaMpi_->world_comm?

	numApps_ = ple_coupling_mpi_set_n_apps(pleSets_);
	appID_ = ple_coupling_mpi_set_get_app_id(pleSets_);

	int local_range[2] = { -1, -1 };
	int distant_range[2] = { -1, -1 };

	// Search for the MPI ranks of CS. 
	// Stop if you find more than one CS instance in the MPI call. 
	int CSNum = 0;
	for (int i = 0; i < numApps_; i++)
	{
		ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(pleSets_, i);
		//std::cout << " LUMA: App ID " << std::string(ai.app_type) << std::endl;
		if (std::string(ai.app_type).find("Code_Saturne") != std::string::npos)
		{
			if (CSNum < 1)
			{ 
				CSRootRank_ = ai.root_rank;
				CSNumRanks_ = ai.n_ranks;
				CSAppID_ = i;
				CSNum++;
			}
			else
			{
				L_ERROR("Multiple instances of Code_Saturne found. LUMA is only programmed to be coupled with a single Code Saturne instance.", GridUtils::logfile);
			}
			
		}
	}

	// Create the intercommunicator between LUMA and Code Saturne
	// I'm only initialising one CS 
	ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD, lumaMpi_->ple_comm, CSRootRank_, &mpiCS_, local_range, distant_range);

	std::cout << "LUMA: I've created the PLE intracomm!" << std::endl;

	// Create the interfaces
	std::cout << "LUMA: Creating PLE interfaces..." << std::endl;
	for (int i = 0; i < interfacesConfig_.size(); i++)
	{
		// Create the data storage for the interface
		// Precice doesn't need to know the offset between LUMA coordinates and world coordinates for the moment. 
		// coordFileName is in world coordinates. I'm not sure about this coordFileName with PLE coupling. 
		exchangeData_->addRepo(interfacesConfig_.at(i).coordFileName, 0, 0, 0);

		// Create space for the data to be written / read for the current repo. 
		for (int j = 0; j < interfacesConfig_.at(i).readData.size(); j++)
		{
			if ((interfacesConfig_.at(i).readData.at(j).find("v") != std::string::npos) || (interfacesConfig_.at(i).readData.at(j).find("V") != std::string::npos))
			{
				exchangeData_->getRepo(i).initialiseVectorData(interfacesConfig_.at(i).readData.at(j));
			}
			else if ((interfacesConfig_.at(i).readData.at(j).find("t") != std::string::npos) || (interfacesConfig_.at(i).readData.at(j).find("T") != std::string::npos))
			{
				exchangeData_->getRepo(i).initialiseScalarData(interfacesConfig_.at(i).readData.at(j));
			}
			else
			{
				L_ERROR("Data type for " + interfacesConfig_.at(i).readData.at(j) + " not recognised. Data in .yml configuration file has to contain v or V for velocity and t or T for temperature.", GridUtils::logfile);
			}
		}
		for (int j = 0; j < interfacesConfig_.at(i).writeData.size(); j++)
		{
			if ((interfacesConfig_.at(i).writeData.at(j).find("v") != std::string::npos) || (interfacesConfig_.at(i).writeData.at(j).find("V") != std::string::npos) )
			{
				exchangeData_->getRepo(i).initialiseVectorData(interfacesConfig_.at(i).writeData.at(j));
			}
			else if ((interfacesConfig_.at(i).writeData.at(j).find("t") != std::string::npos) || (interfacesConfig_.at(i).writeData.at(j).find("T") != std::string::npos))
			{
				exchangeData_->getRepo(i).initialiseScalarData(interfacesConfig_.at(i).writeData.at(j));
			}
			else
			{
				L_ERROR("Data type for " + interfacesConfig_.at(i).writeData.at(j) + " not recognised. Data in .yml configuration file has to contain v or V for velocity and t or T for temperature.", GridUtils::logfile);
			}
		}

		// Create and configure the PLE locator
		addPLELocator();

		std::cout << "Interface created on mesh" << interfacesConfig_.at(i).meshName << std::endl;

	}




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
   
	ple_coupling_mpi_set_destroy(&pleSets_);
	
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

