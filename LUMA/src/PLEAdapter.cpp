
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

	printf("LUMA: Hi after ple synchronise. \n");

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

PLEAdapter::PLEAdapter()
{
}

void PLEAdapter::init(std::string adapterConfigFileName, double timestepSolver, GridManager* lumaGrid, MpiManager* lumaMpi)
{
	adapterConfigFileName_ = adapterConfigFileName;
	LUMAGrid_ = lumaGrid;
	lumaMpi_ = lumaMpi;

    return;
}

ple_lnum_t PLEAdapter::meshExtents(const void * mesh, ple_lnum_t n_max_extents, double tolerance, double extents[])
{

	if (mesh == NULL)
		return 0;

	// In query mode, return maximum extents available  (currently limited to 1)
	if (n_max_extents < 0)
		return 1;

	// Calculate the mesh extents for the current rank
	MpiManager *mpim = MpiManager::getInstance();
	int rank = GridUtils::safeGetRank();

	extents[0] = mpim->rank_core_edge[eXMin][rank];
	extents[1] = mpim->rank_core_edge[eYMin][rank];

#if (L_DIMS == 3)
	extents[2] = mpim->rank_core_edge[eZMin][rank];
	extents[3] = mpim->rank_core_edge[eXMax][rank];
	extents[4] = mpim->rank_core_edge[eYMax][rank];
	extents[5] = mpim->rank_core_edge[eZMax][rank];
#else
	extents[2] = mpim->rank_core_edge[eXMax][rank];
	extents[3] = mpim->rank_core_edge[eYMax][rank];
#endif

	// Add tolerance to the mesh extents. As done by fvm_nodal_extract.c , _elt_extents_finalize in Code_Saturne
	double delta[3];
	for (int i = 0; i < L_DIMS; i++)
		delta[i] = (extents[i + L_DIMS] - extents[i]) * tolerance;

	for (int i = 0; i < L_DIMS; i++)
	{
		extents[i] = extents[i] - delta[i];
		extents[i + L_DIMS] = extents[i + L_DIMS] + delta[i];
	}
	
	// DEBUGGG!!
/* 	for (int i = 0; i < L_DIMS*2; i++)
	{
		std::cout << "LUMA: extents " + std::to_string(i) + " = " + std::to_string(extents[i]) << std::endl;
	} */

	return 1;
}

/// This function is used to mapping from defined values to LUMA mesh
bool PLEAdapter::pointInMesh(const void * mesh, float tolerance_base, float tolerance_fraction, ple_lnum_t n_points, const ple_coord_t point_coords[], const int point_tag[], ple_lnum_t location[], float distance[], int inter_id)
{
	const GridObj* grid = static_cast<const GridObj*>(mesh);
	bool incells = false;
	L_INFO("Number of points to locate " + std::to_string(n_points), GridUtils::logfile);
	
	for (int p = 0; p < n_points; p++)
	{
		// Set the point coordinates in the LUMA coordinate system
		double px, py, pz;

		px = point_coords[3 * p] + L_PLE_OFFSET_X;
		py = point_coords[3 * p + 1] + L_PLE_OFFSET_Y;
		pz = point_coords[3 * p + 2 ] + L_PLE_OFFSET_Z;

		eLocationOnRank loc;
		std::vector<int> pos = { -1, -1, -1 };
		// Check that the point is not in the sending or receiving layers. (I think we don't want the point to be in a sending or receiving halo). 
		if (GridUtils::isOnThisRank(px, py, pz, &loc, grid, &pos))
		{
			// If the point is in the receiving halo, it is counted as unlocated. 
			if (loc == eHalo)
			{
				L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates " + std::to_string(px) +  " " + std::to_string(py) + " " + std::to_string(pz) + " is in the halo of rank " + std::to_string(GridUtils::safeGetRank()), GridUtils::logfile);
				location[p] = -1;
				distance[p] = -1;
			}
			else if (loc == eCore) // Point is in the core of the rank's meshes, including the sending halo. 
			{
				L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates " + std::to_string(px) + " " + std::to_string(py) + " " + std::to_string(pz) + " is inside the mesh of rank " + std::to_string(GridUtils::safeGetRank()), GridUtils::logfile);

				// Flatten the index returned by inOnThisRank
				int id = pos[2] + pos[1] * grid->K_lim + pos[0] * grid->K_lim * grid->M_lim;
				location[p] = id;

				// TODO: To fill the distance array I have to get the LUMA coordinates of id and then calculate the absolute distance between them
				// and the coordinates in point_coords. I think it is just XPos[pos[0]], YPos[pos[1]], ZPos[pos[2]]. Then divide it by dx and this should be smaller than 1
				distance[p] = GridUtils::vecnorm(px - grid->XPos.at(pos[0]),
					py - grid->YPos.at(pos[1]),
					pz - grid->ZPos.at(pos[2]));// / grid->dh;
				
				if (distance[p] > 1)
				{
					L_ERROR("The distance between PLE point with coordinates x = " + std::to_string(px) +
						" y = " + std::to_string(py) +
						" z = " + std::to_string(pz) +
						" and the LUMA mesh position x = " + std::to_string(grid->XPos.at(pos[0])) +
						" y = " + std::to_string(grid->YPos.at(pos[1])) +
						" z = " + std::to_string(grid->ZPos.at(pos[2])) +
						" is " + std::to_string(distance[p]) + "which is bigger than 1. But the point should be inside the cell",
						GridUtils::logfile);
				}
				incells = true;
			}
			else
			{
				L_ERROR("Trying to locate a point from PLE mesh. The point is in the mesh but its location is eNone", GridUtils::logfile);
			}
		}
		else // If the point is not in the current rank, set the distance to -1 and the index of the point
		{
			location[p] = -1;
			distance[p] = -1;

			L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates "+ std::to_string(px) + " " + std::to_string(py) + " " + std::to_string(pz) + " is outside the mesh of rank " + std::to_string(GridUtils::safeGetRank()), GridUtils::logfile);
		}

		L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates " + std::to_string(px) + " " + std::to_string(py) + " " + std::to_string(pz) + " distance " + std::to_string(distance[p]), GridUtils::logfile);
	}
	return incells;
}

/// This function is used to establish mapping from CS to LUMA
void PLEAdapter::pointInMeshMapping(const void * mesh, float tolerance_base, float tolerance_fraction, ple_lnum_t n_points, const ple_coord_t point_coords[], const int point_tag[], ple_lnum_t location[], float distance[])
{
	const GridObj* grid = static_cast<const GridObj*>(mesh);
	L_INFO("Number of points to locate " + std::to_string(n_points), GridUtils::logfile);
	for (int p = 0; p < n_points; p++)
	{
		double px, py, pz;
		// Set the point coordinates in the LUMA coordinate system
		px = point_coords[3 * p + 0] + L_PLE_OFFSET_X;
		py = point_coords[3 * p + 1] + L_PLE_OFFSET_Y;
		pz = point_coords[3 * p + 2] + L_PLE_OFFSET_Z;

		for (int face = 0; face < L_PLE_INTERFACES; face++)
		{
			// ! When pleName the first value is CS_inlet
			if (pleName[face].compare("CS_inlet")) // Result is same, return 0, else return 1
				continue;

			for (int ii = 0; ii < L_COUPLING_ENTITY; ii++)
			{
				double xdist, ydist, zdist;
				// Based on the properity, judge if CS return point is located on face
				xdist = abs(point_coords[3 * p] - plePosX[face][ii]) * (double) L_RESOLUTION;
				ydist = abs(point_coords[3 * p + 1] - plePosY[face][ii]) * (double) L_RESOLUTION;
				zdist = abs(point_coords[3 * p + 2] - plePosZ[face][ii]) * (double) L_RESOLUTION;
				
				if(pleOffset[ii] == eXneg)
					xdist -= 1.0;
				if(pleOffset[ii] == eYneg)
					ydist -= 1.0;
				if(pleOffset[ii] == eZneg)
					zdist -= 1.0;
				
				// Error are supposed as 0.1 lattcie
				if (xdist < 0.1)
				{
					switch (pleOffset[ii])
					{
						case eXpos:
							px += 0.3 / L_RESOLUTION;
							break;
						case eXneg:
							px -= 0.3 / L_RESOLUTION;
							break;
						default:
							break;
					}
				}

				else if (ydist < 0.1)
				{
					switch(pleOffset[ii])
					{
						case eYpos:
							py += 0.3 / L_RESOLUTION;
							break;
						case eYneg:
							py -= 0.3 / L_RESOLUTION;
							break;
						default:
							break;
					}
				}

				else if (zdist < 0.1)
				{
					switch(pleOffset[ii])
					{
						case eZpos:
							pz += 0.3 / L_RESOLUTION;
							break;
						case eZneg:
							pz -= 0.3 / L_RESOLUTION;
							break;
						default:
							break;
					}
			    }
			}
		}
		eLocationOnRank loc;
		std::vector<int> pos = { -1, -1, -1 };
		// Check that the point is not in the sending or receiving layers. (I think we don't want the point to be in a sending or receiving halo). 
		if (GridUtils::isOnThisRank(px, py, pz, &loc, grid, &pos))
		{
			// If the point is in the receiving halo, it is counted as unlocated. 
			if (loc == eHalo)
			{
				L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates " + std::to_string(px) +  " " + std::to_string(py) + " " + std::to_string(pz) + " is in the halo of rank " + std::to_string(GridUtils::safeGetRank()), GridUtils::logfile);
				location[p] = -1;
				distance[p] = -1;
				continue;
			}
			else if (loc == eCore) // Point is in the core of the rank's meshes, including the sending halo. 
			{
				// Pos is the local rank's cell index. I think it returns the index in the grid you pass it. And that is a global index in that grid
				// the indexes in pos go from 0 to Nx (where Nd are the total number of cells in that grid in the d direction, d = x,y,z).
				// ATTENTION! From what Yvan told me I understand that the indices in location should be local to the current rank but I'm not sure
				// how to get them with LUMA. I will leave it like this for the moment because I'm not sure LUMA will be able to work with local rank indices. 

				L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates " + std::to_string(px) + " " + std::to_string(py) + " " + std::to_string(pz) + " is inside the mesh of rank " + std::to_string(GridUtils::safeGetRank()), GridUtils::logfile);

				// Flatten the index returned by inOnThisRank
				int id = pos[2] + pos[1] * grid->K_lim + pos[0] * grid->K_lim * grid->M_lim;
				location[p] = id;

				// TODO: To fill the distance array I have to get the LUMA coordinates of id and then calculate the absolute distance between them
				// and the coordinates in point_coords. I think it is just XPos[pos[0]], YPos[pos[1]], ZPos[pos[2]]. Then divide it by dx and this should be smaller than 1
				distance[p] = GridUtils::vecnorm(px - grid->XPos.at(pos[0]),
					py - grid->YPos.at(pos[1]),
					pz - grid->ZPos.at(pos[2]));// / grid->dh;
				
				if (distance[p] > 1)
				{
					L_ERROR("The distance between PLE point with coordinates x = " + std::to_string(px) +
						" y = " + std::to_string(py) +
						" z = " + std::to_string(pz) +
						" and the LUMA mesh position x = " + std::to_string(grid->XPos.at(pos[0])) +
						" y = " + std::to_string(grid->YPos.at(pos[1])) +
						" z = " + std::to_string(grid->ZPos.at(pos[2])) +
						" is " + std::to_string(distance[p]) + "which is bigger than 1. But the point should be inside the cell",
						GridUtils::logfile);
					break;
				}
				// WARNING: The LUMA search functions are not prepared to work with tolerance. So I think I will set all the LUMA tolerances to 0. 
				// It should work. 
			}
			else
			{
				L_ERROR("Trying to locate a point from PLE mesh. The point is in the mesh but its location is eNone", GridUtils::logfile);
				break;
			}
		}
		else // If the point is not in the current rank, set the distance to -1 and the index of the point
		{
			location[p] = -1;
			distance[p] = -1;

			L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates "+ std::to_string(px) + " " + std::to_string(py) + " " + std::to_string(pz) + " is outside the mesh of rank " + std::to_string(GridUtils::safeGetRank()), GridUtils::logfile);
			continue;
		}

		L_INFO("LUMA: Point " + std::to_string(p) + " with coordinates " + std::to_string(px) + " " + std::to_string(py) + " " + std::to_string(pz) + " distance " + std::to_string(distance[p]), GridUtils::logfile);
	}
}

void PLEAdapter::exchangeMessage(const std::string &messageToSend, std::string *messageReceived)
{
	if (GridUtils::safeGetRank() < 1) {

		MPI_Status status;

		if (messageToSend.size() > 0) {

			if (messageToSend.size() > 32)
				L_WARN("The message to send to PLE is longer than 32 characters and it will be cut. Original message: " + messageToSend, GridUtils::logfile);

			char _op_name_send[33];
			strncpy(_op_name_send, messageToSend.c_str(), 32);
			_op_name_send[32] = '\0';

			/* Exchange command messages */
			if (messageReceived != NULL) 
			{
				char op_name_recv[32];
				MPI_Sendrecv(_op_name_send, 32, MPI_CHAR,
					CSRootRank_, csLumaCouplingTag_,
					op_name_recv, 32, MPI_CHAR,
					CSRootRank_, csLumaCouplingTag_,
					mpiCS_, &status);
				std::string temp(op_name_recv);
				messageReceived = &temp;
			}
			else
			{
				MPI_Send(_op_name_send, 32, MPI_CHAR,
					CSRootRank_, csLumaCouplingTag_,
					mpiCS_);
			}
		}
		else if (messageReceived != NULL) 
		{
			char op_name_recv[32];
			MPI_Recv(op_name_recv, 32, MPI_CHAR,
				CSRootRank_, csLumaCouplingTag_,
				mpiCS_, &status);
			std::string temp(op_name_recv);
			messageReceived = &temp;
		}
	}

	/*if ((messageReceived->size() != 0) && GridUtils::safeGetRank() > -1)
	{
		std::cout << "LUMA: I'm broadcasting op_name_recv and I have no idea why should I do that" << std::endl;
		char op_name_recv[32];
		strncpy(op_name_recv, messageReceived->c_str(), 32);
		MPI_Bcast(op_name_recv, 32, MPI_CHAR, 0, mpiCS_);
		op_name_recv[32] = '\0';
	}*/
}

bool PLEAdapter::configFileRead()
{
	//Read the data from definitions.h //
	
	participantName_ = L_PLE_PARTICIPANT_NAME;
	//std::string app_type = "LUMA";
	//std::string app_name = "LEFT"; // Name of the domain. It is in the mpi command line, or I can also write it in the config file?
								   // I think it is better to read from mpi command line so that there is only one config file. I don't like it but it is 
								   // what CS does.
	tolerance_ = 0.0; // LUMA isOnRank(...) doesn't work with tolerance.

	offset_.push_back(L_PLE_OFFSET_X);
	offset_.push_back(L_PLE_OFFSET_Y);
	offset_.push_back(L_PLE_OFFSET_Z);

	// NOTE: This might be different if using grid refinement!!
	double dh = LUMAGrid_->Grids->dh;
	// These variables are used to reduce the type transfer problem (double->int)
	double sizex, sizey, sizez;

	for (int i = 0; i < L_PLE_INTERFACES; i++)
	{
		InterfaceConfig pleInt;

		pleInt.meshName = pleName[i];

		for (int j = 0; j < L_COUPLING_ENTITY; j++)
		{
			pleInt.position.push_back(plePosX[i][j]);
			pleInt.position.push_back(plePosY[i][j]);
			pleInt.position.push_back(plePosZ[i][j]);

			sizex = pleSizeX[i][j] * (1.0 / dh);
			sizey = pleSizeY[i][j] * (1.0 / dh);
			sizez = pleSizeZ[i][j] * (1.0 / dh);

			pleInt.dimensions.push_back(static_cast<int>(sizex));
			pleInt.dimensions.push_back(static_cast<int>(sizey));
			pleInt.dimensions.push_back(static_cast<int>(sizez));
		}
		std::cout << " Read data " << pleRead[i] << " Write data " << pleWrite[i] << std::endl;
		
		for (char const c : pleRead[i])
		{
			pleInt.readData.push_back(c);
		}
		for (int c = 0; c < pleWrite[i].size(); c++)
		{
			pleInt.writeData.push_back(pleWrite[i][c]);
		}

		/* REMEMBER: to correct this logfile term.pleSizeX[i]*/
		// L_INFO("Interface " + std::to_string(i) + " dx= " + std::to_string(dh) + " Size x " + std::to_string(pleSizeX[i]) + " Nx= " + std::to_string(pleInt.dimensions.at(0)), GridUtils::logfile);

		interfacesConfig_.push_back(pleInt);
	}
    return true;
}

void PLEAdapter::addPLELocator(int i)
{
	//std::cout << "LUMA: Hello before PLE locator creation. " << std::endl;

	if ((i == (coordinates_.size() - 1)) && (i == locators_.size()))  // If this locators_ position doesn't contain any data yet but has a corresponding coordinates.
	{
		L_INFO("Creating the PLE locator number " + std::to_string(i), GridUtils::logfile);

		//-  Create PLE locator
		locators_.push_back(ple_locator_create(mpiCS_, CSNumRanks_, CSRootRank_));
	}
	else if ((i < coordinates_.size()) && (i < locators_.size()))  // If this locators_ position exists
	{
		L_WARN("Regenerating the ple_locator in position " + std::to_string(i), GridUtils::logfile);
		locators_.at(i) = ple_locator_create(mpiCS_, CSNumRanks_, CSRootRank_);
	}
	else
	{
		L_ERROR("PLE locator number " + std::to_string(i) + " can't be created without associated coordinates. Please initialise the coordinates vector before initialising flow data", GridUtils::logfile);
	}

	//std::cout << "LUMA: PLE locator created. " << std::endl;

	int locator_options[PLE_LOCATOR_N_OPTIONS];
	locator_options[PLE_LOCATOR_NUMBERING] = 1;


	// NOTE: I'm not sure if I'll need this vector outside this function. So I'll create it as a local variable for the moment. 
	std::vector<float> lumaToCSdist;
	lumaToCSdist.resize(coordinates_.at(i).size() / L_DIMS);
	//- Call PLE locator set mesh
	ple_locator_set_mesh(locators_.at(i),
		LUMAGrid_->Grids,
		locator_options,
		0.,
		tolerance_,
		L_DIMS,
		coordinates_.at(i).size() / L_DIMS,
		NULL,
		NULL,
		coordinates_.at(i).data(),
		lumaToCSdist.data(),
		&meshExtents,
		&pointInMeshMapping);
	
	for (int l = 0; l < coordinates_.at(i).size() / L_DIMS; l++)
	{
		L_INFO("LUMA to CS dist " + std::to_string(lumaToCSdist.at(l)), GridUtils::logfile);
	}
	
	std::cout << "*-*-*-*-*-I am testing the location from the local mesh*-*-*-*-*-*-*-*-*" << std::endl;
	for (int ii = 0; ii < coordinates_.at(i).size()/L_DIMS;ii++)
	{
		std::cout << coordinates_.at(i)[ii] << " ";
	}

	std::cout << "*-*-*-*-*-I am testing the location form the distant mesh*-*-*-*-*-*-*-*-*" << std::endl;
	size_t nCoupledPoints = ple_locator_get_n_dist_points(locators_.at(i));
	const int* IDCoupledPoints = ple_locator_get_dist_locations(locators_.at(i));
	for (int ii = 0; ii < nCoupledPoints;ii++)
	{
		std::cout << IDCoupledPoints[ii] << " ";
	}
	
	// This is used to obtaine the neighbour position to interpolate
	if (!pleName[i].compare("CS_inlet"))	// compare function: If is CS_inlet, return false
	{
		//const GridObj* grid;
		GridObj* currentGrid = LUMAGrid_->Grids;
		int num_xdiff = currentGrid->M_lim * currentGrid->K_lim;
		int num_ydiff = currentGrid->K_lim;
		int num_zdiff =  1;

		size_t nCoupledPoints = ple_locator_get_n_dist_points(locators_.at(i));
		const int* IDCoupledPoints = ple_locator_get_dist_locations(locators_.at(i));

		std::vector<int> send_id_diff(IDCoupledPoints, IDCoupledPoints + nCoupledPoints);

		for (int ii = 0; ii < nCoupledPoints; ii++)
		{
		//std::cout << positionID_.at(IDCoupledPoints[ii])<<std::endl;

			switch(positionID_.at(IDCoupledPoints[ii]))
			{
				//std::cout << positionID_.at(IDCoupledPoints[ii]) << " ";
				case eXneg:
					send_id_diff[ii] += num_xdiff;
					break;
				case eXpos:
					send_id_diff[ii] -= num_xdiff;
					break;
				case eYneg:
					send_id_diff[ii] += num_ydiff;
					break;
				case eYpos:
					send_id_diff[ii] -= num_ydiff;
					break;
				case eZneg:
					send_id_diff[ii] += num_zdiff;
					break;
				case eZpos:
					send_id_diff[ii] -= num_zdiff;
					break;
				default:
					break;
			}
		}
		send_id_diff_.resize(nCoupledPoints);
		send_id_diff_ = send_id_diff;
	}
	//- Check that all points are effectively located
	ple_lnum_t nExterior = ple_locator_get_n_exterior(locators_.at(i));

	L_INFO("Number of points not located in LUMA:  " + std::to_string(nExterior), GridUtils::logfile);
	//const ple_lnum_t* extPoints = ple_locator_get_exterior_list(locators_.at(i));
	//for (int np = 0; np < nExterior; np++)
	//	std::cout << "My rank: " << GridUtils::safeGetRank() << " Point not located " << extPoints[np] << std::endl;

	//- Check if CS has also located all points. 
	std::string messageToSend, messageReceived;
	if (nExterior > 0)
		messageToSend = "coupling:location:incomplete";
	else
		messageToSend = "coupling:location:ok";

	exchangeMessage(messageToSend, &messageReceived);
	if (messageReceived.find("coupling:location:incomplete") != std::string::npos)
	{
		L_WARN("Coupling with Code_Saturne impossible:\n" + std::to_string(nExterior) + " points from mesh " +
			interfacesConfig_.at(i).meshName + " not located on Code_Saturne mesh.", GridUtils::logfile);

		// Stop coupling
		setSyncFlag(PLE_COUPLING_STOP);

	}

	// If the CS mesh is a boundary, CS wants to know the luma_to_cs_dist
	if (csBoundary[i])
	{
		ple_locator_exchange_point_var(locators_.at(i),
			NULL,
			lumaToCSdist.data(),
			NULL,
			sizeof(float),
			1,
			1);
	}
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define an optional synchronization flag for external couplings.
 *
 * This flag is used by all couplings based on the PLE (Parallel Location
 * and Exchange) group synchronization mechanism, which include couplings
 * with SYRTHES, LUMA, code_saturne, and neptune_cfd.
 *
 * It is defined by a mask, so for example flags f1, f2, and f3 may be
 * combined using the "f1 | f2 | f2" syntax.
 *
 * Note also that for Code_Saturne, in the case of a variable time step,
 * the reference time step is synchronized at the beginning of each
 * iteration, but the actual time step is recomputed later.
 *
 * Possible flags are:
 *  - PLE_COUPLING_TS_MIN        Use smallest time step
 *  - PLE_COUPLING_TS_LEADER     Prescribe time step for the group
 *                               (only one member may set this flag)
 *  - PLE_COUPLING_UNSTEADY      Inform others that this instance is
 *                               using an unsteady solution approach
 *  - PLE_COUPLING_STEADY        Inform others that this instance is
 *                               using a steady solution approach
 *  - PLE_COUPLING_USER_1        User definable flag
 *  - PLE_COUPLING_USER_2        User definable flag
 *  - PLE_COUPLING_USER_3        User definable flag
 *  - PLE_COUPLING_USER_4        User definable flag
 *
 * To force stopping, PLE_COUPLING_STOP may be set. In this case,
 * the calculation will stop at the first synchronization, even if
 * this function is called again with another flag.
 *
 * \param[in]  flag  synchronization flag to apply to couplings
 */
 /*----------------------------------------------------------------------------*/
void PLEAdapter::setSyncFlag(int flag)
{
	int stop_mask = pleCouplingFlag_ & PLE_COUPLING_STOP;

	pleCouplingFlag_ = flag | stop_mask;
}

void PLEAdapter::sendData()
{
	//NOTE: The data is not interpolated before being sent, is interpolated only on receive. 
	// LUMA sends the data from the LUMA cells that are closest to the CS cells. 
	// When CS receives the data, it interpolates it to the actual positions of the CS cells because CS has the distance between each
	// CS and LUMA cell. 
	
	// For all the interfaces
	for (int i = 0; i < interfacesConfig_.size(); i++)
	{
		// Check which ones have data to write to PLE
		for (auto writeJ = interfacesConfig_.at(i).writeData.begin(); writeJ != interfacesConfig_.at(i).writeData.end(); ++writeJ)
		{
			// Get the number of coupled points and their CS cell/face id
			size_t nCoupledPoints_dist = ple_locator_get_n_dist_points(locators_.at(i));
			const int* IDCoupledPoints_dist = ple_locator_get_dist_locations(locators_.at(i));
			// TODO: Check for nested grids and return the correct grid to get the LUMA variables from.
			// pointInMesh should do the check for nested grids too to store the correct indices in IDCoupledPoints.
			// I think there is a function in GridUtils that returns the finest grid in a position.
			// But what if the coupled surface takes more than one grid level?
			// Then I'll have to get the grid for each point in the surface...
			GridObj* currentGrid = LUMAGrid_->Grids;

			if ((*writeJ == 't') || (*writeJ == 'T'))
			{
				// Extract temperature from the LUMA grid. The IDs in IDCoupledPoints should be safe because are given by PLEAdapter::pointInMesh()
				// So it is OK to directly access the temperature array.
				std::vector<double> send_v(nCoupledPoints_dist * L_DIMS, 0.0);
				std::vector<int> send_id(IDCoupledPoints_dist, IDCoupledPoints_dist + nCoupledPoints_dist);

				currentGrid->coupling_extractData("temperature", send_id, send_id_diff_,send_v);
				
				//I think I should add the unit conversion from lattice units(used in LUMA) to physical units(used in CS)
				GridUnits::tlat2phys(send_v);

				ple_locator_exchange_point_var(locators_.at(i), send_v.data(), NULL, NULL, sizeof(double), 1, 0);				
			}
			else if ((*writeJ == 'r') || (*writeJ == 'R'))
			{
				// Extract velocity from the LUMA grid. The IDs in IDCoupledPoints should be safe because are given by PLEAdapter::pointInMesh()
				std::vector<double> send_v(nCoupledPoints_dist, 0.0);
				std::vector<int> send_id(IDCoupledPoints_dist, IDCoupledPoints_dist + nCoupledPoints_dist);

				currentGrid->coupling_extractData("rho", send_id, send_id_diff_,send_v);

				ple_locator_exchange_point_var(locators_.at(i), send_v.data(), NULL, NULL, sizeof(double), 1, 0);

				// TODO: Implement unity conversion from LUMA density to Code_Saturne pressure 

			}
			else if ((*writeJ == 'v') || (*writeJ == 'V'))
			{
				// Extract velocity from the LUMA grid. The IDs in IDCoupledPoints should be safe because are given by PLEAdapter::pointInMesh()
				std::vector<double> send_v(nCoupledPoints_dist * L_DIMS, 0.0);
				std::vector<int> send_id(IDCoupledPoints_dist, IDCoupledPoints_dist + nCoupledPoints_dist);

				currentGrid->coupling_extractData("v", send_id, send_id_diff_, send_v); //send_id_diff_

				// Convert send_v from LUMA units to CS units. 
				GridUnits::ulbm2ud(send_v, currentGrid);
				ple_locator_exchange_point_var(locators_.at(i), send_v.data(), NULL, NULL, sizeof(double), L_DIMS, 0);
			}
			else
			{
				L_ERROR("Data type for " + std::to_string(*writeJ) + " not recognised. The data to transfer to/from PLE has to contain v or V for velocity and t or T for temperature.", GridUtils::logfile);
			}
		}
	}	
	// Extract the data to send from the LUMA grid. 
}

void PLEAdapter::receiveData()
{

	// For all the interfaces
	for (int i = 0; i < interfacesConfig_.size(); i++)
	{
		L_INFO("LUMA: interface " + std::to_string(i) , GridUtils::logfile);

		// Check which ones have data to read from PLE
		for (auto readJ = interfacesConfig_.at(i).readData.begin(); readJ != interfacesConfig_.at(i).readData.end(); ++readJ)
		{
			L_INFO("LUMA: interface " + std::to_string(i) + " has data from PLE.", GridUtils::logfile);

			// Get the number of coupled points and their LUMA cell id
			//size_t nCoupledPoints = ple_locator_get_n_dist_points(locators_.at(i));
			//const int* IDCoupledPoints = ple_locator_get_dist_locations(locators_.at(i));	
			size_t nCoupledPoints = ids_.at(i).size();

			L_INFO("Number of interior points " + std::to_string(nCoupledPoints), GridUtils::logfile);

			const int* IDCoupledPoints = ids_.at(i).data();
			std::cout << nCoupledPoints << " " << *IDCoupledPoints << " " << std::endl;
			//L_INFO("Number of interior points " + std::to_string(nCoupledPoints) + " ID of second point " + std::to_string(IDCoupledPoints[1]), GridUtils::logfile);

			// TODO: Check for nested grids and return the correct grid to get the LUMA variables from. 
			// pointInMesh should do the check for nested grids too to store the correct indices in IDCoupledPoints. 
			// I think there is a function in GridUtils that returns the finest grid in a position.
			// But what if the coupled surface takes more than one grid level?
			// Then I'll have to get the grid for each point in the surface...
			GridObj* currentGrid = LUMAGrid_->Grids;

			// Receive the data from CS
			if ((*readJ == 't') || (*readJ == 'T'))
			{
				// Extract temperature to the LUMA grid. The IDs in IDCoupledPoints should be safe because are given by PLEAdapter::pointInMesh()
				// So it is OK to directly access the temperature array.
				std::vector<double> received_t(nCoupledPoints, 0.0);
				std::vector<int> received_id(IDCoupledPoints, IDCoupledPoints + nCoupledPoints);
				
				//std::cout << "LUMA: Ready to receive temperature. "  << nCoupledPoints << " "  << received_t[nCoupledPoints*L_DIMS-1] << " " << received_t.size() << std::endl;

				// Receive the data from CS
				ple_locator_exchange_point_var(locators_.at(i), NULL, received_t.data(), NULL, sizeof(double), 1, 0);

				//Covert received_t from CS units to LUMA units

				GridUnits::tphys2lat(received_t);

				// Introduce the received temperature to the LUMA grid.
				currentGrid->coupling_addData("temperature", received_id, received_t);
			}
			else if ((*readJ == 'r') || (*readJ == 'R'))
			{
				// Create containers for the received data
				std::vector<double> received_r(nCoupledPoints, 0.0);
				std::vector<int> received_id(IDCoupledPoints, IDCoupledPoints + nCoupledPoints);

				// Receive the data from CS
				ple_locator_exchange_point_var(locators_.at(i), NULL, received_r.data(), NULL, sizeof(double), 1, 0);
				
				// Introduce the received density to the LUMA grid. 
				// The IDs in IDCoupledPoints should be safe because are given by PLEAdapter::pointInMesh()
				// For the moment the interpolation is just nearest neighbour. 
				currentGrid->coupling_addData("rho", received_id, received_r);
		

				// TODO: Implement unity conversion from LUMA density to Code_Saturne pressure 

			}
			else if ((*readJ == 'v') || (*readJ == 'V'))
			{
				// Create containers for the received data
				std::vector<double> received_v(nCoupledPoints * L_DIMS, 0.0);
				std::vector<int> received_id(IDCoupledPoints, IDCoupledPoints + nCoupledPoints);

				// Receive the data from CS
				ple_locator_exchange_point_var(locators_.at(i), NULL, received_v.data(), NULL, sizeof(double), 3, 0);

				// Convert received_v from CS units to LUMA units.
				GridUnits::ud2ulbm(received_v, currentGrid);
				// Introduce the received velocity to the LUMA grid. 
				// The IDs in IDCoupledPoints should be safe because are given by PLEAdapter::pointInMesh()
				currentGrid->coupling_addData("v", received_id, received_v);
				
				//std::cout << "LUMA: Velocity added to LUMA mesh." << std::endl;
			}
			else
			{
				L_ERROR("Data type for " + std::to_string(*readJ) + " not recognised. The data to transfer to/from PLE has to contain v or V for velocity and t or T for temperature.", GridUtils::logfile);
			}
		}
	}
}

bool PLEAdapter::configure()
{
    //Create a PLE locator for each PLE interface in the configuration file. 
	// ple_locator_set_mesh. Do I need PLEInterface? Or just a vector of ple_locators? 
	// I'll start with just a vector of PLE locators

	// Read PLE adapter configuration file.  
	configFileRead(); 
	std::cout << "I've read the configuration file. " << std::endl;

	pleSets_ = ple_coupling_mpi_set_create(pleCouplingFlag_, "LUMA", "LEFT", MPI_COMM_WORLD, lumaMpi_->ple_comm);  //or lumaMpi_->world_comm?

	numApps_ = ple_coupling_mpi_set_n_apps(pleSets_);
	appID_ = ple_coupling_mpi_set_get_app_id(pleSets_);
	//ple_coupling_mpi_set_info_t me = ple_coupling_mpi_set_get_info(pleSets_, appID_);
	//lumaRootRank_ = me.root_rank;
	//lumaNumRanks_ = me.n_ranks;
	//if(GridUtils::safeGetRank() < 1)
	//	L_INFO("LUMA root rank: " + std::to_string(lumaRootRank_) + ". Number of ranks for LUMA: " + std::to_string(lumaNumRanks_), GridUtils::logfile);

	int local_range[2] = { -1, -1 };
	int distant_range[2] = { -1, -1 };

	// Search for the MPI ranks of CS. 
	// Stop if you find more than one CS instance in the MPI call. 
	int CSNum = 0;
	for (int i = 0; i < numApps_; i++)
	{
		ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(pleSets_, i);

		if (std::string(ai.app_type).find("Code_Saturne") != std::string::npos)
		{
			if (CSNum < 1)
			{ 
				CSRootRank_ = ai.root_rank;
				CSNumRanks_ = ai.n_ranks;
				CSAppID_ = i;
				CSNum++;
				L_INFO("CS root rank: " + std::to_string(CSRootRank_) + ". Number of ranks for CS: " + std::to_string(CSNumRanks_), GridUtils::logfile);
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
	int i = 0;
	for (int i = 0; i < interfacesConfig_.size(); i++)
	{

		// Set coordinates in the global coordinate system (same coord system as Code_Saturne)
		setCoordinates(interfacesConfig_.at(i).dimensions, interfacesConfig_.at(i).position, LUMAGrid_->Grids->dh, i);

		/*exchangeData_->addRepo(interfacesConfig_.at(i).dimensions.at(0), interfacesConfig_.at(i).dimensions.at(1), interfacesConfig_.at(i).dimensions.at(2),
			interfacesConfig_.at(i).position.at(0), interfacesConfig_.at(i).position.at(1), interfacesConfig_.at(i).position.at(2),
			LUMAGrid_->Grids->dh);*/

		L_INFO("I'm creating the locators!", GridUtils::logfile);

		// Create space for the data to be written / read for the current repo. 
		for (auto readJ = interfacesConfig_.at(i).readData.begin(); readJ != interfacesConfig_.at(i).readData.end(); ++readJ)
		{
			if ((*readJ == 'v') || (*readJ == 'V'))
			{
				std::cout << "LUMA: " << std::to_string(*readJ) << std::endl;
				initialiseVectorData("v", 0.0, i);
			}
			else if ((*readJ == 't') || (*readJ == 'T'))
			{
				//exchangeData_->getRepo(i).initialiseScalarData(*readJ);
				initialiseScalarData("t", 0.0, i);
			}
			else
			{
				L_ERROR("Data type for " + std::to_string(*readJ) + " not recognised. Data in .yml configuration file has to contain v or V for velocity and t or T for temperature.", GridUtils::logfile);
			}
		}
		for (auto writeJ = interfacesConfig_.at(i).writeData.begin(); writeJ != interfacesConfig_.at(i).writeData.end(); ++writeJ)
		{
			if ((*writeJ == 'v') || (*writeJ == 'V'))
			{
				//exchangeData_->getRepo(i).initialiseVectorData(*writeJ);
				L_INFO("Ready to intialise velocity data", GridUtils::logfile);
				initialiseVectorData("v", 0.0, i);
			}
			else if ((*writeJ == 't') || (*writeJ == 'T'))
			{
				//exchangeData_->getRepo(i).initialiseScalarData(*writeJ);
				initialiseScalarData("t", 0.0, i);
			}
			else
			{
				L_ERROR("Data type for " + std::to_string(*writeJ) + " not recognised. Data in .yml configuration file has to contain v or V for velocity and t or T for temperature.", GridUtils::logfile);
			}
		}

		L_INFO("Creating PLE adapter..", GridUtils::logfile);

		// Create and configure the PLE locator
		addPLELocator(i);

		//L_INFO("Interface created on mesh" + interfacesConfig_.at(i).meshName, GridUtils::logfile);

		// Synchronise with CS after location
		std::string messageToSend, messageReceived;
		messageToSend = "coupling:start";

		exchangeMessage(messageToSend, &messageReceived);
		if (messageReceived.find("coupling:error:location") != std::string::npos)
		{
			L_WARN("Messge received from Code_Saturne: " + messageReceived + " indicates meshes have not been matche correctly. \n " +
				" The calculation will not run.", GridUtils::logfile);

			// Stop coupling
			setSyncFlag(PLE_COUPLING_STOP);

		}
		else
		{
			L_INFO(" Mesh located for Code_Saturne mesh " + interfacesConfig_.at(i).meshName, GridUtils::logfile);
		}

	}





  return true;
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

void PLEAdapter::setCoordinates(std::vector<int> Nxyz, std::vector<double> xyz0, double dx, int inter_id)
{
	/*if (numDataPoints_ > 0)
	{
		L_ERROR("This InOutRepo already contains a coordinates vector", GridUtils::logfile);
	}*/

	L_INFO("I'm in setCoordinates", GridUtils::logfile);

	std::vector<double> coord;
	std::vector<int> id;
	std::vector<double> temp_coord;
	std::vector<int> location;
	std::vector<float> dist;
	std::vector<int> dummy;
	temp_coord.resize(3);
	location.resize(1);
	dist.resize(1);
	dummy.resize(1);
	bool incells;

	for (int en = 0; en < L_COUPLING_ENTITY; en++)
	{
		L_INFO("Nx[" + std::to_string(en) + "]=" + std::to_string(Nxyz.at(L_DIMS*en+0)) + " Ny["+ std::to_string(en) +"]=" + std::to_string(Nxyz.at(L_DIMS*en+1)) 
				+ " Nz["+ std::to_string(en) + "]=" + std::to_string(Nxyz.at(L_DIMS*en+2)), GridUtils::logfile);


		// Fill in the coordinates vector. C style ordering
		for (int i = 0; i < Nxyz.at(L_DIMS*en+0); i++)
		{
			for (int j = 0; j < Nxyz.at(L_DIMS*en+1); j++)
			{
				for (int k = 0; k < Nxyz.at(L_DIMS*en+2); k++)
				{
					// Define on the cell/face center
					temp_coord[0] = xyz0.at(L_DIMS * en + 0) + i * dx + dx / static_cast<double>(2.0);
					temp_coord[1] = xyz0.at(L_DIMS * en + 1) + j * dx + dx / static_cast<double>(2.0);
					temp_coord[2] = xyz0.at(L_DIMS * en + 2) + k * dx + dx / static_cast<double>(2.0);

					// Check if temp_coord is in the current rank and store its LUMA mesh ID.
					incells= pointInMesh(LUMAGrid_->Grids, 0, 0, 1, temp_coord.data(), dummy.data(), location.data(), dist.data(), inter_id);
					if(!incells)
						continue;
					
					// Create the id and position mapping
					if (!pleName[inter_id].compare("CS_inlet"))
						positionID_[location.at(0)] = pleOffset[en];

					// Add the coordinates and the ID to the coordinates and ID arrays if they are in the mesh for the current rank.
					if (location.at(0) != -1)
					{
						for (int d = 0; d < L_DIMS; d++)
							coord.push_back(temp_coord.at(d));
						id.push_back(location.at(0));

						L_INFO("Cell ID " + std::to_string(id.back()), GridUtils::logfile);
					}
					L_INFO("Coord z: " + std::to_string(coord.back()), GridUtils::logfile);
				}
			}
		}
	}

	if (inter_id == coordinates_.size())
	{
		coordinates_.push_back(coord);
		ids_.push_back(id);
	}
	else if (inter_id < coordinates_.size())
	{
		L_INFO("Updating the PLE coupling coordinates for position " + std::to_string(inter_id), GridUtils::logfile);
		coordinates_.at(inter_id) = coord;
		ids_.at(inter_id) = id;
	}
	else
	{
		L_ERROR("The PLE locator vector index " + std::to_string(inter_id) + " is bigger than the size of the PLE locator vector. ", GridUtils::logfile);
	}

	L_INFO("LUMA mesh size " + std::to_string(coordinates_.at(inter_id).size()), GridUtils::logfile);

}

void PLEAdapter::initialiseVectorData(std::string name, double initValue, int i)
{
	L_INFO("Size of coordinate vector " + std::to_string(coordinates_.at(i).size()), GridUtils::logfile);
	
	// Allocate new arrays
	std::map<std::string, std::vector<double>> tempMap;
	std::vector<double> temp;
	temp.resize(coordinates_.at(i).size(),initValue);
	tempMap[name] = temp;

	if ((i == (coordinates_.size() - 1)) && (i == vectorData_.size()))  // If this vectorData_ position doesn't contain any data yet but has a corresponding coordinates.
	{
		L_INFO("Creating a new vectorData_storage for " + name + " at position " + std::to_string(i), GridUtils::logfile);
		vectorData_.push_back(tempMap);
	}
	else if ((i < coordinates_.size()) && (i < vectorData_.size()))  // If this vectorData_ position exists
	{
		L_INFO("Adding " + name + " to vectorData storage at position " + std::to_string(i), GridUtils::logfile);
		vectorData_.at(i)[name] = temp;
	}
	else   // If this vectorData_ position is bigger than the size of the coordinates_ and vectorData_ vectors. 
	{
		L_ERROR("This vectorData doesn't have associated coordinates. Please initialise the coordinates vector before initialising flow data", GridUtils::logfile);
	}

	std::cout << "LUMA: HI! numDataPoints = " << coordinates_.at(i).size() << std::endl;
}

void PLEAdapter::initialiseScalarData(std::string name, double initValue, int i)
{
	// Allocate new arrays
	std::map<std::string, std::vector<double>> tempMap;
	std::vector<double> temp;
	temp.resize(coordinates_.at(i).size() / 3.0, initValue);
	tempMap[name] = temp;

	if ((i == (coordinates_.size() - 1)) && (i == scalarData_.size()))  // If this scalarData_ position doesn't contain any data yet but has a corresponding coordinates.
	{
		L_INFO("Creating a new scalarData_storage for " + name + " at position " + std::to_string(i), GridUtils::logfile);
		scalarData_.push_back(tempMap);
	}
	else if ((i < coordinates_.size()) && (i < scalarData_.size()))  // If this scalarData_ position exists
	{
		L_INFO("Adding " + name + " to scalarData storage at position " + std::to_string(i), GridUtils::logfile);
		scalarData_.at(i)[name] = temp;
	}
	else   // If this scalarData_ position is bigger than the size of the coordinates_ and vectorData_ vectors. 
	{
		L_ERROR("This scalarData_ doesn't have associated coordinates. Please initialise the coordinates vector before initialising flow data", GridUtils::logfile);
	}

	std::cout << "LUMA: HI! numDataPoints = " << coordinates_.at(i).size() << std::endl;
}

PLEAdapter::~PLEAdapter()
{
    teardown();

    return;
}
