/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) 2015, 2016
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * distribution without written consent.
 *
 */

#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"
#include "../inc/PCpts.h"
#include "../inc/GridObj.h"

// ************************************************************************** //
/// \brief Write out position of immersed boundary bodies.
/// \param	timestep	timestep at which the write out is being performed.
void ObjectManager::io_writeBodyPosition(int timestep) {

	int rank = GridUtils::safeGetRank();

	for (size_t ib = 0; ib < iBody.size(); ib++) {


			// Open file for given time step
			std::ofstream jout;
			jout.open(GridUtils::path_str + "/Body_" + std::to_string(ib) + "_position_" + std::to_string(timestep) + 
				"_rank" + std::to_string(rank) + ".out", std::ios::out);
			jout << "x" + std::to_string(timestep) + ", y" + std::to_string(timestep) + ", z" << std::endl;

			// Write out position
			for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
#if (L_DIMS == 3)
				jout	<< iBody[ib].markers[i].position[0] << ", " 
						<< iBody[ib].markers[i].position[1] << ", " 
						<< iBody[ib].markers[i].position[2] << std::endl;
#else
				jout	<< iBody[ib].markers[i].position[0] << ", " 
						<< iBody[ib].markers[i].position[1] << ", " 
						<< 0.0 << std::endl;
#endif
			}
			jout.close();

	}

}


// ************************************************************************** //
/// \brief Write out forces on the markers of immersed boundary bodies.
/// \param	timestep	timestep at which the write out is being performed.
void ObjectManager::io_writeLiftDrag(int timestep) {

	int rank = GridUtils::safeGetRank();

	for (size_t ib = 0; ib < iBody.size(); ib++) {


			// Open file for given time step
			std::ofstream jout;
			jout.open(GridUtils::path_str + "/Body_" + std::to_string(ib) + "_LD_" + std::to_string(timestep) + "_rank" + std::to_string(rank) + ".out", std::ios::out);
			jout << "L" + std::to_string(timestep) + ", D" + std::to_string(timestep) << std::endl;

			// Sum variables
			double Lsum = 0.0, Dsum = 0.0;

			// Compute lift and drag
			for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
				jout	<< iBody[ib].markers[i].force_xyz[0] << ", " 
						<< iBody[ib].markers[i].force_xyz[1] << std::endl;
				Lsum += iBody[ib].markers[i].force_xyz[0];
				Dsum += iBody[ib].markers[i].force_xyz[1];
			}

			jout << "Totals = " << std::endl;
			jout << Lsum << ", " << Dsum << std::endl;
			jout.close();

	}

}

// ************************************************************************** //
/// \brief	Read/write IB body information to restart file.
/// \param	IO_flag	flag indicating write (true) or read (false).
/// \param	level	level of the grid begin written/read
void ObjectManager::io_restart(eIOFlag IO_flag, int level) {

	int rank = GridUtils::safeGetRank();

	if (IO_flag == eWrite) {

		// Output stream
		std::ofstream file;

		if (level == 0) {
			// Overwrite as first to write
			file.open(GridUtils::path_str + "/restart_IBBody_Rnk" + std::to_string(rank) + ".out", std::ios::out);
		} else if (level == 0) {
			// Append
			file.open(GridUtils::path_str + "/restart_IBBody_Rnk" + std::to_string(rank) + ".out", std::ios::out | std::ios::app);
		} else {
			// Must be a subgrid which doesn't own any IB-bodies so return
			return;
		}


		// Counters
		size_t b, m, num_bod = iBody.size();

		// Write out the number of bodies
		file << num_bod << std::endl;

		// Loop over bodies and markers and write out the positions
		for (b = 0; b < num_bod; b++) {

			// Add a separator between bodies
			file << "\t/\t";

			// Number of markers with separator
			file << iBody[b].markers.size() << "\t/\t";

			for (m = 0; m < iBody[b].markers.size(); m++) {

				// Positions of each marker
				file	<< iBody[b].markers[m].position[0] << "\t"
						<< iBody[b].markers[m].position[1] << "\t"
						<< iBody[b].markers[m].position[2] << "\t";


				if (iBody[b].isFlexible) {

					// Old positions of each marker
					file << iBody[b].markers[m].position_old[0] << "\t"
						<< iBody[b].markers[m].position_old[1] << "\t"
						<< iBody[b].markers[m].position_old[2] << "\t";
				}

			}

		}

		// Close file
		file.close();


	} else {

		// Input stream
		std::ifstream file;

		// We only enter this routine if on correct level so no need to check
		file.open("./input/restart_IBBody_Rnk" + std::to_string(rank) + ".out", std::ios::in);

		if (!file.is_open()) {
			L_ERROR("Error opening IBM restart file. Exiting.", GridUtils::logfile);
		}

		// Read in one line of file at a time
		std::string line_in;	// String to store line
		std::istringstream iss;	// Buffer stream

		// Get line up to separator and put in buffer
		std::getline(file,line_in,'/');
		iss.str(line_in);
		iss.seekg(0); // Reset buffer position to start of buffer

		// Counters
		int b, m, num_bod, num_mark;

		// Read in number of bodies
		iss >> num_bod;

		// Check number of bodies is correct
		if (iBody.size() != num_bod) {
			L_ERROR("Number of IBM bodies does not match the number specified in the restart file. Exiting.", GridUtils::logfile);
		}

		// Loop over bodies
		for (b = 0; b < num_bod; b++) {

			// Exit if reached end of file
			if (file.eof()) break;

			// Get next bit up to separator and put in buffer
			std::getline(file,line_in,'/');
			iss.str(line_in);
			iss.seekg(0); // Reset buffer position to start of buffer

			// Read in the number of markers for this body
			iss >> num_mark;

			// Check number of markers the same
			if (iBody[b].markers.size() != num_mark) {
				L_ERROR("Number of IBM markers does not match the number specified for body " +
					std::to_string(b) + " in the restart file. Exiting.", GridUtils::logfile);
			}

			// Read in marker data
			for (m = 0; m < num_mark; m++) {

				// Positions of each marker
				iss		>> iBody[b].markers[m].position[0]
						>> iBody[b].markers[m].position[1]
						>> iBody[b].markers[m].position[2];

				if (iBody[b].isFlexible) {

					// Old positions of each marker
					iss		>> iBody[b].markers[m].position_old[0]
							>> iBody[b].markers[m].position_old[1]
							>> iBody[b].markers[m].position_old[2];
				}
			}

		}

		// Close file
		file.close();

	}

}


// ************************************************************************** //
/// \brief	Write IB body data to VTK file.
///
///			Currently can only write out un-closed bodies like filaments.
///
/// \param	tval	time value at which the write out is being performed.
void ObjectManager::io_vtkIBBWriter(double tval) {

	// Get the rank
	int rank = GridUtils::safeGetRank();

    // Loop through each iBody
    for (size_t ib = 0; ib < iBody.size(); ib++) {

        // Create file name then output file stream
        std::stringstream fileName;
        fileName << GridUtils::path_str + "/vtk_IBout.Body" << iBody[ib].id << "." << std::to_string(rank) << "." << (int)tval << ".vtk";

        std::ofstream fout;
        fout.open( fileName.str().c_str() );

        // Add header information
        fout << "# vtk DataFile Version 3.0f\n";
        fout << "IB Output for body ID " << iBody[ib].id << ", rank " << std::to_string(rank) << " at time t = " << (int)tval << "\n";
        fout << "ASCII\n";
        fout << "DATASET POLYDATA\n";


        // Write out the positions of each Lagrange marker
        fout << "POINTS " << iBody[ib].markers.size() << " float\n";
        for (size_t i = 0; i < iBody[ib].markers.size(); i++) {

        	fout	<< iBody[ib].markers[i].position[0] << " "
					<< iBody[ib].markers[i].position[1] << " "
					<< iBody[ib].markers[i].position[2] << std::endl;
        }


        // Write out the connectivity of each Lagrange marker
        size_t nLines = iBody[ib].markers.size() - 1;

        if (iBody[ib].closed_surface == false)
            fout << "LINES " << nLines << " " << 3 * nLines << std::endl;
        else if (iBody[ib].closed_surface == true)
            fout << "LINES " << nLines + 1 << " " << 3 * (nLines + 1) << std::endl;

        for (size_t i = 0; i < nLines; i++) {
            fout << 2 << " " << i << " " << i + 1 << std::endl;
        }

        // If iBody[ib] is a closed surface then join last point to first point
        if (iBody[ib].closed_surface == true) {
            fout << 2 << " " << nLines << " " << 0 << std::endl;
        }

        fout.close();
    }
}

// ************************************************************************** //
/// \brief	Read in geometry config file.
///
///			Input data must be in correct format as specified in documentation.
///
void ObjectManager::io_readInGeomConfig() {

	// Open config file
	std::ifstream file;
	file.open("./input/geometry.config", std::ios::in);

	// Handle failure to open
	if (!file.is_open()) {
		L_ERROR("Error opening geometry configuration file. Exiting.", GridUtils::logfile);
	}

	// Type of case (the first entry on each line is the keyword describing the body case)
	std::string bodyCase;

	// Start reading in config file
	int bodyID = 0;
	while(file) {

		// Get type of body
		file >> bodyCase;

		// ** READ FROM FILE ** //
		if (bodyCase == "FROM_FILE") {

			// Read in the rest of the data for this case
			std::string boundaryType; file >> boundaryType;
			std::string fileName; file >> fileName;
			int lev; file >> lev;
			int reg; file >> reg;
			double startX; file >> startX;
			double startY; file >> startY;
			double centreZ; file >> centreZ;
			double length; file >> length;
			std::string direction; file >> direction;

			// Get body type
			eObjectType bodyType;
			if (boundaryType == "BBB")
				bodyType = eBBBCloud;
			else if (boundaryType == "BFL")
				bodyType = eBFLCloud;
			else if (boundaryType == "IBM")
				bodyType = eIBBCloud;

			// Get direction
			eCartesianDirection cartDirection;
			if (direction == "eXDirection")
				cartDirection = eXDirection;
			else if (direction == "eYDirection")
				cartDirection = eYDirection;
			else if (direction == "eZDirection")
				cartDirection = eZDirection;

			*GridUtils::logfile << "Initialising Body " << bodyID << " (" << boundaryType << ") from file..." << std::endl;

			// Read in data from point cloud file
			PCpts* _PCpts = NULL;
			_PCpts = new PCpts();
			this->io_readInCloud(_PCpts, bodyType, bodyID, fileName, lev, reg, startX, startY, centreZ, length, cartDirection);
			delete _PCpts;
			*GridUtils::logfile << "Finished creating Body " << bodyID << "..." << std::endl;

		}

		// Incremement body ID and set case to none
		bodyID++;
		bodyCase = "NONE";
	}
	file.close();


//
//#ifdef L_IBM_ON
//
//	*GridUtils::logfile << "Initialising IBM Objects..." << std::endl;
//
//	// Build a body
//	//		body_type == 1 is a rectangle/cuboid with rigid IBM,
//	//		body_type == 2 is a circle/sphere with rigid IBM,
//	//		body_type == 3 is a multi-body test case featuring both the above with rigid IBM
//	//		body_type == 4 is a single inextensible flexible filament with Jacowire IBM
//	//		body_type == 5 is an array of flexible filaments with Jacowire IBM
//	//		body_type == 6 is a plate in 2D with rigid IBM
//	//		body_type == 7 is the same as the previous case but with a Jacowire flexible flap added to the trailing edge
//	//		body_type == 8 is a plate in 3D with rigid IBM
//	//		body_type == 9 is the same as the previous case but with a rigid but moving filament array commanded by a single 2D Jacowire filament
//
//#if defined L_INSERT_RECTANGLE_CUBOID
//	this->ibm_buildBody(1);
//	*GridUtils::logfile << "Case: Rectangle/Cuboid using IBM" << std::endl;
//
//#elif defined L_INSERT_CIRCLE_SPHERE
//	this->ibm_buildBody(2);
//	*GridUtils::logfile << "Case: Circle/Sphere using IBM" << std::endl;
//
//#elif defined L_INSERT_BOTH
//	this->ibm_buildBody(3);
//	*GridUtils::logfile << "Case: Rectangle/Cuboid + Circle/Sphere using IBM" << std::endl;
//
//#elif defined L_INSERT_FILAMENT
//	this->ibm_buildBody(4);
//	*GridUtils::logfile << "Case: Single 2D filament using Jacowire IBM" << std::endl;
//
//#elif defined L_INSERT_FILARRAY
//	this->ibm_buildBody(5);
//	*GridUtils::logfile << "Case: Array of filaments using Jacowire IBM" << std::endl;
//
//#elif defined L_2D_RIGID_PLATE_IBM
//	this->ibm_buildBody(6);
//	*GridUtils::logfile << "Case: 2D rigid plate using IBM" << std::endl;
//
//#elif defined L_2D_PLATE_WITH_FLAP
//	this->ibm_buildBody(7);
//	*GridUtils::logfile << "Case: 2D rigid plate using IBM with flexible flap" << std::endl;
//
//#elif defined L_3D_RIGID_PLATE_IBM
//	this->ibm_buildBody(8);
//	*GridUtils::logfile << "Case: 3D rigid plate using IBM" << std::endl;
//
//#elif defined L_3D_PLATE_WITH_FLAP
//	this->ibm_buildBody(9);
//	*GridUtils::logfile << "Case: 3D rigid plate using IBM with flexible 2D flap" << std::endl;
//
//#endif
//
//
//
//
//#endif // End IBM_ON

}


// ************************************************************************** //
/// \brief	Read in point cloud data.
///
///			Input data must be in tab separated, 3-column format in the input 
///			directory.
///
/// \param	_PCpts	pointer to empty point cloud data container.
/// \param	objtype	type of object to be read in.
void ObjectManager::io_readInCloud(PCpts* _PCpts, eObjectType objtype, int bodyID, std::string fileName, int on_grid_lev, int on_grid_reg,
		double body_start_x, double body_start_y, double body_centre_z, double body_length, eCartesianDirection scale_direction) {

	// Temporary variables
	double tmp_x, tmp_y, tmp_z;
	int a = 0;

	// Case-specific variables
	GridObj* g = NULL;

	// Open input file
	std::ifstream file;
	file.open("./input/" + fileName, std::ios::in);

	// Handle failure to open
	if (!file.is_open()) {
		L_ERROR("Error opening cloud input file. Exiting.", GridUtils::logfile);
	}

	// Get grid pointer
	GridUtils::getGrid(_Grids, on_grid_lev, on_grid_reg, g);

	// Return if this process does not have this grid
	if (g == NULL) return;

	// Get rank for debugging
	int rank = GridUtils::safeGetRank();

	// Round bounding box dimensions to nearest voxel edge on this grid
	body_start_x = std::round(body_start_x / g->dh) * g->dh;
	body_start_y = std::round(body_start_y / g->dh) * g->dh;
	body_centre_z = std::round(body_centre_z / g->dh) * g->dh + g->dh / 2.0;	// Centre shifted to voxel centre
	body_length = std::round(body_length / g->dh) * g->dh;

	// Loop over lines in file
	while (!file.eof()) {

		// Read in one line of file at a time
		std::string line_in;	// String to store line in
		std::istringstream iss;	// Buffer stream to store characters

		// Get line up to new line separator and put in buffer
		std::getline(file, line_in, '\n');
		iss.str(line_in);	// Put line in the buffer
		iss.seekg(0);		// Reset buffer position to start of buffer

		// Add coordinates to data store
		iss >> tmp_x;
		iss >> tmp_y;
		iss >> tmp_z;

		_PCpts->x.push_back(tmp_x);
		_PCpts->y.push_back(tmp_y);

		// If running a 2D calculation, only read in x and y coordinates and force z coordinates to match the domain
#if (L_DIMS == 3)
		_PCpts->z.push_back(tmp_z);
#else
		_PCpts->z.push_back(0);
#endif

		// Insert the ID of the point within this point cloud (needed later for assigning marker IDs)
		_PCpts->id.push_back(_PCpts->id.size());

	}
	file.close();

	// Error if no data
	if (_PCpts->x.empty() || _PCpts->y.empty() || _PCpts->z.empty()) {
		L_ERROR("Failed to read object data from cloud input file.", GridUtils::logfile);
	}
	else {
		*GridUtils::logfile << "Successfully acquired object data from cloud input file." << std::endl;
	}


	// Rescale coordinates to fit into size required
#ifdef L_CLOUD_DEBUG
	*GridUtils::logfile << "Rescaling..." << std::endl;
#endif

	double scale_factor;
	// Scale slightly smaller (to voxel centres) to ensure symmetrical distribution of voxels
	if (scale_direction == eXDirection) {
		scale_factor = (body_length - g->dh) /
			std::fabs(*std::max_element(_PCpts->x.begin(), _PCpts->x.end()) - *std::min_element(_PCpts->x.begin(), _PCpts->x.end()));
	}
	else if (scale_direction == eYDirection) {
		scale_factor = (body_length - g->dh) /
			std::fabs(*std::max_element(_PCpts->y.begin(), _PCpts->y.end()) - *std::min_element(_PCpts->y.begin(), _PCpts->y.end()));
	}
	else if (scale_direction == eZDirection) {
		scale_factor = (body_length - g->dh) /
		std::fabs(*std::max_element(_PCpts->z.begin(), _PCpts->z.end()) - *std::min_element(_PCpts->z.begin(), _PCpts->z.end()));
	}

	// Shift to voxel centre not to edge
	double shift_x = (body_start_x + g->dh / 2.0) - scale_factor * *std::min_element(_PCpts->x.begin(), _PCpts->x.end());
	double shift_y = (body_start_y + g->dh / 2.0) - scale_factor * *std::min_element(_PCpts->y.begin(), _PCpts->y.end());

	// z-shift based on centre of object
	double scaled_body_length = scale_factor * (*std::max_element(_PCpts->z.begin(), _PCpts->z.end()) - *std::min_element(_PCpts->z.begin(), _PCpts->z.end()));
	double z_start_position = body_centre_z - (scaled_body_length / 2);
	double shift_z = z_start_position - scale_factor * *std::min_element(_PCpts->z.begin(), _PCpts->z.end());

	// Apply to each point to convert to global positions
	for (a = 0; a < static_cast<int>(_PCpts->x.size()); a++) {
		_PCpts->x[a] *= scale_factor; _PCpts->x[a] += shift_x;
		_PCpts->y[a] *= scale_factor; _PCpts->y[a] += shift_y;
		_PCpts->z[a] *= scale_factor; _PCpts->z[a] += shift_z;
	}

	// Write out the points after scaling and shifting
#ifdef L_CLOUD_DEBUG
	if (!_PCpts->x.empty()) {
		if (rank == 0) {
			std::ofstream fileout;
			fileout.open(GridUtils::path_str + "/CloudPtsPreFilter_Body" + std::to_string(bodyID) + "_Rank" + std::to_string(rank) + ".out", std::ios::out);
			for (size_t i = 0; i < _PCpts->x.size(); i++) {
				fileout << std::to_string(_PCpts->x[i]) + '\t' + std::to_string(_PCpts->y[i]) + '\t' + std::to_string(_PCpts->z[i]) + '\t' + std::to_string(_PCpts->id[i]);
				fileout << std::endl;
			}
			fileout.close();
		}
	}
#endif

	// Declare local indices
	std::vector<int> ijk;
	eLocationOnRank loc = eNone;

	// Exclude points which are not on this rank
#ifdef L_CLOUD_DEBUG
	*GridUtils::logfile << "Filtering..." << std::endl;
#endif
	a = 0;
	do {

		// If on this rank get its indices
		if (GridUtils::isOnThisRank(_PCpts->x[a], _PCpts->y[a], _PCpts->z[a], &loc, g))
		{
			// Increment counter
			a++;
		}
		// If not, erase
		else {
			_PCpts->x.erase(_PCpts->x.begin() + a);
			_PCpts->y.erase(_PCpts->y.begin() + a);
			_PCpts->z.erase(_PCpts->z.begin() + a);
			_PCpts->id.erase(_PCpts->id.begin() + a);

		}

	} while (a < static_cast<int>(_PCpts->x.size()));

	// Write out the points remaining in for debugging purposes
#ifdef L_CLOUD_DEBUG
	*GridUtils::logfile << "There are " << std::to_string(_PCpts->x.size()) << " points on this rank." << std::endl;
	*GridUtils::logfile << "Writing to file..." << std::endl;
	if (!_PCpts->x.empty()) {
		std::ofstream fileout;
		fileout.open(GridUtils::path_str + "/CloudPts_Body" + std::to_string(bodyID) + "_Rank" + std::to_string(rank) + ".out",std::ios::out);
		for (size_t i = 0; i < _PCpts->x.size(); i++) {
			fileout << std::to_string(_PCpts->x[i]) + '\t' + std::to_string(_PCpts->y[i]) + '\t' + std::to_string(_PCpts->z[i]) + '\t' + std::to_string(_PCpts->id[i]);
			fileout << std::endl;
		}
		fileout.close();
	}
#endif

	// If there are points left
	if (!_PCpts->x.empty())	{


		// Perform a different post-processing action depending on the type of body
		switch (objtype)
		{

		case eBBBCloud:

#ifdef L_CLOUD_DEBUG
			*GridUtils::logfile << "Labelling..." << std::endl;
#endif

			// Label the grid sites
			for (a = 0; a < static_cast<int>(_PCpts->x.size()); a++) {

				// Get indices if on this rank
				if (GridUtils::isOnThisRank(_PCpts->x[a], _PCpts->y[a], _PCpts->z[a], &loc, g, &ijk))
				{
					// Update Typing Matrix
					if (g->LatTyp(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim) == eFluid)
					{
						g->LatTyp(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim) = eSolid;
					}
				}
			}
			break;

		case eBFLCloud:

#ifdef L_CLOUD_DEBUG
			*GridUtils::logfile << "Building..." << std::endl;
#endif
			// Call BFL body builder
			bfl_buildBody(_PCpts, bodyID);
			break;

		case eIBBCloud:

#ifdef L_CLOUD_DEBUG
			*GridUtils::logfile << "Building..." << std::endl;
#endif
			// Call IBM body builder
			ibm_buildBody(_PCpts, g, bodyID);
			break;

		}
	}
}
// *****************************************************************************
/// \brief	Write out the forces on a solid object.
///
///			Writes out the forces on solid objects in the domain computed using
///			momentum exchange. Each rank writes its own file. Output is a CSV file.
/// \param	tval	time value at which write out is taking place.
void ObjectManager::io_writeForceOnObject(double tval) {
	
	int rank = GridUtils::safeGetRank();

	// Get grid on which object resides
	GridObj *g = NULL;
	GridUtils::getGrid(_Grids, L_OBJECT_ON_GRID_LEV, L_OBJECT_ON_GRID_REG, g);
	// If this grid exists on this process
	if (g != NULL)
	{
		// Create stream
		std::ofstream fout;
		fout.precision(L_OUTPUT_PRECISION);
		// Filename
		std::stringstream fileName;
		fileName << GridUtils::path_str + "/LiftDrag" << "Rnk" << rank << ".csv";
		// Open file
		fout.open(fileName.str().c_str(), std::ios::out | std::ios::app);
		// Write out the header (first time step only)
		if (static_cast<int>(tval) == 0) fout << "Time,Fx,Fy,Fz" << std::endl;

		fout << std::to_string(tval) << ","
			<< std::to_string(forceOnObjectX / pow(2, L_OBJECT_ON_GRID_LEV)) << ","
			<< std::to_string(forceOnObjectY / pow(2, L_OBJECT_ON_GRID_LEV)) << ","
#if (dims == 3)
			<< std::to_string(forceOnObjectZ / pow(2, L_OBJECT_ON_GRID_LEV))
#else
			<< std::to_string(0.0)
#endif
			<< std::endl;

		fout.close();
	}

}
// *****************************************************************************
