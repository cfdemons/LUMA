/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) The University of Manchester 2017
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * further distribution commericially or otherwise without written consent.
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
/// \brief	Read/write body information to restart file.
/// \param	IO_flag	flag indicating write (true) or read (false).
/// \param	level	level of the grid begin written/read
void ObjectManager::io_restart(eIOFlag IO_flag, int level) {

	int rank = GridUtils::safeGetRank();

	if (IO_flag == eWrite) {

		// Output stream
		std::ofstream file;

		// IB BODIES //

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


		// BFL BODIES //

		// TODO


	} else {

		// Input stream
		std::ifstream file;

		// IB BODIES //

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


		// BFL BODIES //

		// TODO

	}

}


// ************************************************************************** //
/// \brief	Wrapper for writing body position data to VTK file.
///
/// \param	tval	time value at which the write out is being performed.
void ObjectManager::io_vtkBodyWriter(int tval)
{

	// Get the rank
	int rank = GridUtils::safeGetRank();

    // Loop through each iBody
	for (IBBody& body : iBody)
	{
		// Call the writer
		body.writeVtkPosition(tval);

	}

    // Loop through each BFL Body
	for (BFLBody& body : pBody)
	{
		// Call the writer
		body.writeVtkPosition(tval);
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


	// Increment offset counter until first valid line is reached
	int fileOffset;
	std::string line;
	file.seekg(std::ios::beg);
	do {

		// Get the current position within the file
		fileOffset = static_cast<int>(file.tellg());

		// Get the whole line
		getline(file, line);

	} while (line[0] == '#');

	// Reset file position to the start of the last read line
	file.seekg(fileOffset, std::ios::beg);

	// Type of case (the first entry on each line is the keyword describing the body case)
	std::string bodyCase;

	// Start reading in config file
	int bodyID = 0;
	while(!file.eof()) {

		// Get type of body
		file >> bodyCase;

		// ** READ FROM FILE ** //
		if (bodyCase == "FROM_FILE")
		{

			// Read in the rest of the data for this case
			std::string boundaryType; file >> boundaryType;
			std::string fileName; file >> fileName;
			int lev; file >> lev;
			int reg; file >> reg;
			std::string xRefType; file >> xRefType;
			double xRef; file >> xRef;
			std::string yRefType; file >> yRefType;
			double yRef; file >> yRef;
			std::string zRefType; file >> zRefType;
			double zRef; file >> zRef;
			double length; file >> length;
			std::string direction; file >> direction;
			std::string flex_rigid; file >> flex_rigid;
			std::string BC; file >> BC;

			bool xRefCen = GeomPacked::interpretRef(xRefType);
			bool yRefCen = GeomPacked::interpretRef(yRefType);
			bool zRefCen = GeomPacked::interpretRef(zRefType);

			L_INFO("Initialising Body " + std::to_string(bodyID) + " (" + boundaryType + ") from file...", GridUtils::logfile);
			
			// Get body type
			eObjectType bodyType;
			if (boundaryType == "BBB")
				bodyType = eBBBCloud;
			else if (boundaryType == "BFL")
				bodyType = eBFLCloud;
			else if (boundaryType == "IBM")
				bodyType = eIBBCloud;

			// Get direction
			eCartesianDirection scaleDirection;
			if (direction == "X")
				scaleDirection = eXDirection;
			else if (direction == "Y")
				scaleDirection = eYDirection;
			else if (direction == "Z")
				scaleDirection = eZDirection;

			// Check if flexible (note: BFL is always rigid no matter what the input is)
			eMoveableType moveProperty;
			if (flex_rigid == "FLEXIBLE") {
				moveProperty = eFlexible;
				hasMovingBodies = true;
			}
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
				hasMovingBodies = true;
			}
			else if (flex_rigid == "RIGID")
				moveProperty = eRigid;
			else
				moveProperty = eRigid;

			// Get fixed BC
			bool clamped;
			if (BC == "CLAMPED")
				clamped = true;
			else if (BC == "SUPPORTED")
				clamped = false;
			else
				clamped = false;

			// Packed information into a geometry instance
			GeomPacked *geom = new GeomPacked(
				bodyType, bodyID, fileName, lev, reg, 
				xRefCen, xRef, yRefCen, yRef, zRefCen, zRef,
				length, scaleDirection, moveProperty, clamped
				);

			// Read in data from point cloud file
			PCpts* _PCpts = NULL;
			_PCpts = new PCpts();
			
			L_INFO("Reading in point cloud...", GridUtils::logfile); 
			this->io_readInCloud(_PCpts, geom);
			delete _PCpts;
			delete geom;
			*GridUtils::logfile << "Finished creating Body " << bodyID << "..." << std::endl;
		}

		// ** INSERT FILAMENT ** //
		else if (bodyCase == "FILAMENT")
		{

			// Read in the rest of the data for this case
			std::string boundaryType; file >> boundaryType;
			int lev; file >> lev;
			int reg; file >> reg;
			double startX; file >> startX;
			double startY; file >> startY;
			double startZ; file >> startZ;
			double length; file >> length;
			double angleVert; file >> angleVert;
			double angleHorz; file >> angleHorz;
			std::string flex_rigid; file >> flex_rigid;
			std::string BC; file >> BC;

			*GridUtils::logfile << "Initialising Body " << bodyID << " (" << boundaryType << ") as a filament..." << std::endl;

			// Sort data
			std::vector<double> start_position, angles;
			start_position.push_back(startX);
			start_position.push_back(startY);
			start_position.push_back(startZ);
			angles.push_back(angleVert);
			angles.push_back(angleHorz);

			// Check if flexible (note: BFL is always rigid no matter what the input is)
			eMoveableType moveProperty;
			if (flex_rigid == "FLEXIBLE") {
				moveProperty = eFlexible;
				hasMovingBodies = true;
			}
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
				hasMovingBodies = true;
			}
			else if (flex_rigid == "RIGID")
				moveProperty = eRigid;

			// Get fixed BC
			bool clamped;
			if (BC == "CLAMPED")
				clamped = true;
			else if (BC == "SUPPORTED")
				clamped = false;

			// Get grid pointer
			GridObj* g = NULL;
			GridUtils::getGrid(_Grids, lev, reg, g);

			// If rank has grid
			if (g != NULL) {

				// Build either BFL or IBM body constructor (note: most of the actual building takes place in the base constructor)
				if (boundaryType == "IBM") {
					iBody.emplace_back(g, bodyID, start_position, length, angles, moveProperty, clamped);

					// If no markers then get rid of the body
					if (iBody.back().markers.size() == 0)
						iBody.erase(iBody.end());
				}
				else if (boundaryType == "BFL") {
					pBody.emplace_back(g, bodyID, start_position, length, angles);

					// If no markers then get rid of the body
					if (pBody.back().markers.size() == 0)
						pBody.erase(pBody.end());
				}
			}
			*GridUtils::logfile << "Finished creating Body " << bodyID << "..." << std::endl;
		}

		// ** INSERT CIRCLE/SPHERE ** //
		else if (bodyCase == "CIRCLE_SPHERE")
		{

			// Read in the rest of the data for this case
			std::string boundaryType; file >> boundaryType;
			int lev; file >> lev;
			int reg; file >> reg;
			double centreX; file >> centreX;
			double centreY; file >> centreY;
			double centreZ; file >> centreZ;
			double radius; file >> radius;
			std::string flex_rigid; file >> flex_rigid;

			*GridUtils::logfile << "Initialising Body " << bodyID << " (" << boundaryType << ") as a circle/sphere..." << std::endl;

			// Sort data
			std::vector<double> centre_point, angles;
			centre_point.push_back(centreX);
			centre_point.push_back(centreY);
			centre_point.push_back(centreZ);

			// Check if flexible (note: BFL is always rigid no matter what the input is)
			eMoveableType moveProperty;
			if (flex_rigid == "FLEXIBLE")
				L_ERROR("Circle/sphere cannot be flexible. Exiting.", GridUtils::logfile);
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
				hasMovingBodies = true;
			}
			else if (flex_rigid == "RIGID")
				moveProperty = eRigid;

			// Get grid pointer
			GridObj* g = NULL;
			GridUtils::getGrid(_Grids, lev, reg, g);

			// If rank has grid
			if (g != NULL) {

				// Build either BFL or IBM body constructor (note: most of the actual building takes place in the base constructor)
				if (boundaryType == "IBM") {
					iBody.emplace_back(g, bodyID, centre_point, radius, moveProperty);

					// If no markers then get rid of the body
					if (iBody.back().markers.size() == 0)
						iBody.erase(iBody.end());
				}
				else if (boundaryType == "BFL") {
					pBody.emplace_back(g, bodyID, centre_point, radius);

					// If no markers then get rid of the body
					if (pBody.back().markers.size() == 0)
						pBody.erase(pBody.end());
				}
			}
			*GridUtils::logfile << "Finished creating Body " << bodyID << "..." << std::endl;
		}

		// ** INSERT SQUARE/CUBE ** //
		else if (bodyCase == "SQUARE_CUBE")
		{

			// Read in the rest of the data for this case
			std::string boundaryType; file >> boundaryType;
			int lev; file >> lev;
			int reg; file >> reg;
			double centreX; file >> centreX;
			double centreY; file >> centreY;
			double centreZ; file >> centreZ;
			double length; file >> length;
			double height; file >> height;
			double depth; file >> depth;
			double angleVert; file >> angleVert;
			double angleHorz; file >> angleHorz;
			std::string flex_rigid; file >> flex_rigid;

			*GridUtils::logfile << "Initialising Body " << bodyID << " (" << boundaryType << ") as a square/cube..." << std::endl;

			// Sort data
			std::vector<double> centre_point, dimensions, angles;
			centre_point.push_back(centreX);
			centre_point.push_back(centreY);
			centre_point.push_back(centreZ);
			dimensions.push_back(length);
			dimensions.push_back(height);
			dimensions.push_back(depth);
			angles.push_back(angleVert);
			angles.push_back(angleHorz);

			// Check if flexible (note: BFL is always rigid no matter what the input is)
			eMoveableType moveProperty;
			if (flex_rigid == "FLEXIBLE")
				L_ERROR("Circle/sphere cannot be flexible. Exiting.", GridUtils::logfile);
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
				hasMovingBodies = true;
			}
			else if (flex_rigid == "RIGID")
				moveProperty = eRigid;

			// Get grid pointer
			GridObj* g = NULL;
			GridUtils::getGrid(_Grids, lev, reg, g);

			// If rank has grid
			if (g != NULL) {

				// Build either BFL or IBM body constructor (note: most of the actual building takes place in the base constructor)
				if (boundaryType == "IBM") {
					iBody.emplace_back(g, bodyID, centre_point, dimensions, angles, moveProperty);

					// If no markers then get rid of the body
					if (iBody.back().markers.size() == 0)
						iBody.erase(iBody.end());
				}
				else if (boundaryType == "BFL") {
					pBody.emplace_back(g, bodyID, centre_point, dimensions, angles);

					// If no markers then get rid of the body
					if (pBody.back().markers.size() == 0)
						pBody.erase(pBody.end());
				}
			}
			*GridUtils::logfile << "Finished creating Body " << bodyID << "..." << std::endl;
		}

		// Increment body ID and set case to none
		bodyID++;
		bodyCase = "NONE";
	}
	file.close();
}


// ************************************************************************** //
/// \brief	Read in point cloud data.
///
///			Input data must be in tab separated, 3-column format in the input 
///			directory.
///
/// \param	_PCpts	reference to pointer to empty point cloud data container.
/// \param	geom	structure containing object data as parsed from the config file.
void ObjectManager::io_readInCloud(PCpts*& _PCpts, GeomPacked *geom)
{

	// Temporary variables
	double tmp_x, tmp_y, tmp_z;
	int a = 0;

	// Case-specific variables
	GridObj* g = NULL;

	// Open input file
	std::ifstream file;
	file.open("./input/" + geom->fileName, std::ios::in);

	// Handle failure to open
	if (!file.is_open()) {
		L_ERROR("Error opening cloud input file. Exiting.", GridUtils::logfile);
	}

	// Get grid pointer
	GridUtils::getGrid(_Grids, geom->onGridLev, geom->onGridReg, g);

	// Return if this process does not have this grid
	if (g == NULL) return;

	// Get rank for debugging
	int rank = GridUtils::safeGetRank();

	// Round reference values to voxel centres
	double bodyRefX = std::round(geom->bodyRefX / g->dh) * g->dh;
	double bodyRefY = std::round(geom->bodyRefY / g->dh) * g->dh;
	double bodyRefZ = std::round(geom->bodyRefZ / g->dh) * g->dh;
	double bodyLength = std::round(geom->bodyLength / g->dh) * g->dh;

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
		_PCpts->id.push_back(static_cast<int>(_PCpts->id.size()));

	}
	file.close();

	// Error if no data
	if (_PCpts->x.empty() || _PCpts->y.empty() || _PCpts->z.empty()) {
		L_ERROR("Failed to read object data from cloud input file.", GridUtils::logfile);
	}
	else {
		L_INFO("Successfully acquired object data from cloud input file.", GridUtils::logfile);
	}


	// Rescale coordinates to fit into size required
#ifdef L_CLOUD_DEBUG
	*GridUtils::logfile << "Rescaling..." << std::endl;
#endif

	double scale_factor;
	// Scale slightly smaller as distribution of voxels would be asymmetric if points sit on edge
	if (geom->scaleDirection == eXDirection)
	{
		scale_factor = (bodyLength - 2 * L_SMALL_NUMBER * g->dh) /
			std::fabs(*std::max_element(_PCpts->x.begin(), _PCpts->x.end()) - *std::min_element(_PCpts->x.begin(), _PCpts->x.end()));
	}
	else if (geom->scaleDirection == eYDirection)
	{
		scale_factor = (bodyLength - 2 * L_SMALL_NUMBER * g->dh) /
			std::fabs(*std::max_element(_PCpts->y.begin(), _PCpts->y.end()) - *std::min_element(_PCpts->y.begin(), _PCpts->y.end()));
	}
	else if (geom->scaleDirection == eZDirection)
	{
		scale_factor = (bodyLength - 2 * L_SMALL_NUMBER * g->dh) /
		std::fabs(*std::max_element(_PCpts->z.begin(), _PCpts->z.end()) - *std::min_element(_PCpts->z.begin(), _PCpts->z.end()));
	}

	// If reference is a centre, shift to centre of voxel
	double shiftX, shiftY, shiftZ, scaledDistance, startPos;
	if (geom->isRefXCentre)
	{
		bodyRefX += (g->dh / 2.0);
		scaledDistance = scale_factor * 
			std::fabs(*std::max_element(_PCpts->x.begin(), _PCpts->x.end()) - *std::min_element(_PCpts->x.begin(), _PCpts->x.end()));
		scaledDistance = std::round(scaledDistance / g->dh) * g->dh;	// Round to nearest voxel multiple
		startPos = bodyRefX - (scaledDistance / 2.0);
		shiftX = startPos - scale_factor * *std::min_element(_PCpts->x.begin(), _PCpts->x.end());
	}
	else
	{
		shiftX = (bodyRefX + L_SMALL_NUMBER * g->dh) - scale_factor * *std::min_element(_PCpts->x.begin(), _PCpts->x.end());
	}

	if (geom->isRefYCentre)
	{
		bodyRefY += (g->dh / 2.0);
		scaledDistance = scale_factor *
			std::fabs(*std::max_element(_PCpts->y.begin(), _PCpts->y.end()) - *std::min_element(_PCpts->y.begin(), _PCpts->y.end()));
		scaledDistance = std::round(scaledDistance / g->dh) * g->dh;
		startPos = bodyRefY - (scaledDistance / 2.0);
		shiftY = startPos - scale_factor * *std::min_element(_PCpts->y.begin(), _PCpts->y.end());
	}
	else
	{
		shiftY = (bodyRefY + L_SMALL_NUMBER * g->dh) - scale_factor * *std::min_element(_PCpts->y.begin(), _PCpts->y.end());
	}

	if (geom->isRefZCentre)
	{
		bodyRefZ += (g->dh / 2.0);
		scaledDistance = scale_factor *
			std::fabs(*std::max_element(_PCpts->z.begin(), _PCpts->z.end()) - *std::min_element(_PCpts->z.begin(), _PCpts->z.end()));
		scaledDistance = std::round(scaledDistance / g->dh) * g->dh;
		startPos = bodyRefZ - (scaledDistance / 2.0);
		shiftZ = startPos - scale_factor * *std::min_element(_PCpts->z.begin(), _PCpts->z.end());
	}
	else
	{
		shiftZ = (bodyRefZ + L_SMALL_NUMBER * g->dh) - scale_factor * *std::min_element(_PCpts->z.begin(), _PCpts->z.end());
	}

	// Declare local indices
	std::vector<int> ijk;
	eLocationOnRank loc = eNone;

	// Filter: erase is O(n^2) so propose creating a copy instead
	PCpts *_filtered = new PCpts();

	// Apply shift and scale to each point to convert to global positions
	for (a = 0; a < static_cast<int>(_PCpts->x.size()); a++)
	{
		_PCpts->x[a] *= scale_factor; _PCpts->x[a] += shiftX;
		_PCpts->y[a] *= scale_factor; _PCpts->y[a] += shiftY;
#if (L_DIMS == 3)
		_PCpts->z[a] *= scale_factor; _PCpts->z[a] += shiftZ;
#endif

		// Apply a rnak filter at the same time
		if (GridUtils::isOnThisRank(_PCpts->x[a], _PCpts->y[a], _PCpts->z[a], &loc, g))
		{
			_filtered->x.push_back(_PCpts->x[a]);
			_filtered->y.push_back(_PCpts->y[a]);
			_filtered->z.push_back(_PCpts->z[a]);
			_filtered->id.push_back(_PCpts->id[a]);
		}

	}

	// Free old array and assign new array to pointer which will be passed back out
	delete _PCpts;
	_PCpts = _filtered;

	// Write out the points after scaling, shifting and filtering
#ifdef L_CLOUD_DEBUG
	if (!_PCpts->x.empty()) {
		if (rank == 0) {
			std::ofstream fileout;
			fileout.open(GridUtils::path_str + "/CloudPts_Body" + std::to_string(geom->bodyID) + "_Rank" + std::to_string(rank) + ".out", std::ios::out);
			for (size_t i = 0; i < _PCpts->x.size(); i++) {
				fileout << std::to_string(_PCpts->x[i]) + '\t' + std::to_string(_PCpts->y[i]) + '\t' + std::to_string(_PCpts->z[i]) + '\t' + std::to_string(_PCpts->id[i]);
				fileout << std::endl;
			}
			fileout.close();
		}
	}
#endif	

	// If there are points left
	if (!_PCpts->x.empty())
	{

#ifdef L_CLOUD_DEBUG
		*GridUtils::logfile << "Building..." << std::endl;
#endif

		// Perform a different post-processing action depending on the type of body
		switch (geom->objtype)
		{

		case eBBBCloud:

			// Call labeller for BBB
			addBouncebackObject(g, geom, _PCpts);
			break;

		case eBFLCloud:

			// Call constructor to build BFL body
			pBody.emplace_back(g, geom->bodyID, _PCpts);
			break;

		case eIBBCloud:

			// Call constructor to build IBM body
			iBody.emplace_back(g, geom->bodyID, _PCpts, geom->moveProperty, geom->isClamped);
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
void ObjectManager::io_writeForcesOnObjects(double tval) {
	
	// Declarations
	int rank = GridUtils::safeGetRank();
	std::ofstream fout;
	fout.precision(L_OUTPUT_PRECISION);
	std::stringstream fileName;

	// BB OBJECTS //

	// Get grid on which object resides
	GridObj *g = NULL;
	GridUtils::getGrid(_Grids, bbbOnGridLevel, bbbOnGridReg, g);
	// If this grid exists on this process
	if (g != NULL)
	{
		// Filename
		fileName << GridUtils::path_str + "/LiftDragBBB_Rnk" << rank << ".csv";

		// Open file
		fout.open(fileName.str().c_str(), std::ios::out | std::ios::app);

		// Write out the header (first time step only)
		if (static_cast<int>(tval) == L_OUT_EVERY_FORCES) fout << "Time,Fx,Fy,Fz" << std::endl;

		// Scaled with respect to refinement ratio
		fout << std::to_string(tval) << ","
			<< std::to_string(bbbForceOnObjectX * g->refinement_ratio) << ","
			<< std::to_string(bbbForceOnObjectY * g->refinement_ratio) << ","
#if (dims == 3)
			<< std::to_string(bbbForceOnObjectZ * g->refinement_ratio)
#else
			<< std::to_string(0.0)
#endif
			<< std::endl;

		fout.close();
	}

	// BFL OBJECTS //
	ObjectManager *objman = ObjectManager::getInstance();
	
	// Search for bodies on this rank
	for (BFLBody& body : objman->pBody)
	{
		// Filename
		fileName << GridUtils::path_str + "/LiftDragBFL_Rnk" << rank << ".csv";

		// Open file
		fout.open(fileName.str().c_str(), std::ios::out | std::ios::app);

		// Write out the header (first time step only)
		if (static_cast<int>(tval) == L_OUT_EVERY_FORCES) fout << "Time,Fx,Fy,Fz" << std::endl;

		// Summation required
		double bodyForceX = 0.0;
		double bodyForceY = 0.0;
		double bodyForceZ = 0.0;

		for (BFLMarker& marker : body.markers)
		{
			// Sum marker contributions
			bodyForceX += marker.forceX;
			bodyForceY += marker.forceY;
			bodyForceZ += marker.forceZ;

		}

		// Scaled with respect to refinement ratio
		fout << std::to_string(tval) << ","
			<< std::to_string(bodyForceX * body._Owner->refinement_ratio) << ","
			<< std::to_string(bodyForceY * body._Owner->refinement_ratio) << ","
#if (dims == 3)
			<< std::to_string(bodyForceZ * body._Owner->refinement_ratio)
#else
			<< std::to_string(0.0)
#endif
			<< std::endl;

		fout.close();
	}

}
// *****************************************************************************
