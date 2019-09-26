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

#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"
#include "../inc/PCpts.h"
#include "../inc/GridObj.h"


// *****************************************************************************
///	\brief	Write out position of immersed boundary bodies
///
///	\param	timestep		current time step
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


// *****************************************************************************
///	\brief	Write out forces on the markers of immersed boundary bodies
void ObjectManager::io_writeLiftDrag() {

	// Force conversion
	double forceScaling, volWidth, volDepth;

	// Write out lift and drag
	int rank = GridUtils::safeGetRank();

	// Loop through all bodies
	for (size_t ib = 0; ib < iBody.size(); ib++) {

		// If owning and flexible then write or if body is rigid and this rank has some markers
		if ((iBody[ib].isFlexible && iBody[ib].owningRank == rank) || (!iBody[ib].isFlexible && iBody[ib].validMarkers.size() > 0)) {

			// Sort out which markers to write out
			std::vector<int> validMarkers;

			// If flexible then write out all markers, if rigid only markers on this rank
			if (iBody[ib].isFlexible && iBody[ib].owningRank == rank)
				validMarkers = GridUtils::onespace(0, static_cast<int>(iBody[ib].markers.size())-1);
			else if (!iBody[ib].isFlexible && iBody[ib].validMarkers.size() > 0)
				validMarkers = iBody[ib].validMarkers;

			// Open file for given time step
			std::ofstream jout;
			jout.open(GridUtils::path_str + "/Body_" + std::to_string(iBody[ib].id) + "_LD_rank" + std::to_string(rank) + ".out", std::ios::app);
			jout.precision(L_OUTPUT_PRECISION);

			// If first time then write out header
			if (_Grids->t == 0)
				jout << "Timestep\tTime (s)\tDrag (N)\tLift (N)" << std::endl;

			// Convert force per volume to force
#if (L_DIMS == 2)
			forceScaling = iBody[ib]._Owner->dm * iBody[ib]._Owner->dh / SQ(iBody[ib]._Owner->dt) * 1.0 / iBody[ib]._Owner->dh;
#elif (L_DIMS == 3)
			forceScaling = iBody[ib]._Owner->dm * iBody[ib]._Owner->dh / SQ(iBody[ib]._Owner->dt);
#endif

			// Write out timestep data
			jout << _Grids->t << "\t" << _Grids->t * _Grids->dt;

			// Compute lift and drag
			for (auto m : validMarkers) {

				// Get volume scaling
				volWidth = iBody[ib].markers[m].epsilon;
#if (L_DIMS == 2)
				volDepth = 1.0;
#elif (L_DIMS == 3)
				volDepth = iBody[ib].markers[m].ds;
#endif

				// Write out force on markers
				for (int dir = 0; dir < L_DIMS; dir++)
					jout << "\t" << iBody[ib].markers[m].force_xyz[dir] * volWidth * volDepth * iBody[ib].markers[m].ds * forceScaling;
			}

			// Next line and close file
			jout << std::endl;
			jout.close();
		}
	}
}


// *****************************************************************************
///	\brief	Read/write body information to restart file
///
///	\param	IO_flag		flag indicating write (true) or read (false)
///	\param	level		current grid level
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
					file << iBody[b].markers[m].position0[0] << "\t"
						<< iBody[b].markers[m].position0[1] << "\t"
						<< iBody[b].markers[m].position0[2] << "\t";
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
					iss		>> iBody[b].markers[m].position0[0]
							>> iBody[b].markers[m].position0[1]
							>> iBody[b].markers[m].position0[2];
				}
			}

		}

		// Close file
		file.close();


		// BFL BODIES //

		// TODO

	}

}


// *****************************************************************************
///	\brief	Wrapper for writing body position data to VTK file
///
///	\param	tval		time value at which the write out is being performed
void ObjectManager::io_vtkBodyWriter(int tval)
{

	// Get the rank
	int rank = GridUtils::safeGetRank();

    // Loop through each iBody
	for (IBBody& body : iBody) {

		// Call the writer
		if (body.owningRank == rank)
			body.writeVtkPosition(tval);
	}

    // Loop through each BFL Body
	for (BFLBody& body : pBody) {

		// Call the writer
		if (body.owningRank == rank)
			body.writeVtkPosition(tval);
	}
}


// *****************************************************************************
///	\brief	Read in geometry config file
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
	int iBodyID = 0, pBodyID = 0;
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

			L_INFO("Initialising Body " + std::to_string(iBodyID + pBodyID) + " (" + boundaryType + ") from file...", GridUtils::logfile);

			// Get body type
			eObjectType bodyType;
			if (boundaryType == "BBB")
				bodyType = eBBBCloud;
			else if (boundaryType == "BFL")
				bodyType = eBFLCloud;
			else if (boundaryType == "IBM") {
				bodyType = eIBBCloud;
				hasIBMBodies[lev] = true;
			}

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
				hasFlexibleBodies[lev] = true;
			}
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
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
				bodyType, iBodyID + pBodyID, fileName, lev, reg,
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
			*GridUtils::logfile << "Finished creating Body " << iBodyID + pBodyID << "..." << std::endl;

			// Increment counter
			if (bodyType == eBFLCloud)
				pBodyID++;
			else if (bodyType == eIBBCloud)
				iBodyID++;
		}

		// ** INSERT FILAMENT ** //
		else if (bodyCase == "FILAMENT_ARRAY")
		{

			// Read in the rest of the data for this case
			std::string boundaryType; file >> boundaryType;
			int lev; file >> lev;
			int reg; file >> reg;
			int nFil; file >> nFil;
			double startX; file >> startX;
			double startY; file >> startY;
			double startZ; file >> startZ;
			double spaceX; file >> spaceX;
			double spaceY; file >> spaceY;
			double spaceZ; file >> spaceZ;
			double length; file >> length;
			double height; file >> height;
			double depth; file >> depth;
			double angleVert; file >> angleVert;
			double angleHorz; file >> angleHorz;
			std::string flex_rigid; file >> flex_rigid;
			std::string nElementsString; file >> nElementsString;
			std::string BC; file >> BC;
			double density; file >> density;
			double YoungMod; file >> YoungMod;

			// Need to shift the body if using walls
			double shiftX = 0.0, shiftY = 0.0, shiftZ = 0.0;
			if (L_WALL_LEFT == eSolid)
				shiftX = L_WALL_THICKNESS_LEFT;
			if (L_WALL_BOTTOM == eSolid)
				shiftY = L_WALL_THICKNESS_BOTTOM;
			if (L_WALL_FRONT == eSolid)
				shiftZ = L_WALL_THICKNESS_FRONT;

			// Sort data
			std::vector<double> start_position, spacing, angles;
			start_position.push_back(startX + shiftX);
			start_position.push_back(startY + shiftY);
			start_position.push_back(startZ + shiftZ);
			spacing.push_back(spaceX);
			spacing.push_back(spaceY);
			spacing.push_back(spaceZ);
			angles.push_back(angleVert);
			angles.push_back(angleHorz);

			// Check if flexible (note: BFL is always rigid no matter what the input is)
			eMoveableType moveProperty;
			if (flex_rigid == "FLEXIBLE") {
				moveProperty = eFlexible;
				hasFlexibleBodies[lev] = true;
			}
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
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
			int nElements;
			if (g != NULL) {

				// Get number of elements to use in FEM
				if (nElementsString == "CONFORMING") {
					nElements = static_cast<int>(std::floor(length / g->dh));
				}
				else {
					nElements = static_cast<int>(std::stod(nElementsString));
				}
			}

			// Loop through and build bodies
			std::vector<double> position(3, 0.0);
			for (int i = 0; i < nFil; i++) {

				*GridUtils::logfile << "Initialising Body " << iBodyID + pBodyID << " (" << boundaryType << ") as a filament..." << std::endl;

				// Only build if this rank has this grid
				if (g != NULL) {

					// Loop through dimensions
					for (int d = 0; d < L_DIMS; d++)
						position[d] = start_position[d] + i * spacing[d];

					// Build either BFL or IBM body constructor (note: most of the actual building takes place in the base constructor)
					if (boundaryType == "IBM") {
						hasIBMBodies[lev] = true;
						iBody.emplace_back(g, iBodyID + pBodyID, position, length, height, depth, angles, moveProperty, nElements, clamped, density, YoungMod);
					}
					else if (boundaryType == "BFL") {
						pBody.emplace_back(g, iBodyID + pBodyID, position, length, angles);
					}
				}

				*GridUtils::logfile << "Finished creating Body " << iBodyID + pBodyID << "..." << std::endl;

				// Increment counter
				if (boundaryType == "IBM")
					iBodyID++;
				else if (boundaryType == "BFL")
					pBodyID++;
			}
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

			*GridUtils::logfile << "Initialising Body " << iBodyID + pBodyID << " (" << boundaryType << ") as a circle/sphere..." << std::endl;

			// Need to shift the body if using walls
			double shiftX = 0.0, shiftY = 0.0, shiftZ = 0.0;
			if (L_WALL_LEFT == eSolid)
				shiftX = L_WALL_THICKNESS_LEFT;
			if (L_WALL_BOTTOM == eSolid)
				shiftY = L_WALL_THICKNESS_BOTTOM;
			if (L_WALL_FRONT == eSolid)
				shiftZ = L_WALL_THICKNESS_FRONT;
			
			// Sort data
			std::vector<double> centre_point, angles;
			centre_point.push_back(centreX + shiftX);
			centre_point.push_back(centreY + shiftY);
			centre_point.push_back(centreZ + shiftZ);

			// Check if flexible (note: BFL is always rigid no matter what the input is)
			eMoveableType moveProperty;
			if (flex_rigid == "FLEXIBLE")
				L_ERROR("Circle/sphere cannot be flexible. Exiting.", GridUtils::logfile);
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
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
					hasIBMBodies[lev] = true;
					iBody.emplace_back(g, iBodyID + pBodyID, centre_point, radius, moveProperty);
				}
				else if (boundaryType == "BFL") {
					pBody.emplace_back(g, iBodyID + pBodyID, centre_point, radius);
				}
			}
			*GridUtils::logfile << "Finished creating Body " << iBodyID + pBodyID << "..." << std::endl;

			// Increment counter
			if (boundaryType == "IBM")
				iBodyID++;
			else if (boundaryType == "BFL")
				pBodyID++;
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

			*GridUtils::logfile << "Initialising Body " << iBodyID + pBodyID << " (" << boundaryType << ") as a square/cube..." << std::endl;

			// Need to shift the body if using walls
			double shiftX = 0.0, shiftY = 0.0, shiftZ = 0.0;
			if (L_WALL_LEFT == eSolid)
				shiftX = L_WALL_THICKNESS_LEFT;
			if (L_WALL_BOTTOM == eSolid)
				shiftY = L_WALL_THICKNESS_BOTTOM;
			if (L_WALL_FRONT == eSolid)
				shiftZ = L_WALL_THICKNESS_FRONT;

			// Sort data
			std::vector<double> centre_point, dimensions, angles;
			centre_point.push_back(centreX + shiftX);
			centre_point.push_back(centreY + shiftY);
			centre_point.push_back(centreZ + shiftZ);
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
					hasIBMBodies[lev] = true;
					iBody.emplace_back(g, iBodyID + pBodyID, centre_point, dimensions, angles, moveProperty);
				}
				else if (boundaryType == "BFL") {
					pBody.emplace_back(g, iBodyID + pBodyID, centre_point, dimensions, angles);
				}
			}
			*GridUtils::logfile << "Finished creating Body " << iBodyID + pBodyID << "..." << std::endl;

			// Increment counter
			if (boundaryType == "IBM")
				iBodyID++;
			else if (boundaryType == "BFL")
				pBodyID++;
		}


		// ** INSERT PLATE ** //
		else if (bodyCase == "PLATE")
		{

			// Can't build plate in 2D
#if (L_DIMS == 2)
			L_ERROR("Can't build plate in 2D. Exiting.", GridUtils::logfile);
#endif

			// Read in the rest of the data for this case
			std::string boundaryType; file >> boundaryType;
			int lev; file >> lev;
			int reg; file >> reg;
			double centreX; file >> centreX;
			double centreY; file >> centreY;
			double centreZ; file >> centreZ;
			double length; file >> length;
			double width; file >> width;
			double angleX; file >> angleX;
			double angleY; file >> angleY;
			double angleZ; file >> angleZ;
			std::string flex_rigid; file >> flex_rigid;

			*GridUtils::logfile << "Initialising Body " << iBodyID + pBodyID << " (" << boundaryType << ") as a plate..." << std::endl;

			// Need to shift the body if using walls
			double shiftX = 0.0, shiftY = 0.0, shiftZ = 0.0;
			if (L_WALL_LEFT == eSolid)
				shiftX = L_WALL_THICKNESS_LEFT;
			if (L_WALL_BOTTOM == eSolid)
				shiftY = L_WALL_THICKNESS_BOTTOM;
			if (L_WALL_FRONT == eSolid)
				shiftZ = L_WALL_THICKNESS_FRONT;

			// Sort data
			std::vector<double> centre_point, angles;
			centre_point.push_back(centreX + shiftX);
			centre_point.push_back(centreY + shiftY);
			centre_point.push_back(centreZ + shiftZ);
			angles.push_back(angleX);
			angles.push_back(angleY);
			angles.push_back(angleZ);

			// Check if flexible (note: BFL is always rigid no matter what the input is)
			eMoveableType moveProperty;
			if (flex_rigid == "FLEXIBLE")
				L_ERROR("Plate cannot be flexible. Exiting.", GridUtils::logfile);
			else if (flex_rigid == "MOVABLE") {
				moveProperty = eMovable;
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
					hasIBMBodies[lev] = true;
					iBody.emplace_back(g, iBodyID + pBodyID, centre_point, length, width, angles, moveProperty);
				}
				else if (boundaryType == "BFL") {
					pBody.emplace_back(g, iBodyID + pBodyID, centre_point, length, width, angles);
				}
			}
			*GridUtils::logfile << "Finished creating Body " << iBodyID + pBodyID << "..." << std::endl;

			// Increment counter
			if (boundaryType == "IBM")
				iBodyID++;
			else if (boundaryType == "BFL")
				pBodyID++;
		}

		// Increment body ID and set case to none
		bodyCase = "NONE";
	}
	file.close();

	L_INFO("Finalising bodies...", GridUtils::logfile);

	// Do some more IBM setup required after reading all bodies
	ibm_finaliseReadIn(iBodyID);

	L_INFO("Finalising Complete.", GridUtils::logfile);
}



// *****************************************************************************
/// \brief	Read in point cloud data
///
///			Input data must be in tab separated, 3-column format in the input
///			directory.
///
///	\param	_PCpts		reference to pointer to empty point cloud data container
///	\param	geom		structure containing object data as parsed from the config file
void ObjectManager::io_readInCloud(PCpts*& _PCpts, GeomPacked *geom)
{

	// Temporary variables
	double tmp_x, tmp_y, tmp_z, dCell;
	int a = 0;

	// Case-specific variables
	GridObj* g = NULL;

	// Open input file
	std::ifstream file;
	file.open("./input/" + geom->fileName, std::ios::in);

	// Handle failure to open
	if (!file.is_open())
		L_ERROR("Error opening cloud input file: " + geom->fileName + ". Exiting.", GridUtils::logfile);

	// If the level is set to -1 then object can span levels
	if (geom->onGridLev < 0)
	{
		// For scaling use the finest grid scale
		dCell = _Grids[0].dh / pow(2, L_NUM_LEVELS);

		// For range checking use the coarsest grid
		g = _Grids;
	}
	else
	{
		// Get required grid pointer
		GridUtils::getGrid(_Grids, geom->onGridLev, geom->onGridReg, g);

		// Return if this process does not have this grid
		if (g == NULL) return;

		// Set scaling
		dCell = g->dh;
	}

	// Round reference values to complete number of voxels as measured from origin
	double bodyRefX = std::round(geom->bodyRefX / dCell) * dCell;
	double bodyRefY = std::round(geom->bodyRefY / dCell) * dCell;
	double bodyRefZ = std::round(geom->bodyRefZ / dCell) * dCell;
	double bodyLength = std::round(geom->bodyLength / dCell) * dCell;

	// Write the scaled reference data
#ifdef L_CLOUD_DEBUG
	std::string msg("Scaled reference values are:");
	msg += " X = " + std::to_string(bodyRefX);
	msg += " Y = " + std::to_string(bodyRefY);
	msg += " Z = " + std::to_string(bodyRefZ);
	msg += " Length = " + std::to_string(bodyLength);
	L_DEBUG(msg, GridUtils::logfile);
#endif

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
	if (_PCpts->x.empty() || _PCpts->y.empty() || _PCpts->z.empty())
		L_ERROR("Failed to read object data from cloud input file.", GridUtils::logfile);
	else
		L_INFO("Successfully acquired object data from cloud input file.", GridUtils::logfile);


	// Rescale coordinates to fit into size required
#ifdef L_CLOUD_DEBUG
	L_DEBUG("Rescaling...", GridUtils::logfile);
#endif

	double scale_factor;
	// Scale slightly smaller as distribution of voxels would be asymmetric if points sit on edge
	if (geom->scaleDirection == eXDirection)
	{
		scale_factor = (bodyLength - 2 * L_SMALL_NUMBER * dCell) /
			std::fabs(*std::max_element(_PCpts->x.begin(), _PCpts->x.end()) - *std::min_element(_PCpts->x.begin(), _PCpts->x.end()));
	}
	else if (geom->scaleDirection == eYDirection)
	{
		scale_factor = (bodyLength - 2 * L_SMALL_NUMBER * dCell) /
			std::fabs(*std::max_element(_PCpts->y.begin(), _PCpts->y.end()) - *std::min_element(_PCpts->y.begin(), _PCpts->y.end()));
	}
	else if (geom->scaleDirection == eZDirection)
	{
		scale_factor = (bodyLength - 2 * L_SMALL_NUMBER * dCell) /
		std::fabs(*std::max_element(_PCpts->z.begin(), _PCpts->z.end()) - *std::min_element(_PCpts->z.begin(), _PCpts->z.end()));
	}

	// If reference is a centre, shift to centre of voxel
	double shiftX, shiftY, shiftZ, scaledDistance, startPos;
	if (geom->isRefXCentre)
	{
		bodyRefX += (dCell / 2.0);
		scaledDistance = scale_factor * 
			std::fabs(*std::max_element(_PCpts->x.begin(), _PCpts->x.end()) - *std::min_element(_PCpts->x.begin(), _PCpts->x.end()));
		scaledDistance = std::round(scaledDistance / dCell) * dCell;	// Round to nearest voxel multiple
		startPos = bodyRefX - (scaledDistance / 2.0);
		shiftX = startPos - scale_factor * *std::min_element(_PCpts->x.begin(), _PCpts->x.end());
	}
	else
	{
		shiftX = (bodyRefX + L_SMALL_NUMBER * dCell) - scale_factor * *std::min_element(_PCpts->x.begin(), _PCpts->x.end());
	}

	if (geom->isRefYCentre)
	{
		bodyRefY += (dCell / 2.0);
		scaledDistance = scale_factor *
			std::fabs(*std::max_element(_PCpts->y.begin(), _PCpts->y.end()) - *std::min_element(_PCpts->y.begin(), _PCpts->y.end()));
		scaledDistance = std::round(scaledDistance / dCell) * dCell;
		startPos = bodyRefY - (scaledDistance / 2.0);
		shiftY = startPos - scale_factor * *std::min_element(_PCpts->y.begin(), _PCpts->y.end());
	}
	else
	{
		shiftY = (bodyRefY + L_SMALL_NUMBER * dCell) - scale_factor * *std::min_element(_PCpts->y.begin(), _PCpts->y.end());
	}

	if (geom->isRefZCentre)
	{
		bodyRefZ += (dCell / 2.0);
		scaledDistance = scale_factor *
			std::fabs(*std::max_element(_PCpts->z.begin(), _PCpts->z.end()) - *std::min_element(_PCpts->z.begin(), _PCpts->z.end()));
		scaledDistance = std::round(scaledDistance / dCell) * dCell;
		startPos = bodyRefZ - (scaledDistance / 2.0);
		shiftZ = startPos - scale_factor * *std::min_element(_PCpts->z.begin(), _PCpts->z.end());
	}
	else
	{
		shiftZ = (bodyRefZ + L_SMALL_NUMBER * dCell) - scale_factor * *std::min_element(_PCpts->z.begin(), _PCpts->z.end());
	}

	// Declare local indices
	std::vector<int> ijk;
	eLocationOnRank loc = eNone;

	// Filter: erase is O(n^2) so create a copy instead
	PCpts *_filtered = new PCpts();

	// Apply shift and scale to each point to convert to global positions
	for (a = 0; a < static_cast<int>(_PCpts->x.size()); a++)
	{
		_PCpts->x[a] *= scale_factor; _PCpts->x[a] += shiftX;
		_PCpts->y[a] *= scale_factor; _PCpts->y[a] += shiftY;
#if (L_DIMS == 3)
		_PCpts->z[a] *= scale_factor; _PCpts->z[a] += shiftZ;
#endif

		// Apply a rank filter at the same time
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
	if (!_PCpts->x.empty())
	{
		std::ofstream fileout;
		fileout.open(GridUtils::path_str + "/CloudPts_Body" + std::to_string(geom->bodyID) + "_Rank" + std::to_string(GridUtils::safeGetRank()) + ".out", std::ios::out);
		for (size_t i = 0; i < _PCpts->x.size(); i++)
		{
			fileout << std::to_string(_PCpts->x[i]) + '\t' + std::to_string(_PCpts->y[i]) + '\t' + std::to_string(_PCpts->z[i]) + '\t' + std::to_string(_PCpts->id[i]);
			fileout << std::endl;
		}
		fileout.close();
	}
#endif	

	// If there are points left
	if (!_PCpts->x.empty() && geom->objtype != eIBBCloud)
	{

		L_INFO("Building body on this rank...", GridUtils::logfile);

		// Perform a different post-processing action depending on the type of body
		switch (geom->objtype)
		{

		case eBBBCloud:

			// Call labeller for BBB
			if (geom->onGridLev < 0)
				addBouncebackObject(geom, _PCpts);		// Can cross over grid levels
			else
				addBouncebackObject(g, geom, _PCpts);
			break;

		case eBFLCloud:

			// Call constructor to build BFL body
			pBody.emplace_back(g, geom->bodyID, _PCpts);
			break;
		}
	}
	else if (geom->objtype == eIBBCloud) {

		// Call constructor to build IBM body
		iBody.emplace_back(g, geom->bodyID, _PCpts, geom->moveProperty);
	}
}


// *****************************************************************************
/// \brief	Write out the forces on a solid object
///
///			Writes out the forces on solid objects in the domain computed using
///			momentum exchange. Each rank writes its own file. Output is a CSV file.
///
///	\param	tval		time value at which write out is taking place
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
		if (static_cast<int>(tval) == L_EXTRA_OUT_FREQ) fout << "Time,Fx,Fy,Fz" << std::endl;

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
		if (static_cast<int>(tval) == L_EXTRA_OUT_FREQ) fout << "Time,Fx,Fy,Fz" << std::endl;

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
/// \brief	Write out FEM body to VTK
///
///	\param	tval		time value at which write out is taking place
void ObjectManager::io_vtkFEMWriter(int tval)
{

	// If there are no flexible bodies on this rank then exit
	if (idxFEM.size() == 0)
		return;

	// Get the rank
	int rank = GridUtils::safeGetRank();

	// Create string and file streams
	std::stringstream fileName;
	fileName << GridUtils::path_str + "/vtk_out.FEM" << std::to_string(rank) << "." << (int)tval << ".vtk";
	std::ofstream fout;
	fout.open(fileName.str().c_str());

	// Number of nodes and lines
	int nNodes = 0;
	std::vector<int> nLinesBody;

    // Loop through all bodies
	for (size_t ib = 0; ib < idxFEM.size(); ib++) {

		// Increment number of nodes and lines
		nNodes += static_cast<int>(iBody[idxFEM[ib]].fBody->nodes.size());
		nLinesBody.push_back(static_cast<int>(iBody[idxFEM[ib]].fBody->nodes.size()) - 1);
	}

	// Add header information
	fout << "# vtk DataFile Version 3.0f\n";
	fout << "Output for rank " << std::to_string(rank) << " at time t = " << (int)tval << "\n";
	fout << "ASCII\n";
	fout << "DATASET POLYDATA\n";

	// Write out the positions of each Lagrange marker
	fout << "POINTS " << nNodes << " float\n";

    // Loop through all bodies and there nodes
	for (size_t ib = 0; ib < idxFEM.size(); ib++) {
		for (size_t m = 0; m < iBody[idxFEM[ib]].fBody->nodes.size(); m++) {

			// Write out positions
			fout << iBody[idxFEM[ib]].fBody->nodes[m].position[eXDirection] << " "
				 << iBody[idxFEM[ib]].fBody->nodes[m].position[eYDirection] << " "
				 << iBody[idxFEM[ib]].fBody->nodes[m].position[eZDirection] << std::endl;
		}
	}

	// Write out the connectivity of each Lagrange marker
	fout << "LINES " << std::accumulate(nLinesBody.begin(), nLinesBody.end(), 0) << " " << 3 * std::accumulate(nLinesBody.begin(), nLinesBody.end(), 0) << std::endl;

	// Loop through number of lines
	int count = 0;
	for (size_t i = 0; i < nLinesBody.size(); i++) {
		for (int j = 0; j < nLinesBody[i]; j++) {
			fout << 2 << " " << count << " " << count + 1 << std::endl;
			count++;
		}
		count++;
	}

	// Close file
	fout.close();
}


// *****************************************************************************
/// \brief	Write out tip positions of flexible filaments
///
///	\param	tval		time value at which write out is taking place
void ObjectManager::io_writeTipPositions(int tval) {

	// Loop through FEM bodies which this rank owns
	for (auto ib : idxFEM) {

		// Create string and file streams
		std::ofstream fout;
		fout.open(GridUtils::path_str + "/Body_" + std::to_string(iBody[ib].id) + "_TipPositions.out", std::ios::app);
		fout.precision(L_OUTPUT_PRECISION);

		// Write out header
		if (tval == 0)
			fout << "Timestep\tt\tTipX\tTipY\tTipZ" << std::endl;

		// Index of last markers
		int idx = static_cast<int>(iBody[ib].markers.size()) - 1;

		// Write out data
		fout << tval << "\t" << tval * _Grids->dt << "\t" << iBody[ib].markers[idx].position[eXDirection] << "\t"
														  << iBody[ib].markers[idx].position[eYDirection] << "\t"
														  << iBody[ib].markers[idx].position[eZDirection] << std::endl;

		// Close file
		fout.close();
	}
}
