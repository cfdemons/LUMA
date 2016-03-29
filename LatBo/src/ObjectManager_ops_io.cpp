#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"
#include "../inc/MpiManager.h"
#include <iostream>
#include <fstream>
#include <sstream>

// ***************************************************************************************************
// Routine to write out the coordinates of IBBodies at a given time step
void ObjectManager::io_write_body_pos(int timestep) {

	for (size_t ib = 0; ib < iBody.size(); ib++) {


			// Open file for given time step
			std::ofstream jout;
			jout.open(GridUtils::path_str + "/Body_" + std::to_string(ib) + "_position_" + std::to_string(timestep) + 
				"_rank" + std::to_string(MpiManager::my_rank) + ".out", std::ios::out);
			jout << "x" + std::to_string(timestep) + ", y" + std::to_string(timestep) + ", z" << std::endl;

			// Write out position
			for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
#if (dims == 3)
				jout << iBody[ib].markers[i].position[0] << ", " << iBody[ib].markers[i].position[1] << ", " << iBody[ib].markers[i].position[2] << std::endl;
#else
				jout << iBody[ib].markers[i].position[0] << ", " << iBody[ib].markers[i].position[1] << ", " << 0.0 << std::endl;
#endif
			}
			jout.close();

	}

}


// ***************************************************************************************************
// Routine to write out the coordinates of IBbodies at a given time step
void ObjectManager::io_write_lift_drag(int timestep) {

	for (size_t ib = 0; ib < iBody.size(); ib++) {


			// Open file for given time step
			std::ofstream jout;
			jout.open(GridUtils::path_str + "/Body_" + std::to_string(ib) + "_LD_" + std::to_string(timestep) + "_rank" + std::to_string(MpiManager::my_rank) + ".out", std::ios::out);
			jout << "L" + std::to_string(timestep) + ", D" + std::to_string(timestep) << std::endl;

			// Sum variables
			double Lsum = 0.0, Dsum = 0.0;

			// Compute lift and drag
			for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
				jout << iBody[ib].markers[i].force_xyz[0] << ", " << iBody[ib].markers[i].force_xyz[1] << std::endl;
				Lsum += iBody[ib].markers[i].force_xyz[0];
				Dsum += iBody[ib].markers[i].force_xyz[1];
			}

			jout << "Totals = " << std::endl;
			jout << Lsum << ", " << Dsum << std::endl;
			jout.close();

	}

}

// ***************************************************************************************************
void ObjectManager::io_restart(bool IO_flag, int level) {

	if (IO_flag) {

		// Output stream
		std::ofstream file;

		if (MpiManager::my_rank == 0 && level == 0) { // Overwrite as first to write
			file.open(GridUtils::path_str + "/restart_IBBody.out", std::ios::out);
		} else if (level == 0) { // Append
			file.open(GridUtils::path_str + "/restart_IBBody.out", std::ios::out | std::ios::app);
		} else { // Must be a subgrid which doesn't own any IB-bodies so return
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

				// Old positions of each marker
				file	<< iBody[b].markers[m].position_old[0] << "\t"
						<< iBody[b].markers[m].position_old[1] << "\t"
						<< iBody[b].markers[m].position_old[2] << "\t";

			}

		}

		// Close file
		file.close();


	} else {

		// Input stream
		std::ifstream file;

		// Only level 0 grids can own IB-bodies
		if (level == 0) {
			file.open(GridUtils::path_str + "/restart_IBBody.out", std::ios::in);
		}

		if (!file.is_open()) {
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "Error opening IBM restart file. Exiting." << std::endl;
			exit(EXIT_FAILURE);
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
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "Number of IBM bodies does not match the number specified in the restart file. Exiting." << std::endl;
			exit(EXIT_FAILURE);
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
				std::cout << "Error: See Log File" << std::endl;
				*GridUtils::logfile << "Number of IBM markers does not match the number specified for body " <<
					b << " in the restart file. Exiting." << std::endl;
				exit(EXIT_FAILURE);
			}

			// Read in marker data
			for (m = 0; m < num_mark; m++) {

				// Positions of each marker
				iss		>> iBody[b].markers[m].position[0]
						>> iBody[b].markers[m].position[1]
						>> iBody[b].markers[m].position[2];

				// Old positions of each marker
				iss		>> iBody[b].markers[m].position_old[0]
						>> iBody[b].markers[m].position_old[1]
						>> iBody[b].markers[m].position_old[2];

			}

		}

		// Close file
		file.close();

	}

}

// ***************************************************************************************************
// Routine to write out the vtk (position) for each IB body at time step t (current capability is for unclosed objects only)
void ObjectManager::io_vtk_IBwriter(double tval) {

    // Loop through each iBody
    for (size_t ib = 0; ib < iBody.size(); ib++) {

        // Create file name then output file stream
        std::stringstream fileName;
        fileName << GridUtils::path_str + "/vtk_IBout.Body" << ib << "." << (int)tval << ".vtk";

        std::ofstream fout;
        fout.open( fileName.str().c_str() );

        // Add header information
        fout << "# vtk DataFile Version 3.0f\n";
        fout << "IB Output for body ID " << ib << " at time t = " << (int)tval << "\n";
        fout << "ASCII\n";
        fout << "DATASET POLYDATA\n";


        // Write out the positions of each Lagrange marker
        fout << "POINTS " << iBody[ib].markers.size() << " float\n";
        for (size_t i = 0; i < iBody[ib].markers.size(); i++) {

#if (dims == 3)
				fout << iBody[ib].markers[i].position[0] << " " << iBody[ib].markers[i].position[1] << " " << iBody[ib].markers[i].position[2] << std::endl;
#else
				fout << iBody[ib].markers[i].position[0] << " " << iBody[ib].markers[i].position[1] << " " << 1.0 << std::endl; // z = 1.0 as fluid ORIGIN is at z = 1.0
#endif
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


// ***************************************************************************************************
// Routine to read in point cloud data in tab separated, 3-column format from the input directory
void ObjectManager::readInPCData(PCpts* _PCpts) {

#ifdef BFL_ON

	// Temporary variables
	double tmp_x, tmp_y, tmp_z;

	// Open input file
	std::ifstream file;
	file.open("./input/bfl_input.in", std::ios::in);
	
	// Handle failure to open
	if (!file.is_open()) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Error opening BFL input file. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Get grid pointer
	GridObj* g;	GridUtils::getGrid(_Grids, bfl_on_grid_lev, bfl_on_grid_reg, g);

	// Loop over lines in file
	while (!file.eof()) {

		// Read in one line of file at a time
		std::string line_in;	// String to store line in
		std::istringstream iss;	// Buffer stream to store characters

		// Get line up to new line separator and put in buffer
		std::getline(file,line_in,'\n');
		iss.str(line_in);	// Put line in the buffer
		iss.seekg(0);		// Reset buffer position to start of buffer

		// Add coordinates to data store
		iss >> tmp_x;
		iss >> tmp_y;
		iss >> tmp_z;
			
		_PCpts->x.push_back(tmp_x);
		_PCpts->y.push_back(tmp_y);
			
		// If running a 2D calculation, only read in x and y coordinates and force z coordinates to match the domain
#if (dims == 3)
		_PCpts->z.push_back(tmp_z);
#else
		_PCpts->z.push_back(0);
#endif

		}
	file.close();

	// Warn of no data case
	if (_PCpts->x.empty() || _PCpts->y.empty() || _PCpts->z.empty()) {
		std::cout << "Warning: Rank " + std::to_string(MpiManager::my_rank) + " contains no BFL data." << std::endl;
	}

	
	// Rescale coordinates and shift
	double scale_factor = bfl_length_x / 
		std::fabs(*std::max_element(_PCpts->x.begin(), _PCpts->x.end()) - *std::min_element(_PCpts->x.begin(), _PCpts->x.end()));
	double shift_x =  std::floor( start_bfl_x - scale_factor * *std::min_element(_PCpts->x.begin(), _PCpts->x.end()) );
	double shift_y =  std::floor( start_bfl_y - scale_factor * *std::min_element(_PCpts->y.begin(), _PCpts->y.end()) );
	double shift_z =  std::floor( start_bfl_z - scale_factor * *std::min_element(_PCpts->z.begin(), _PCpts->z.end()) );

	// Apply
	for (size_t a = 0; a < _PCpts->x.size(); a++) {
		_PCpts->x[a] *= scale_factor; _PCpts->x[a] += shift_x;
		_PCpts->y[a] *= scale_factor; _PCpts->y[a] += shift_y;
		_PCpts->z[a] *= scale_factor; _PCpts->z[a] += shift_z;
	}

	// Exclude points which are not on this rank
	int a = 0;
	do {

		// Delete points not on this rank
		if ( GridUtils::isOnThisRank( 
			BFLBody::getVoxInd(_PCpts->x[a]), BFLBody::getVoxInd(_PCpts->y[a]), BFLBody::getVoxInd(_PCpts->z[a]), *g) ) a++;
		else {
			_PCpts->x.erase(_PCpts->x.begin() + a);
			_PCpts->y.erase(_PCpts->y.begin() + a);
			_PCpts->z.erase(_PCpts->z.begin() + a);
		}

	} while (a < (int)_PCpts->x.size());


	// Write out the points remaining in for debugging purposes
#ifdef BFL_DEBUG
	if (!_PCpts->x.empty()) {
		std::ofstream fileout;
		fileout.open(GridUtils::path_str + "/BFLpts_rank" + std::to_string(MpiManager::my_rank) + ".out",std::ios::out);
		for (size_t i = 0; i < _PCpts->x.size(); i++) {
			fileout << std::to_string(_PCpts->x[i]) + '\t' + std::to_string(_PCpts->y[i]) + '\t' + std::to_string(_PCpts->z[i]);
			fileout << std::endl;
		}
		fileout.close();
	}
#endif


#endif
}

// ***************************************************************************************************
// Routine to read in point cloud data in tab separated, 3-column format from the input directory
void ObjectManager::readInPointData(PCpts* _PCpts) {

#ifdef SOLID_FROM_FILE

	// Temporary variables
	double tmp_x, tmp_y, tmp_z;

	// Open input file
	std::ifstream file;
	file.open("./input/pointcloud.in", std::ios::in);
	
	// Handle failure to open
	if (!file.is_open()) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Error opening Point Cloud input file. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Get grid pointer
	GridObj* g = NULL;	GridUtils::getGrid(_Grids, object_on_grid_lev, object_on_grid_reg, g);
	if (g == NULL) return; // If this grid does not exist, then exit

	// Loop over lines in file
	while (!file.eof()) {

		// Read in one line of file at a time
		std::string line_in;	// String to store line in
		std::istringstream iss;	// Buffer stream to store characters

		// Get line up to new line separator and put in buffer
		std::getline(file,line_in,'\n');
		iss.str(line_in);	// Put line in the buffer
		iss.seekg(0);		// Reset buffer position to start of buffer

		// Add coordinates to data store
		iss >> tmp_x;
		iss >> tmp_y;
		iss >> tmp_z;
			
		_PCpts->x.push_back(tmp_x);
		_PCpts->y.push_back(tmp_y);
			
		// If running a 2D calculation, only read in x and y coordinates and force z coordinates to match the domain
#if (dims == 3)
		_PCpts->z.push_back(tmp_z);
#else
		_PCpts->z.push_back(0);
#endif

		}
	file.close();

	// Warn of no data case
	if (_PCpts->x.empty() || _PCpts->y.empty() || _PCpts->z.empty()) {
		std::cout << "Warning: Rank " + std::to_string(MpiManager::my_rank) + " contains no object data." << std::endl;
	}

	
	// Rescale coordinates and shift
	double scale_factor = object_length_x / 
		std::fabs(*std::max_element(_PCpts->x.begin(), _PCpts->x.end()) - *std::min_element(_PCpts->x.begin(), _PCpts->x.end()));
	double shift_x =  std::floor( start_object_x - scale_factor * *std::min_element(_PCpts->x.begin(), _PCpts->x.end()) );
	double shift_y =  std::floor( start_object_y - scale_factor * *std::min_element(_PCpts->y.begin(), _PCpts->y.end()) );
	double shift_z =  std::floor( start_object_z - scale_factor * *std::min_element(_PCpts->z.begin(), _PCpts->z.end()) );

	// Apply
	for (size_t a = 0; a < _PCpts->x.size(); a++) {
		_PCpts->x[a] *= scale_factor; _PCpts->x[a] += shift_x;
		_PCpts->y[a] *= scale_factor; _PCpts->y[a] += shift_y;
		_PCpts->z[a] *= scale_factor; _PCpts->z[a] += shift_z;
	}


	// Write out
	/*std::ofstream fileout;
	fileout.open("./Out_" + std::to_string(MpiManager::my_rank) + ".out",std::ios::out);
	fileout << _PCpts->x.size() << " points read in." << std::endl;
	for (size_t a = 0; a < _PCpts->x.size(); a++) {
		fileout << _PCpts->x[a] << " " << _PCpts->y[a] << " " << _PCpts->z[a] << std::endl;
	}
	fileout.close();*/

	// Declare local variables
	int local_i, local_j, local_k;

	// Exclude points which are not on this rank
	int a = 0;
	do {

		local_i = BFLBody::getVoxInd(_PCpts->x[a]);
		local_j = BFLBody::getVoxInd(_PCpts->y[a]);
		local_k = BFLBody::getVoxInd(_PCpts->z[a]);

		// Delete points not on this rank
		if ( GridUtils::isOnThisRank( local_i, local_j, local_k, *g) ) {

			// Update Typing Matrix
			if ( g->LatTyp(local_i, local_j, local_k, g->YInd.size(), g->ZInd.size()) == 1 ) g->LatTyp(local_i, local_j, local_k, g->YInd.size(), g->ZInd.size()) = 0;

			// Increment counter
			a++;


		} else {

			// Delete point from store as not on this rank
			_PCpts->x.erase(_PCpts->x.begin() + a);
			_PCpts->y.erase(_PCpts->y.begin() + a);
			_PCpts->z.erase(_PCpts->z.begin() + a);

		}

	} while (a < (int)_PCpts->x.size());


#endif
}