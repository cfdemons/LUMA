// TecplotLiteSorter.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

//#define TECPLOT_DEBUG
#define dims 2
#define output_precision 10
#define num_grids 1
#define max_time 1000
#define time_step 1000
#define num_ranks 8


void sortrows(std::vector<std::vector<double>>& matrix, int col) {

	// Loop over rows with sort
    std::stable_sort(
		matrix.begin(), matrix.end(),

		// Lambda expression to decide whether line goes before another
        [col](const std::vector<double>& row1, const std::vector<double>& row2) {
            return row1[col] < row2[col];
              }
	);
}

void unitTest() {
	
	// Container
	std::vector< std::vector<double> > data;

	// Random data
	data.resize(4);
	data[0].push_back(4); data[0].push_back(5); data[0].push_back(6);
	data[1].push_back(0); data[1].push_back(2); data[1].push_back(8);
	data[2].push_back(1); data[2].push_back(2); data[2].push_back(3);
	data[3].push_back(6); data[3].push_back(7); data[3].push_back(9);

	// Sort
	sortrows(data,3);
	sortrows(data,2);
	sortrows(data,1);

	// Output
	for (int i = 0; i < 4; i++) {
		std::cout << std::endl;
		for (int j = 0 ; j < 3; j++) {
			std::cout << data[i][j] << "\t";
		}
	}


}

int main(int argc, char* argv)
{

#if (dims == 3 && defined TECPLOT_DEBUG)
	int num_vars = 47;
#elif (dims == 2 && defined TECPLOT_DEBUG)
	int num_vars = 27;
#else
	int num_vars = 19;
#endif

	// Files
	std::ofstream tecfile;
	std::ifstream srcfile;
	std::string filename;
	std::string line_in;	// String to store line
	std::istringstream iss;	// Buffer stream
	int line_counter;

#ifdef TECPLOT_DEBUG
	std::string filetype = "tecplotdebug";
#else
	std::string filetype = "io_lite";
#endif

	// Text variables
	std::string time_string;
	std::string variables; variables.resize(num_vars);
	double tmp;
	double type;
	double rank;

	// Variable containers
	std::vector< std::vector<double> > data;
	std::vector<double> row; row.resize(num_vars);

	// Item counter
	unsigned int item_number = 0;

	// Loop over time
	for (int time = 0; time <= max_time; time += time_step) {

		// Reset item numbers and flush stores
		item_number = 0;
		data.clear();
		row.clear(); row.resize(num_vars-1);	// Do not store the rank				

		// Loop over levels
		for (int level = 0; level < num_grids; level++) {

#ifndef TECPLOT_DEBUG
			// Loop over ranks
			for (int rank = 0; rank < num_ranks; rank++)
#endif
			
			{

				// Reset line counter
				line_counter = 0;
				iss.clear(); // Flush the string stream so eofbit is reset
				
				// Try to open file
#ifdef TECPLOT_DEBUG
				filename = "./process/" + filetype + ".Lev" + std::to_string(level) + 
					".Reg0." + std::to_string(time) + ".dat";
#else
				filename = "./process/" + filetype + ".Lev" + std::to_string(level) + ".Reg0.Rnk" + 
					std::to_string(rank) + "." + std::to_string(time) + ".dat";
#endif
				srcfile.open(filename, std::ios::in);
				if (!srcfile.is_open()) continue;


#ifdef TECPLOT_DEBUG
				std::cout << "Reading in Grid Level " << level << " of " << num_grids-1 << 
					", Time " << time << " of " << max_time << std::endl;
#else
				std::cout << "Reading in Rank " << rank << " of " << num_ranks-1 << 
					", Grid Level " << level << " of " << num_grids-1 << 
					", Time " << time << " of " << max_time << std::endl;
#endif

				// Read in one line at a time
				while( !srcfile.eof() ) {

					// Get line and put in buffer
					std::getline(srcfile, line_in, '\n');
					iss.str(line_in);
					iss.seekg(0); // Reset buffer position to start of buffer

					if (line_counter  < 2) {

						// Ignore header line
						line_counter++;
						continue;

					} else if (line_counter == 1) {

						 // Read in time
						time_string = line_in;
						line_counter++;
						continue;

					} else if (line_counter == 2) {

						// Read in variable titles
						variables = line_in;
						line_counter++;
						continue;
					}

					// Else continue...

					// Increment line counter
					line_counter++;

					// Read in rank (omit from data store as it messes up diff comprisons)
					iss >> rank;

					// Read in type
					iss >> type;

					// Don't read in certain site types
					if (type == 2 || type == 3 || type == 5) continue;

					// Increment item number store
					item_number++;

					// Store type
					row[0] = type;
					
					// Read in rest of data
					for (int item = 1; item < num_vars-1; item++) {
						iss >> tmp; row[item] = tmp;
					}
					
					// Assign to data store
					data.push_back(row);

				}

				// Reached end of file so close file
				srcfile.close();

				// Delete erroneous space at end of file
				data.erase(data.end()-1,data.end()); item_number--;


			}

		}

		// Break if the time step doesn't exist
		if (item_number == 0) continue;


		// Sort items
		std::cout << "Sorting data..." << std::endl;

		// X Y Z variables are in columns 1 2 3
		sortrows(data,3);
		sortrows(data,2);
		sortrows(data,1);



		// Write Data //
		std::cout << "Writing data..." << std::endl;

		// Open file
		tecfile.open("./out/tecplot." + std::to_string(time) + ".dat");

		tecfile.precision(output_precision);
		tecfile.setf(std::ios::fixed);
		tecfile.setf(std::ios::showpoint);

		// Write out tecplot header
		tecfile << "TITLE = All grid quantities excluding refined regions" << std::endl;
		tecfile << "FILETYPE = FULL" << std::endl;
#ifdef TECPLOT_DEBUG
		tecfile << "VARIABLES = \"TYPE\" \"X\" \"Y\" \"Z\" \"F\" \"FEQ\" \"RHO\" \"UX\" \"UY\" \"UZ\"" << std::endl;
#else
		tecfile << "VARIABLES = \"TYPE\" \"X\" \"Y\" \"Z\" \"RHO\" \"UX\" \"UY\" \"UZ\" \"TA_RHO\" \"TA_UX\" \"TA_UY\" \"TA_UZ\" "<<
			"\"TA_UXUX\" \"TA_UXUY\" \"TA_UXUZ\" \"TA_UYUY\" \"TA_UYUZ\" \"TA_UZUZ\"" << std::endl;
#endif
		tecfile << "ZONE" << std::endl;
		tecfile << "I = " << item_number << std::endl;
		tecfile << "DATAPACKING = POINT" << std::endl;
		tecfile << "SOLUTIONTIME = " << std::to_string(time) << std::endl;


		// Write out items (rank not included)
		for (size_t i = 0; i < data.size(); i++) {
			for (int j = 0; j < num_vars-1; j++) {
				tecfile << data[i][j] << "\t";
			}
			tecfile << std::endl;
		}
		
		// Time step complete
		tecfile.close();
	}


	return 0;
}

