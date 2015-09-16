#pragma once

/** GridUtils Class is a utility class to hold all the general 
 *  methods used by the GridObj and others
 */

#include <vector>
#include <iostream>
#include <fstream>

class GridUtils {

public:
	std::ofstream* logfile;		// Handle to output file


// Constructor and Destructor
public:
	GridUtils();
	~GridUtils();

	// Public Utility Methods
	std::vector<int> onespace(int min, int max);						// Function: onespace
	std::vector<double> linspace(double min, double max, int n);		// Function: linspace
	double vecnorm(double vec[]);										// Function: vecnorm + overloads
	double vecnorm(double val1, double val2);
	double vecnorm(double val1, double val2, double val3);
	double vecnorm(std::vector<double>& vec);
	std::vector<int> indmapref(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start); // Function: indmapref
	double dotprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: dotprod
	std::vector<double> matrix_multiply(std::vector< std::vector<double> >& A, std::vector<double>& x);	// Function: matrix_multiply
	size_t getOpposite(size_t direction);	// Function: getOpposite

	// Set log file method
	void setLogFile(std::ofstream* logfile);
};

