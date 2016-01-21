#pragma once

/** GridUtils Class is a utility class to hold all the general
 *  methods used by the GridObj and others. Everything about this is static
 *  as no need to instantiate it for every grid on a process.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include "definitions.h"

class GridObj;

class GridUtils {

	// Properties //

public:
	static std::ofstream* logfile;			// Handle to output file
	static std::string path_str;            // Static string representing output path

	// Methods //

private:
	// Since class is static make constructors private to prevent instantiation
	GridUtils();
	~GridUtils();

public:
	// IO utilities
	static int createOutputDirectory(std::string path_str);		// Output directory creator

	// Mathematical and numbering utilities
	static std::vector<int> onespace(int min, int max);						// Function: onespace
	static std::vector<double> linspace(double min, double max, int n);		// Function: linspace
	static double vecnorm(double vec[]);										// Function: vecnorm + overloads
	static double vecnorm(double val1, double val2);
	static double vecnorm(double val1, double val2, double val3);
	static double vecnorm(std::vector<double>& vec);
	static std::vector<int> indmapref(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start); // Function: indmapref
	static std::vector<int> revindmapref(int fine_i, int x_start, int fine_j, int y_start, int fine_k, int z_start); // Function: revindmapref
	static double dotprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: dotprod
	static std::vector<double> matrix_multiply(std::vector< std::vector<double> >& A, std::vector<double>& x);	// Function: matrix_multiply

	// LBM-specific utilities
	static size_t getOpposite(size_t direction);	// Function: getOpposite

	// MPI-related utilities
	static bool isOnEdge(unsigned int i, unsigned int j, unsigned int k, GridObj& pGrid);	// Function: isOnEdge
	static bool isOverlapPeriodic(unsigned int i, unsigned int j, unsigned int k,
		GridObj& pGrid,	unsigned int lattice_dir);	// Function: isOverlapPeriodic
	static bool isOnThisRank(unsigned int gi, unsigned int gj, unsigned int gk, GridObj& pGrid);	// Function: isOnThisRank + overloads
	static bool isOnThisRank(unsigned int gl, unsigned int xyz, GridObj& pGrid);
	static bool hasThisSubGrid(GridObj& pGrid, int RegNum);	// Function: hasThisSubGrid
	static bool isOnSenderLayer(double site_position, char dir, char* maxmin);
	static bool isOnRecvLayer(double site_position, char dir, char* maxmin);

	// General Utilities
	static void getGrid(GridObj*& Grids, int level, int region, GridObj*& ptr);

};

