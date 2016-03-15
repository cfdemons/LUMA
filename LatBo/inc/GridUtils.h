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
	static std::vector<double> subtract(std::vector<double> a, std::vector<double> b);			// Function: subtract
	static std::vector<double> add(std::vector<double> a, std::vector<double> b);				// Function: add
	static std::vector<double> vecmultiply(double scalar, std::vector<double> vec);				// Function: multiply
	static std::vector<double> crossprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: crossprod
	static std::vector<double> matrix_multiply(std::vector< std::vector<double> >& A, std::vector<double>& x);	// Function: matrix_multiply

	// Templated functions //

	// Function: vecnorm + overload
	template <typename NumType>
	static NumType vecnorm(NumType a1, NumType a2, NumType a3) {
		return (NumType)sqrt( a1*a1 + a2*a2 + a3*a3 );
	}
	template <typename NumType>
	static NumType vecnorm(NumType a1, NumType a2) {
		return (NumType)sqrt( a1*a1 + a2*a2 );
	}
	
	// Function: upToZero
	template <typename NumType>
	static NumType upToZero(NumType x) {	
		if (x < 0) return 0;
		else return x;		
	};

	// Function: factorial
	template <typename NumType>
	static NumType factorial(NumType n) {
		if (n == 0) 
			return 1;
		else
			return n * GridUtils::factorial(n - 1);
	};

	// LBM-specific utilities
	static size_t getOpposite(size_t direction);	// Function: getOpposite

	// MPI-related utilities
	static bool isOverlapPeriodic(int i, int j, int k, GridObj& pGrid);	// Function: isOverlapPeriodic
	static bool isOnThisRank(int gi, int gj, int gk, GridObj& pGrid);	// Function: isOnThisRank + overloads
	static bool isOnThisRank(int gl, int xyz, GridObj& pGrid);
	static bool hasThisSubGrid(GridObj& pGrid, int RegNum);	// Function: hasThisSubGrid
	// The following supercede the old isOnEdge function to allow for different sized overlaps produced by different refinement levels.
	static bool isOnSenderLayer(double pos_x, double pos_y, double pos_z);		// Is site on any sender layer
	static bool isOnRecvLayer(double pos_x, double pos_y, double pos_z);		// Is site on any recv layer
	static bool isOnSenderLayer(double site_position, std::string dir, std::string maxmin);		// Is site on specified sender layer
	static bool isOnRecvLayer(double site_position, std::string dir, std::string maxmin);		// Is site on speicfied recv layer

	// General Utilities
	static void getGrid(GridObj*& Grids, int level, int region, GridObj*& ptr);	// Function to get pointer to grid in hierarchy

};

