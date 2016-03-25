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
	static std::vector<int> getFineIndices(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start); // Function: getFineIndices
	static std::vector<int> getCoarseIndices(int fine_i, int x_start, int fine_j, int y_start, int fine_k, int z_start); // Function: getCoarseIndices
	static double dotprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: dotprod
	static std::vector<double> subtract(std::vector<double> a, std::vector<double> b);			// Function: subtract
	static std::vector<double> add(std::vector<double> a, std::vector<double> b);				// Function: add
	static std::vector<double> vecmultiply(double scalar, std::vector<double> vec);				// Function: multiply
	static std::vector<double> crossprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: crossprod
	static std::vector<double> matrix_multiply(std::vector< std::vector<double> >& A, std::vector<double>& x);	// Function: matrix_multiply

	// LBM-specific utilities
	static size_t getOpposite(size_t direction);	// Function: getOpposite
	static void getGrid(GridObj*& Grids, int level, int region, GridObj*& ptr);		// Function to get pointer to grid in hierarchy

	// MPI-related utilities
	static bool isOverlapPeriodic(int i, int j, int k, GridObj& pGrid);	// Function: isOverlapPeriodic
	static bool isOnThisRank(int gi, int gj, int gk, GridObj& pGrid);	// Function: isOnThisRank + overloads
	static bool isOnThisRank(int gl, int xyz, GridObj& pGrid);
	static bool hasThisSubGrid(GridObj& pGrid, int RegNum);	// Function: hasThisSubGrid
	// The following supercede the old isOnEdge function to allow for different sized overlaps produced by different refinement levels.
	static bool isOnSenderLayer(double pos_x, double pos_y, double pos_z);		// Is site on any sender layer
	static bool isOnRecvLayer(double pos_x, double pos_y, double pos_z);		// Is site on any recv layer
	static bool isOnSenderLayer(double site_position, int dir, int maxmin);		// Is site on specified sender layer
	static bool isOnRecvLayer(double site_position, int dir, int maxmin);		// Is site on speicfied recv layer


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


	// ************************************************************************
	// Map global to local indices where locals is an empty vector container
	template <typename NumType>
	static void global_to_local(int i, int j, int k, GridObj* g, std::vector<NumType>& locals) {

		// Find indices in arrays
		auto loc_i = std::find(g->XInd.begin(), g->XInd.end(), i);
		auto loc_j = std::find(g->YInd.begin(), g->YInd.end(), j);
		auto loc_k = std::find(g->ZInd.begin(), g->ZInd.end(), k);

		// Put indices in locals
		if (loc_i != g->XInd.end()) locals.push_back( static_cast<NumType>(std::distance(g->XInd.begin(),loc_i)) );
		if (loc_j != g->YInd.end()) locals.push_back( static_cast<NumType>(std::distance(g->YInd.begin(),loc_j)) );
		if (loc_k != g->ZInd.end()) locals.push_back( static_cast<NumType>(std::distance(g->ZInd.begin(),loc_k)) );

		return;	// Will return vector with size() < 3 if global index not found on this grid


	}

	// ************************************************************************
	// Map local to global indices where globals is an empty vector container
	template <typename NumType>
	static void local_to_global(int i, int j, int k, GridObj* g, std::vector<NumType>& globals) {

		// Put indices in globals
		globals.push_back( static_cast<NumType>(g->XInd[i]) );
		globals.push_back( static_cast<NumType>(g->YInd[j]) );
		globals.push_back( static_cast<NumType>(g->ZInd[k]) );

		return;

	}




};

