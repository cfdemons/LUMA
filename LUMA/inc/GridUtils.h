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

#pragma once

/// \enum eCartesianDirection
/// \brief Enumeration for directional options.
enum eCartesianDirection
{
	eXDirection,	///< X-direction
	eYDirection,	///< Y-direction
	eZDirection		///< Z-direction
};

/// \enum eMinMax
/// \brief	Enumeration for minimum and maximum.
///
///			Some utility methods need to know whether they should be looking
///			at or for a maximum or minimum edge of a grid so we use this 
///			enumeration to specify.
enum eMinMax
{
	eMinimum,		///< Minimum
	eMaximum		///< Maximum
};

#include "stdafx.h"
#include "definitions.h"
#include "GridObj.h"
#include "hdf5luma.h"


/// \brief	Grid utility class.
///
///			Class provides grid utilities including commonly used logical tests. 
///			This is a static class and so there is no need to instantiate it.
class GridUtils {

	// Properties //

public:
	static std::ofstream* logfile;			///< Handle to output file
	static std::string path_str;            ///< Static string representing output path
	static const int dir_reflect[L_dims * 2][L_nVels];	///< Array with hardcoded direction numbering for specular reflection

	// Methods //

private:
	/// Private constructor since class is static
	GridUtils();
	/// Private destructor
	~GridUtils();

public:
	// IO utilities
	static int createOutputDirectory(std::string path_str);		// Output directory creator

	// Mathematical and numbering utilities
	static std::vector<int> onespace(int min, int max);						// Function: onespace
	static std::vector<double> linspace(double min, double max, int n);		// Function: linspace
	static double vecnorm(double vec[]);									// Function: vecnorm + overloads
	static double vecnorm(double val1, double val2);
	static double vecnorm(double val1, double val2, double val3);
	static double vecnorm(std::vector<double> vec);
	static std::vector<int> getFineIndices(
		int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start); // Function: getFineIndices
	static std::vector<int> getCoarseIndices(
		int fine_i, int x_start, int fine_j, int y_start, int fine_k, int z_start); // Function: getCoarseIndices
	static double indexToPosition(int index, double dx);	// Function: indexToPosition
	static double dotprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: dotprod
	static std::vector<double> subtract(std::vector<double> a, std::vector<double> b);			// Function: subtract
	static std::vector<double> add(std::vector<double> a, std::vector<double> b);				// Function: add
	static std::vector<double> vecmultiply(double scalar, std::vector<double> vec);				// Function: multiply
	static std::vector<double> crossprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: crossprod
	static std::vector<double> matrix_multiply(const std::vector< std::vector<double> >& A, const std::vector<double>& x);	// Function: matrix_multiply

	// LBM-specific utilities
	static int getOpposite(int direction);	// Function: getOpposite
	static void getGrid(GridObj*& Grids, int level, int region, GridObj*& ptr);		// Function to get pointer to grid in hierarchy

	// MPI-related utilities
	static bool isOverlapPeriodic(int i, int j, int k, const GridObj& pGrid);	// Function: isOverlapPeriodic
	static bool isOnThisRank(int gi, int gj, int gk, const GridObj& pGrid);	// Function: isOnThisRank + overloads
	static bool isOnThisRank(int gl, enum eCartesianDirection xyz, const GridObj& pGrid);
	static bool hasThisSubGrid(const GridObj& pGrid, int RegNum);	// Function: hasThisSubGrid
	// The following supercede the old isOnEdge function to allow for different sized overlaps produced by different refinement levels.
	static bool isOnSenderLayer(double pos_x, double pos_y, double pos_z);		// Is site on any sender layer
	static bool isOnRecvLayer(double pos_x, double pos_y, double pos_z);		// Is site on any recv layer
	static bool isOnSenderLayer(double site_position, enum eCartesianDirection xyz, enum eMinMax minmax);	// Is site on specified sender layer
	static bool isOnRecvLayer(double site_position, enum eCartesianDirection xyz, enum eMinMax minmax);		// Is site on speicfied recv layer
	static bool isOffGrid(int i, int j, int k, GridObj& g);						// Is site off supplied grid


	// Templated functions //

	/// \brief Computes the L2-norm.
	/// \param a1 first component of the vector
	/// \param a2 second component of the vector
	/// \param a3 third component of the vector
	/// \return NumType scalar quantity
	template <typename NumType>
	static NumType vecnorm(NumType a1, NumType a2, NumType a3) {
		return (NumType)sqrt( a1*a1 + a2*a2 + a3*a3 );
	}
	/// \brief Computes the L2-norm.
	/// \param a1 first component of the vector
	/// \param a2 second component of the vector
	/// \return NumType scalar quantity
	template <typename NumType>
	static NumType vecnorm(NumType a1, NumType a2) {
		return (NumType)sqrt( a1*a1 + a2*a2 );
	}
	
	/// \brief	Rounds a negative value up to zero.
	///
	///			If value is positive, return the value unchanged.
	///
	/// \param x value to be rounded
	/// \return NumType rounded value
	template <typename NumType>
	static NumType upToZero(NumType x) {	
		if (x < static_cast<NumType>(0)) return static_cast<NumType>(0);
		else return x;		
	};

	/// \brief	Rounds a value greater than a limit down to this value.
	///
	///			If value is less than or equal to the limit, return the value 
	///			unchanged.
	///
	/// \param x value to be rounded
	/// \param limit value to be rounded down to
	/// \return NumType rounded value
	template <typename NumType>
	static NumType downToLimit(NumType x, NumType limit) {	
		if (x > limit) return limit;
		else return x;		
	};

	/// \brief	Computes the factorial of the supplied value.
	///
	///			If n == 0 then returns 1.
	///
	/// \param n factorial
	/// \return NumType n factorial
	template <typename NumType>
	static NumType factorial(NumType n) {
		if (n == static_cast<NumType>(0)) 
			return static_cast<NumType>(1);
		else
			return n * GridUtils::factorial(n - 1);
	};

	/// \brief	Performs a strided memcpy.
	///
	///			Memcpy() is designed to copy blocks of contiguous memory.
	///			Strided copy copies a pattern of contiguous blocks.
	///
	/// \param dest pointer to start of destination memory
	/// \param src pointer to start of source memory
	/// \param block size of contiguous block
	/// \param offset offset from the start of the soruce array
	/// \param stride number of elements between start of first block and start of second
	/// \param count number of blocks in pattern
	/// \param buf_offset offset from start of destination buffer to start writing. Default is zero if not supplied.
	template <typename NumType>
	static void stridedCopy(NumType *dest, NumType *src, size_t block, 
		size_t offset, size_t stride, size_t count,
		size_t buf_offset = 0) {

		for (size_t i = 0; i < count; i++) {
			memcpy(dest + buf_offset, src + offset + (stride * i), block * sizeof(NumType));
			buf_offset += block;
		}
	};


	// ************************************************************************
	/// \brief	Maps global indices to local indices.
	///
	///			Takes a vector container and populates it with the 
	///			local indices where the supplied global site can be found
	///			on the grid supplied. If global indicies are not found on the
	///			supplied grid then local index of -1 is returned.
	///
	/// \param i global index
	/// \param j global index
	/// \param k global index
	/// \param g grid on which local indices are required
	/// \param[out] locals vector container for local indices
	template <typename NumType>
	static void global_to_local(int i, int j, int k, GridObj* g, std::vector<NumType>& locals) {

		// Clear array if not empty
		if (locals.size()) locals.clear();

		// Find indices in arrays
		auto loc_i = std::find(g->XInd.begin(), g->XInd.end(), i);
		auto loc_j = std::find(g->YInd.begin(), g->YInd.end(), j);
		auto loc_k = std::find(g->ZInd.begin(), g->ZInd.end(), k);

		// Put indices in locals
		if (loc_i != g->XInd.end()) locals.push_back( static_cast<NumType>(std::distance(g->XInd.begin(),loc_i)) );
		else locals.push_back(-1);
		if (loc_j != g->YInd.end()) locals.push_back( static_cast<NumType>(std::distance(g->YInd.begin(),loc_j)) );
		else locals.push_back(-1);
		if (loc_k != g->ZInd.end()) locals.push_back( static_cast<NumType>(std::distance(g->ZInd.begin(),loc_k)) );
		else locals.push_back(-1);

		return;	// Directions not on the grid are returned with index -1


	}

	// ************************************************************************
	/// \brief	Maps local indices to global indices.
	///
	///			Takes a vector container and populates it with the 
	///			global indices of the supplied local site
	///
	/// \param i local index
	/// \param j local index
	/// \param k local index
	/// \param g grid on which global indices are required
	/// \param[out] globals vector container for global indices
	template <typename NumType>
	static void local_to_global(int i, int j, int k, GridObj* g, std::vector<NumType>& globals) {

		// Clear array if not empty
		if (globals.size()) globals.clear();

		// Put indices in globals
		globals.push_back( static_cast<NumType>(g->XInd[i]) );
		globals.push_back( static_cast<NumType>(g->YInd[j]) );
		globals.push_back( static_cast<NumType>(g->ZInd[k]) );

		return;

	}




};

