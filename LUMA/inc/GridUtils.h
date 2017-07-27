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

#ifndef GRIDUTILS_H
#define GRIDUTILS_H

#include "stdafx.h"
#include "GridObj.h"

/// \brief	Grid utility class.
///
///			Class provides grid utilities including commonly used logical tests. 
///			This is a static class and so there is no need to instantiate it.
class GridUtils {

	// Properties //

public:
	static std::ofstream* logfile;			///< Handle to output file
	static std::string path_str;            ///< Static string representing output path
	static const int dir_reflect[L_DIMS][L_NUM_VELS];	///< Array with hardcoded direction numbering for specular reflection
	static const int dir_opposites[L_NUM_VELS];			///< Array with hardcoded direction numbering for bounce-back opposites

	// Methods //

private:
	/// Private constructor since class is static
	GridUtils();
	/// Private destructor
	~GridUtils();

public:
	// IO utilities
	static void createOutputDirectory(std::string path_str);		// Output directory creator
	static void readVelocityFromFile(std::string path_str, std::vector<double>& x_coord, std::vector<double>& y_coord, std::vector<double>& z_coord, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& uz);  //Reads coordinates and velocity data from file_name. Stores the coordinates of each point in the vectors x, y and z and the velocity components in the vectors ux, uy and uz. It expects the file to have a column for uz even with L_DIMS = 2 

	// Mathematical and numbering utilities
	static std::vector<int> onespace(int min, int max);						// Function: onespace
	static std::vector<double> linspace(double min, double max, int n);		// Function: linspace
	static double vecnorm(double vec[L_DIMS]);								// Function: vecnorm + overloads
	static double vecnorm(double val1, double val2);
	static double vecnorm(double val1, double val2, double val3);
	static double vecnorm(std::vector<double> vec);
	static std::vector<int> getFineIndices(
		int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start); // Function: getFineIndices
	static std::vector<int> getCoarseIndices(
		int fine_i, int x_start, int fine_j, int y_start, int fine_k, int z_start);				// Function: getCoarseIndices
	static double dotprod(std::vector<double> vec1, std::vector<double> vec2);					// Function: dotprod
	static std::vector<double> subtract(std::vector<double> a, std::vector<double> b);			// Function: subtract
	static std::vector<double> add(std::vector<double> a, std::vector<double> b);				// Function: add
	static std::vector<double> vecmultiply(double scalar, std::vector<double> vec);				// Function: multiply
	static std::vector<double> crossprod(std::vector<double> vec1, std::vector<double> vec2);	// Function: crossprod
	static std::vector<double> matrix_multiply(const std::vector< std::vector<double> >& A, const std::vector<double>& x);	// Function: matrix_multiply
	static std::vector<double> divide(std::vector<double> vec1, double scalar);					// Divide vector by a scalar

	// LBM-specific utilities
	static int getOpposite(int direction);	// Function: getOpposite
	static int getReflect(int direction, eCartesianDirection plane);
	static void getGrid(int level, int region, GridObj*& ptr);							// Wrapper using default hierarchy to get grid pointer
	static void getGrid(GridObj* const Grids, int level, int region, GridObj*& ptr);	// Function to get pointer to grid in hierarchy
	static void getFinestGrid(GridObj*& ptr);											// Get a pointer to finest grid in hierarchy
	static double normaliseToLink(double value, int v);									// Normalise value wrt to the lattice link length
	static eType setBCPrecedence(eType currentBC, eType desiredBC);						// Determine BC based on any existing BC

	// MPI-related utilities
	static bool isOverlapPeriodic(int i, int j, int k, GridObj const & pGrid);		// Is this halo periodically connected to neighbour
	static bool isOnThisRank(double x, double y, double z, 
		eLocationOnRank *loc = nullptr, GridObj const * const grid = nullptr,
		std::vector<int> *pos = nullptr);											// Is a site on this MPI rank nad if so, where is it?
	static bool isOnThisRank(double xyz, eCartesianDirection dir, 
		eLocationOnRank *loc = nullptr, GridObj const * const grid = nullptr, int *pos = nullptr);
	static bool intersectsRefinedRegion(GridObj const & pGrid, int RegNum);	// Does the refined region interesect the current rank.
	// The following supercede the old isOnEdge function to allow for different sized overlaps produced by different refinement levels.
	static bool isOnSenderLayer(double pos_x, double pos_y, double pos_z);			// Is site on any sender layer
	static bool isOnRecvLayer(double pos_x, double pos_y, double pos_z);			// Is site on any recv layer
	static bool isOnSenderLayer(double site_position, eCartMinMax edge);			// Is site on specified sender layer
	static bool isOnRecvLayer(double site_position, eCartMinMax edge);				// Is site on specified recv layer
	static int getMpiDirection(int offset_vector[]);								// Get MPI direction from vector
	static int safeGetRank();														// Parallel/Serial safe method to get rank

	// Coordinate Management
	static bool isOffGrid(int i, int j, int k, GridObj const * const g);					// Is site off supplied grid
	static void getEnclosingVoxel(double x, double y, double z, GridObj const * const g, std::vector<int> *ijk);	// Take a position and a get a local ijk
	static void getEnclosingVoxel(double x, GridObj const * const g, eCartesianDirection dir, int *ijk);
	static bool isOnTransitionLayer(double pos_x, double pos_y, double pos_z, GridObj const * const grid);	// Is site on any TL to upper
	static bool isOnTransitionLayer(double position, eCartMinMax edge, GridObj const * const grid);			// Is site on specified TL to upper
	static bool isWithinDomainWall(double posX, double posY, double posZ, std::vector<int>* normalVector = nullptr);
	static bool isWithinDomainWall(double posX, double posY, double posZ,
		std::vector<int>* normalVector, eCartesianDirection* normalDirection, unsigned int* edgeCount);	// Is site within a wall region

	// Templated functions //

	/// \brief	Multiplies a scalar by a vector.
	/// \param	vec1	a vector double.
	/// \param	vec2	another vector double.
	/// \return a vector which is an elementwise mulitplication of the two vectors.
	template <typename T, typename T1, typename Ret>
	static std::vector<Ret> vecmultiply(std::vector<T> vec1, std::vector<T1> vec2)
	{
		assert(vec1.size() == vec2.size());
		std::vector<Ret> result(vec1.size());
		for (size_t i = 0; i < vec1.size(); i++)
		{
			result[i] = vec1[i] * vec2[i];
		}
		return result;
	}

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
			return n * factorial(n - 1);
	};

	/// \brief	Performs a strided memcpy.
	///
	///			Memcpy() is designed to copy blocks of contiguous memory.
	///			Strided copy copies a pattern of contiguous blocks.
	///
	/// \param dest		pointer to start of destination memory.
	/// \param src		pointer to start of source memory.
	/// \param block	size of contiguous block.
	/// \param offset	offset from the start of the soruce array.
	/// \param stride	number of elements between start of first block and start of second.
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

	/// \brief	Function capable of extrapolating any quantity.
	///
	///	\param	grid		grid on which quantity resides.
	/// \param	quantity	quantity to search.
	/// \param	direction	direction in which to source extrapolation data.
	///	\param	order		order of extrapolation.
	/// \param	i			target site x index.
	/// \param	j			target site y index.
	/// \param	k			target site z index.
	/// \param	p			target site extra index.
	/// \param	max			extra index max size for flattening.
	///	\returns			extrapolated value.
	template <typename NumType>
	static NumType extrapolate(GridObj const & grid, IVector<NumType> const & quantity,
		const std::vector<int> direction, int order, 
		const int i, const int j, const int k,
		const int p = NULL, const int max = 1)
	{
		int status = 0;
		NumType result;

		if (order == 0)
		{
			int i1 = i + direction[eXDirection];
			int j1 = j + direction[eYDirection];
			int k1 = k + direction[eZDirection];

			if (GridUtils::isOffGrid(i1, j1, k1, &grid))
				status = -1;
			else
				result = quantity[flatten(i1, j1, k1, p, max, grid)];

		}
		else if (order == 1)
		{
			int i1 = i + direction[eXDirection];
			int i2 = i1 + direction[eXDirection];
			int j1 = j + direction[eYDirection];
			int j2 = j1 + direction[eYDirection];
			int k1 = k + direction[eZDirection];
			int k2 = k1 + direction[eZDirection];

			if (GridUtils::isOffGrid(i1, j1, k1, &grid) || GridUtils::isOffGrid(i2, j2, k2, &grid))
				status = -1;
			else
			{
				result = 2.0 * quantity[flatten(i1, j1, k1, p, max, grid)] -
					quantity[flatten(i2, j2, k2, p, max, grid)];
			}
		}
		else
			L_ERROR("Invalid order of extrapolation requested.", GridUtils::logfile);

		if (status != 0)
		{
			std::string msg;
			msg += "Extrapolation does not have enough available cells.";
			L_ERROR(msg, GridUtils::logfile);
		}

		return result;
	};


private:
	/// 2D/3D/4D flattener
	static int flatten(const int i, const int j, const int k, const int v, const int max, GridObj const & grid)
	{
		return (v + k * max + j * max * grid.K_lim + i * max * grid.K_lim * grid.M_lim);
	}

};

#endif