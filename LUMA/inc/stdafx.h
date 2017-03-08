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

// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef STDAFX_H
#define STDAFX_H

#ifdef _DEBUG
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>	
	#ifndef DBG_NEW
		#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
		#define new DBG_NEW
	#endif
#endif  // _DEBUG

// Deprecated macros depending on platform
#ifdef __GNUC__
	#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
	#define DEPRECATED __declspec(deprecated)
#else
	#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
	#define DEPRECATED
#endif

// Frequently used headers (speeds up compilation in VS if put in the pre-compiled header module)
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <assert.h>

// Check OS is Windows or not
#ifdef _WIN32

#include <SDKDDKVer.h>
#define NOMINMAX		// Stop Windows.h redefining min/max
#include <Windows.h>
#include <tchar.h>

#else // Compiling with gcc through Code::Blocks / Eclipse on Linux

#include <stdlib.h> // Includes exit() function
#include <cstring>

#endif
#include <stdio.h>


/****************************************************/
// Enumerations //
/****************************************************/

#include "Enumerations.h"


/****************************************************/
// Our definitions //
/****************************************************/

#ifdef _WIN32
#define L_IS_NAN _isnan			///< Not a Number declaration (Windows)
#else
#define L_IS_NAN std::isnan		///< Not a Number declaration (Unix)
#endif

// Squared operator
#define SQ(x) ((x) * (x))

// Small number for comparing floats to zero
#define L_SMALL_NUMBER 1e-8

// Failure coede
#define LUMA_FAILED 12345


/****************************************************/
// Our headers (include after the enumerations) //
/****************************************************/

// Include definitions, singletons and headers to be made available everywhere for convenience.
#include "definitions.h"
#include "GridManager.h"
#include <mpi.h>
#include "MpiManager.h"
#include "GridUtils.h"
#include "GridUnits.h"


/****************************************************/
// Our global functions //
/****************************************************/

// Global variable references
extern const int c[3][L_NUM_VELS];				///< Lattice velocities
extern const int c_opt[L_NUM_VELS][3];			///< Lattice velocities optimised arrangement
extern const double w[L_NUM_VELS];				///< Quadrature weights
extern const double cs;							///< Lattice sound speed

// Debug stuff -- maybe I should put all these debug statements into some static
// class and just compile blank functions if not in debugging mode?
#ifdef L_IBM_DEBUG
#define L_DACTION_WRITE_OUT_FORCES \
std::ofstream testout; \
testout.open(GridUtils::path_str + "/force_i_LB.out", std::ios::app); \
testout << "\nNEW TIME STEP" << std::endl; \
for (size_t j = 1; j < M_lim - 1; j++) { \
	for (size_t i = 0; i < N_lim; i++) { \
		for (size_t v = 0; v < L_NUM_VELS; v++) { \
			testout << force_i(i, j, 0, v, M_lim, K_lim, L_NUM_VELS) << "\t"; \
		} \
		testout << std::endl; \
	} \
	testout << std::endl; \
} \
testout.close(); \

#else
	#define L_DACTION_WRITE_OUT_FORCES
#endif


/// Error definition
#define L_ERROR errorfcn	///< Error function shorthand
/// \brief Fatal Error function.
///
///			Writes error to the user and further information to the supplied logfile.
///			Inlined since this header is included everywhere.
///
///	\param	msg			string to be printed to the log file.
///	\param	logfile		pointer to the logfile where the message is to be written.
inline void errorfcn(const std::string &msg, std::ofstream *logfile)
{
#ifdef L_BUILD_FOR_MPI
	std::cout << "Rank " + std::to_string(MpiManager::getInstance()->my_rank);
#endif

	std::cout << " Error: See Log File" << std::endl;
	*logfile << "**ERROR**: " << msg << std::endl;
	logfile->close();

#ifdef L_BUILD_FOR_MPI
	MPI_Finalize();
#endif
	exit(LUMA_FAILED);
}

/// Regular writer
#define L_INFO infofcn	///< Info function shorthand
/// \brief Info / logger function.
///
///			Writes string to the supplied logfile.
///			Inlined since this header is included everywhere.
///
///	\param	msg			string to be printed to the log file.
///	\param	logfile		pointer to the logfile where the message is to be written.
inline void infofcn(const std::string &msg, std::ofstream *logfile)
{

	*logfile << "Info: " << msg << std::endl;
}

/// Regular writer
#define L_WARN warnfcn	///< Warning function shorthand
/// \brief Warning function.
///
///			Writes string to the supplied logfile.
///			Inlined since this header is included everywhere.
///
///	\param	msg			string to be printed to the log file.
///	\param	logfile		pointer to the logfile where the message is to be written.
inline void warnfcn(const std::string &msg, std::ofstream *logfile)
{

	*logfile << "WARNING: " << msg << std::endl;
}

#endif
