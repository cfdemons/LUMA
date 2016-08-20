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


#ifdef L_BUILD_FOR_MPI
	#define H5_HAVE_PARALLEL
#endif

#include "H5Cpp.h"

#define H5_BUILT_AS_DYNAMIC_LIB
#define HDF5_EXT_ZLIB
#define HDF5_EXT_SZIP

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif

