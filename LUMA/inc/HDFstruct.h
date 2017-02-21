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
#ifndef HDFSTRUCT_H
#define HDFSTRUCT_H

#include "stdafx.h"
/// \struct HDFstruct
/// \brief	Structure for storing halo information for HDF5.
///
///			Structure also stores the amount of writable data on the grid.
struct HDFstruct {

	int i_start;	///< Starting i-index for writable region
	int i_end;		///< Ending i-index for writable region
	int j_start;	///< Starting j-index for writable region
	int j_end;		///< Ending j-index for writable region
	int k_start;	///< Starting k-index for writable region
	int k_end;		///< Ending k-index for writable region

	// Identifiers
	int level;		///< Grid level to which these data correspond
	int region;		///< Region number to which these data correspond

	/// Writable data count
	unsigned int writable_data_count = 0;
};

#endif