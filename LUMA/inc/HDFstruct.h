/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
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