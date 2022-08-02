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

/* This file holds all the code for the core LBM operations including collision,
streaming and macroscopic calulcation.
*/

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/IVector.h"
#include "../inc/ObjectManager.h"

using namespace std;
// *****************************************************************************
/// \brief	Add an external data array to the LUMA flow fields.
///
///			Overwrites a LUMA macroscopic field at set locations with the values provided and recalculates the populations
///
/// \param name		Name of the field. 
/// \param cellIDs	Cell indices at which to overwrite the data. The indices are local to the current grid and process.
/// \param data		Data to write in the cellIDs indices. 
// TODO: add an option to interpolate the data.
void GridObj::coupling_addData(std::string name, const std::vector<int>& cellIDs, const std::vector<double>& data)
{
	// Chek which data to add. 
	if ((name.find("v") != std::string::npos) || (name.find("V") != std::string::npos))
	{
		// Check that cellIDs and data are the same size and that the size is > 0
		if ((cellIDs.size() == 0) || (cellIDs.size() != data.size()/L_DIMS))
		{
			L_ERROR("LUMA tried to read in external data but the data vector and cellID vector sizes do not match. Data name " +
				name + " cell ID size: " + std::to_string(cellIDs.size()) +
				"data size: " + std::to_string(data.size()), GridUtils::logfile);
		}
		std::cout << "This is id from addData" << std::endl;
		for (int i = 0; i < cellIDs.size(); i++)
		{
			u[0 + cellIDs[i] * L_DIMS] = data[0 + i * L_DIMS];
			u[1 + cellIDs[i] * L_DIMS] = data[1 + i * L_DIMS];
#if (L_DIMS == 3)
			u[2 + cellIDs[i] * L_DIMS] = data[2 + i * L_DIMS];
#endif
			for (int v = 0; v < L_NUM_VELS; ++v)
			{
				f[v + cellIDs[i] * L_NUM_VELS] =
					_LBM_equilibrium_opt(cellIDs[i], v);
			}
		}
	} 
	else if ((name.find("rho") != std::string::npos) || (name.find("RHO") != std::string::npos))
	{
		// Check that cellIDs and data are the same size and that the size is > 0
		if ((cellIDs.size() == 0) || (cellIDs.size() != data.size()))
		{
			L_ERROR("LUMA tried to read in external data but the data vector and cellID vector sizes do not match. Data name " +
				name + " cell ID size: " + std::to_string(cellIDs.size()) +
				"data size: " + std::to_string(data.size()), GridUtils::logfile);
		}

		for (int i = 0; i < cellIDs.size(); i++)
		{
			rho[cellIDs[i]] = data[i];
		}
	}
	else if ((name.find("temperature") != std::string::npos) || (name.find("TEMPERATURE") != std::string::npos))
	{
		for (int i = 0; i < cellIDs.size(); i++)
		{
			T[cellIDs[i]] = data[i];
		}
	}
	else
	{
		L_ERROR(" LUMA tried to read in external data but LUMA doesn't contain a solved variable called " + name, GridUtils::logfile);
	}
}

// *****************************************************************************
/// \brief	Extract data from LUMA and add it to an external data array
///
///			Writes the LUMA macroscopic field at set locations into the provided data array
///
/// \param name		Name of the field. 
/// \param cellIDs	Cell indices at which to extract the field data. The indices are local to the current grid and process.
/// \param data		Vector in which to store the data. 
// TODO: add an option to interpolate the data.
void GridObj::coupling_extractData(std::string name, const std::vector<int>& cellIDs, const std::vector<int> cellIDs_diff, std::vector<double>& data)
{
	// Chek which data to add. 
	if ((name.find("v") != std::string::npos) || (name.find("V") != std::string::npos))
	{
		// Resise the vector data to the size of cellIDs (vector data)
		data.resize(cellIDs.size() * L_DIMS);
 		for (int i = 0; i < cellIDs.size(); i++)
		{
			data[0 + i * L_DIMS] = (u[0 + cellIDs[i] * L_DIMS]+u[0 + cellIDs_diff[i] * L_DIMS])/2.0; //u[0 + cellIDs[i] * L_DIMS]; //
			data[1 + i * L_DIMS] = (u[1 + cellIDs[i] * L_DIMS]+u[1 + cellIDs_diff[i] * L_DIMS])/2.0; //u[1 + cellIDs[i] * L_DIMS]; //
#if (L_DIMS == 3)
			data[2 + i * L_DIMS] = (u[2 + cellIDs[i] * L_DIMS]+u[2 + cellIDs_diff[i] * L_DIMS])/2.0; //u[2 + cellIDs[i] * L_DIMS]; //
#endif
			//std::cout << cellIDs[i] << " " << cellIDs_diff[i]<<" ";
		} 
	}
	else if ((name.find("rho") != std::string::npos) || (name.find("RHO") != std::string::npos))
	{
		// Resise the vector data to the size of cellIDs
		data.resize(cellIDs.size());
		for (int i = 0; i < cellIDs.size(); i++)
		{
			data[i] = rho[cellIDs[i]];
		}
	}
	else if((name.find("temperature") != std::string::npos) || (name.find("TEMPERATURE") != std::string::npos))
	{
		
		// Resise the scalar data to the size of cellIDs
		data.resize(cellIDs.size());
		for (int i = 0; i < cellIDs.size(); i++)
		{
			data[i] = T[cellIDs[i]];
		}
	
	}
	else
	{
		L_ERROR(" LUMA tried to read in external data but LUMA doesn't contain a solved variable called " + name, GridUtils::logfile);
	}
}

