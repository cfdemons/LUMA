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

#pragma once
#include <vector>
#include <memory>

#include "InOutRepo.h"
#include "GridUtils.h"

class InOutData
{
public:
	InOutData();
	~InOutData();

	// Adds repo to the data
	int addRepo(std::unique_ptr<InOutRepo>&& repo);

	// Creates a new repo with the coordFile and the position of the LBM domain in world coordinates. 
	int addRepo(std::string coordFile, double x0, double y0, double z0);

	// Creates a new repo with coordinates_ with Nx, Ny and Nz cells separated dx and the position (x, y, z) of the LB domain in world coordinates
	int addRepo(int Nx, int Ny, int Nz, double x, double y, double z, double dx);

	// Set / unset all the repos to ready to write its data to the LBM simulation. 
	void setAllReadyToWriteToLBM( bool write);

	bool removeRepo(int repoNumber);
	
	InOutRepo& getRepo(int i) { return *data_[i]; }

	int getNumRepos() { return data_.size(); }

private:
	std::vector<std::unique_ptr<InOutRepo>> data_;
};
