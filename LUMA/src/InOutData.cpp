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

#include "../inc/InOutData.h"
#include "stddef.h"
#include <memory>

InOutData::InOutData()
{
}

InOutData::~InOutData()
{
}

int InOutData::addRepo(std::unique_ptr<InOutRepo>&& repo)
{
	data_.push_back(std::move(repo));
	return static_cast<int>(data_.size()) - 1;
}

int InOutData::addRepo(std::string coordFile, double x, double y, double z)
{
	data_.emplace_back(std::make_unique<InOutRepo>(coordFile, x, y, z));
	return static_cast<int>(data_.size()) - 1;
}

int InOutData::addRepo(int Nx, int Ny, int Nz, double x, double y, double z, double dx)
{
	data_.emplace_back(std::make_unique<InOutRepo>(Nx, Ny, Nz, x, y, z, dx));
	return static_cast<int>(data_.size()) - 1;
}

void InOutData::setAllReadyToWriteToLBM(bool write)
{
	if (data_.size() != 0)
	{
		for (int i = 0; i < data_.size(); i++)
			data_.at(i)->bReadyToWriteToLBM = write;
	}
}

bool InOutData::removeRepo(int repoNumber)
{
	if ((size_t) repoNumber >= data_.size()) return false;
	data_.erase(data_.begin() + repoNumber);
	return true;
}

