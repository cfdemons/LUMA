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

#include "../inc/InOutRepo.h"
//#include "Utils.h"
//#include "BoundingBox.h"

void InOutRepo::setCoordinates(std::string fileName, double x, double y, double z)
{
	if (numDataPoints_ > 0)
	{
		L_ERROR("This InOutRepo already contains a coordinates vector", GridUtils::logfile);
	}

	x0_ = x;
	y0_ = y;
	z0_ = z;

	// Open input file
	std::ifstream file;
	file.open(fileName.c_str());
	if (!file.is_open())
	{
		L_ERROR("File: " + fileName + " not found!", GridUtils::logfile);
	}

	// Read the file and fill a vector with the coordinates (I don't need to know how long the vectors will be)
	double component;
	while (file >> component)
	{
		coordinates_.push_back(component);
	}

	file.close();

	numDataPoints_ = coordinates_.size() / 3;

	L_INFO("Input data file size: " + std::to_string(coordinates_.size() / 3), GridUtils::logfile);
}

void InOutRepo::setCoordinates(InOutRepo* other)
{
	if (numDataPoints_ > 0)
	{
		L_ERROR("This InOutRepo already contains a coordinates vector", GridUtils::logfile);
	}

	if (other->getinitElements() < 1)
	{
		L_ERROR("The InOutRepo used to initialise the current InOutRepo contains no coordinates information. Aborting simulation.", GridUtils::logfile);
	}

	x0_ = other->getx0();
	y0_ = other->gety0();
	z0_ = other->getz0();

	numDataPoints_ = other->getinitElements();

	// Copy other's coordinates to coordinates_
	for (int i = 0; i < (3 * numDataPoints_); i++)
	{
		coordinates_.push_back(other->getCoordinates()[i]);
	}
}

void InOutRepo::setCoordinates(int Nx, int Ny, int Nz, double x, double y, double z, double dx)
{
	if (numDataPoints_ > 0)
	{
		L_ERROR("This InOutRepo already contains a coordinates vector", GridUtils::logfile);
	}

	x0_ = x;
	y0_ = y;
	z0_ = z;

	// Fill in the coordinates vector. C style ordering
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				coordinates_.push_back(x0_ + i * dx + dx / static_cast<double>(2.0));
				coordinates_.push_back(y0_ + j * dx + dx / static_cast<double>(2.0));
				coordinates_.push_back(z0_ + k * dx + dx / static_cast<double>(2.0));
			}
		}
	}

	numDataPoints_ = coordinates_.size() / 3;

	// Fill in the ordered data.
	numX_ = Nx;
	numY_ = Ny;
	numZ_ = Nz;
	C_ = true;
	ordered_ = true;
	dx_ = dx;
}

// Read inlet conditions from file and fill in the repo with the information
// It can read two types of files: 
// - SEM files: First line: SEM, Second line: U V W uu uv uw vv vw ww epsilon, Rest of lines: data in columns in the order of the second line. 
// - Velocity files: First line: velocity, Second line: u v w, Rest of lines: data in columns in the order of the second line.
void InOutRepo::readInletFromFile(std::string fname)
{
	if (numDataPoints_ < 1)
	{
		L_ERROR("This InOutData doesn't contain coordinates information. Please add coordinates before trying to add flow data", GridUtils::logfile);
	}


	L_INFO("Reading in inlet BC from file...", GridUtils::logfile);

	// Load the file
	std::ifstream file;
	file.open(fname, std::ios::in);
	if (!file.is_open())
	{
		L_ERROR("BC file cannot be opened.", GridUtils::logfile);
	}
	
	// Read the file header
	std::string line_in;	// String to store line in
	GridUtils::readLine(file, line_in);

	// If the file contains SEM data
	if (line_in.compare("SEM") == 0)
	{
		// Create space for the data in  
		initialiseVectorData("uMean");
		initialiseVectorData("ReDiag");
		initialiseVectorData("ReUpperDiag");
		initialiseScalarData("epsIn");

		// Check that the second line of the header contains the correct data. 
		std::string line2;
		GridUtils::readLine(file, line2);
		std::string header = "U V W uu uv uw vv vw ww epsilon";
		GridUtils::checkHeader(header, line2, fname);
	}
	else if (line_in.compare("velocity") == 0) // If the file contains velocity data only
	{
		initialiseVectorData("velocity");

		// Check that the second line of the header contains the correct data. 
		std::string line2;
		GridUtils::readLine(file, line2);
		std::string header = "u v w";
		GridUtils::checkHeader(header, line2, fname);
	}
	else
	{
		L_ERROR("Input file data type " + line_in + ", not recognised. Stopping simulation. ", GridUtils::logfile);
	}


	// Loop over the number of elements in InOutData
	for (int i = 0; i < getinitElements(); i++)
	{
		if (!file.eof())
		{
			// Read in one line of file at a time
			std::string line_in;	// String to store line in
			std::istringstream iss;	// Buffer stream to store characters

									// Get line up to new line separator and put in buffer
			GridUtils::readLine(file, line_in);
			iss.str(line_in);	// Put line in the buffer
			iss.seekg(0);		// Reset buffer position to start of buffer

			if (getVectorData("uMean") != nullptr)
			{
				iss >> getVectorData("uMean")[3 * i]
					>> getVectorData("uMean")[3 * i + 1]
					>> getVectorData("uMean")[3 * i + 2];
			}

			if (getVectorData("velocity") != nullptr)
			{
				iss >> getVectorData("velocity")[3 * i]
					>> getVectorData("velocity")[3 * i + 1]
					>> getVectorData("velocity")[3 * i + 2];
			}

			if ((getVectorData("ReDiag") != nullptr) && (getVectorData("ReUpperDiag") != nullptr))
			{
				iss >> getVectorData("ReDiag")[3 * i]
					>> getVectorData("ReUpperDiag")[3 * i]
					>> getVectorData("ReUpperDiag")[3 * i + 1]
					>> getVectorData("ReDiag")[3 * i + 1]
					>> getVectorData("ReUpperDiag")[3 * i + 2]
					>> getVectorData("ReDiag")[3 * i + 2];
			}

			if (getScalarData("epsIn") != nullptr)
			{
				iss >> getScalarData("epsIn")[i];
			}
		}
		else
		{
			L_ERROR("File " + fname + ", contains less elements than the capacity of the data to fill. Data capacity: "
				+ std::to_string(getinitElements()) + ". Stopping simulation. ", GridUtils::logfile);
		}
	}
}

bool InOutRepo::addDummyData(const std::vector<std::string>& variables)
{
	for (size_t v = 0; v < variables.size(); v++)
	{
		if (variables[v].compare("uMean") == 0)
		{
			// Create the variable
			initialiseVectorData(variables[v], 0.0);

			double a = coordinates_[2];
			double b = coordinates_[3*(numDataPoints_ -1) + 1];
			double p = static_cast<T>(0.5) * (a + b);
			double q = b - p;

			int ind = 0;
			for (int i = 0; i < numX_; i++)
			{
				for (int j = 0; j < numY_; j++)
				{
					
					for (int k = 0; k < numZ_; k++)
					{
						double yPos = coordinates_[3 * ind + 1];

						getVectorData(variables[v])[3*ind] = 1.5 * (1.0 - std::pow((yPos - p) / q, 2));
						ind++;
					}
				}
			}
		}
		else if (variables[v].compare("ReDiag") == 0)
		{
			// Create the variable
			initialiseVectorData(variables[v], 0.0);

			double a = coordinates_[2];
			double b = coordinates_[3 * (numDataPoints_ - 1) + 1];
			double p = 0.5 * (a + b);
			double q = b - p;

			int ind = 0;
			for (int i = 0; i < numX_; i++)
			{
				for (int j = 0; j < numY_; j++)
				{
					for (int k = 0; k < numZ_; k++)
					{
						double yPos = coordinates_[3 * ind + 1];
						getVectorData(variables[v])[3 * ind] =     1.5 - (1.5 * (1.0 - std::pow((yPos - p) / q, 2))) / 100.0;
						getVectorData(variables[v])[3 * ind + 1] = 1.5 - (1.5 * (1.0 - std::pow((yPos - p) / q, 2))) / 1000.0;
						getVectorData(variables[v])[3 * ind + 2] = 1.5 - (1.5 * (1.0 - std::pow((yPos - p) / q, 2))) / 2000.0;
						ind++;
					}
				}
			}
		}
		else if (variables[v].compare("ReUpperDiag") == 0)
		{
			// Create the variable
			initialiseVectorData(variables[v], 0.0);

			int ind = 0;
			for (int i = 0; i < numX_; i++)
			{
				for (int j = 0; j < numY_; j++)
				{
					for (int k = 0; k < numZ_; k++)
					{
						double yPos = coordinates_[3 * ind + 1];
						getVectorData(variables[v])[3 * ind] = 0.003175 * yPos - 0.003317;		
						ind++;
					}
				}
			}
		}
		else if (variables[v].compare("epsIn") == 0)
		{
			// Create the variable
			initialiseScalarData(variables[v], 0.0);

			double a = coordinates_[2];
			double b = coordinates_[3 * (numDataPoints_ - 1) + 1];
			double p = 0.5 * (a + b);
			double q = b - p;

			int ind = 0;
			for (int i = 0; i < numX_; i++)
			{
				for (int j = 0; j < numY_; j++)
				{
					for (int k = 0; k < numZ_; k++)
					{
						double yPos = coordinates_[3 * ind + 1];
						getScalarData(variables[v])[ind] = 1.0 - (1.0 * (1.0 - std::pow((yPos - p) / q, 2))) / 100.0;
					}
				}
			}
		}
		else
		{
			L_INFO("Variable " + variables[v] + " not recognised. Cannot create dummy data for it.", GridUtils::logfile);
			return false;
		}
	}
	
	return true;


}

bool InOutRepo::ordered(int numX, int numY, int numZ, double dx, bool C)
{
	if (numDataPoints_ < 1)
	{
		L_ERROR("This InOutRepo does not contain coordinates. Aborting simulation", GridUtils::logfile);
	}
	
	if ((numX * numY * numZ) != numDataPoints_)
	{
		return false;
	}

	// I have to implement the part that actually checks if the data is ordered!!!!!!!!
	/*double coordX = initCoord_[0];
	double coordY = initCoord_[1];
	double coordZ = initCoord_[2];

	if (C)
	{
		for (int i = 0; i < numX; i++)
		{
			for (int j = 0; j < numY; j++)
			{
				for (int k = 0; k < numZ; k++)
				{
					if (coordX - initCoord_[3 * ])

				}
			}
		}
	} */

	numX_ = numX;
	numY_ = numY;
	numZ_ = numZ;
	C_ = C;
	ordered_ = true;
	dx_ = dx;

	return true;
}

bool InOutRepo::isSEM()
{
	if (vectorData_.find("uMean") == vectorData_.end())
	{
		return false;
	}
	else if (vectorData_.find("ReDiag") == vectorData_.end())
	{
		return false;
	}
	else if (vectorData_.find("ReUpperDiag") == vectorData_.end())
	{
		return false;
	}
	else if (scalarData_.find("epsIn") == scalarData_.end())
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool InOutRepo::getIndices(double x, double y, double z, double *i, double *j, double *k)
{
	if (!ordered_)
	{
		return false;
	}

	double minx = coordinates_[0];
	double miny = coordinates_[1];
	double minz = coordinates_[2];

	// Check that the point x,y,z is inside of initCoord_
	/*if ((x < minx) || (x > maxx) || (y < miny) || (y > maxy) || (z < minz) || (z > maxz))
	{
		return false;
	}*/

	*i = (x - minx) / dx_;
	*j = (y - miny) / dx_;
	*k = (z - minz) / dx_;

	return true;
}

bool InOutRepo::getElementNum(int i, int j, int k, int * elementNum)
{
	if (!ordered_)
	{
		return false;
	}

	if (C_)
	{
		*elementNum = k + j * numZ_ + i * numZ_ * numY_;
	}
	else
	{
		*elementNum = i + j * numX_ + k * numX_ * numY_;
	}

	return true;
}

double* InOutRepo::getVectorData(std::string name)
{
	if (vectorData_.find(name) != vectorData_.end())
	{
		return vectorData_[name].data();
	}
	else
	{
		return nullptr;
	}
}

double * InOutRepo::getScalarData(std::string name)
{
	if (scalarData_.find(name) != scalarData_.end())
	{
		return scalarData_[name].data();
	}
	else
	{
		return nullptr;
	}
}


void InOutRepo::setWorldOffset(double x0, double y0, double z0)
{
	x0_ = x0;
	y0_ = y0;
	z0_ = z0;
}

void InOutRepo::initialiseVectorData(std::string name, double initValue)
{
	if (numDataPoints_ < 1)
	{
		L_ERROR("This InOutRepo doesn't contain coordinates. Please initialise the coordinates vector before initialising flow data", GridUtils::logfile);
	}

	// Free the previous arrays.
	if (vectorData_.find(name) != vectorData_.end())
	{
		vectorData_.erase(name);
	}

	// Allocate new arrays
	std::vector<double> temp;
	temp.resize(3 * numDataPoints_);
	vectorData_.emplace(name, temp);

	// Set all the data to the init value
	for (int i = 0; i < numDataPoints_; i++)
	{
		vectorData_[name].at(3 * i) = initValue; 
		vectorData_[name].at(3 * i + 1) = initValue;
		vectorData_[name].at(3 * i + 2) = initValue;
	}
}

void InOutRepo::initialiseScalarData(std::string name, double initValue)
{
	if (numDataPoints_ < 1)
	{
		L_ERROR("This InOutRepo doesn't contain coordinates. Please initialise the coordinates vector before initialising flow data", GridUtils::logfile);
	}

	// Free the previous arrays.
	if (scalarData_.find(name) != scalarData_.end())
	{
		scalarData_.erase(name);
	}

	// Allocate new arrays
	std::vector<double> temp;
	temp.resize(numDataPoints_);
	scalarData_.emplace(name, temp);

	// Set all the data to initValue
	for (int i = 0; i < numDataPoints_; i++)
	{
		scalarData_[name].at(i) = initValue;
	}
}

InOutRepo::InOutRepo()
{
	numDataPoints_ = -1;
}

InOutRepo::InOutRepo(std::string fileName, double x, double y, double z)
{
	setCoordinates(fileName, x, y, z);
}

InOutRepo::InOutRepo(int Nx, int Ny, int Nz, double x, double y, double z, double dx)
{
	setCoordinates(Nx, Ny, Nz, x, y, z, dx);
}

InOutRepo::~InOutRepo()
{


}

