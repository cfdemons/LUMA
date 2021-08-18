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
#include "stdafx.h"

// Class represents a repository for data which is to be communicated into and out of the underlying LBM simulation
class InOutRepo
{

public:
	InOutRepo();
	InOutRepo(std::string fileName, double x, double y, double z);
	InOutRepo(int Nx, int Ny, int Nz, double x, double y, double z, double dx);
	~InOutRepo();

	// Flags to indicate whether action needs to be taken on the data it holds
	bool bReadyToWriteToLBM = false;
	bool bReadyToReadFromLBM = false;

	// Set the coordinates
	void setCoordinates(std::string fileName, double x, double y, double z);
	void setCoordinates(InOutRepo* other);
	void setCoordinates(int Nx, int Ny, int Nz, double x, double y, double z, double dx);

	// Read inlet conditions from file and fill in the repo with the information
	// It can read two types of files: 
	// - SEM files: First line: SEM, Second line: U V W uu uv uw vv vw ww epsilon, Rest of lines: data in columns in the order of the second line. 
	// - Velocity files: First line: velocity, Second line: u v w, Rest of lines: data in columns in the order of the second line.
	void readInletFromFile(std::string fname);

	// Adds dummy data for each variable in the variables string. Returns false if the variable is not recognised
	bool addDummyData(const std::vector<std::string>& variables);

	// Checks if the data is ordered in an equispaced grid with numX, numY, numZ cells in x,y,z directions separated dx in each direction
    // C is true if the order of variation is x,y,z, false if it is z,y,x
	// TODO: This function is not finished. For the moment it only checks if the number of elements in the repo is = numX * numY * numZ
	// and sets numX_, numY_, numZ.
	bool ordered(int numX, int numY, int numZ, double dx, bool C);

	// Checks if the repo contains the necessary information to be used as a mean flow for SEM. 
	bool isSEM();

	// Get the indices from a coordinate. 
	// Returns false if the data is not ordered. 
	// It doesn't check if the point is inside the coordinates_ range. 
	bool getIndices(double x, double y, double z, double* i, double* j, double* k);

	// Returns the element number for a i,j,k index
	// Returns false if the data is not ordered. 
	// It doesn't check if the point is inside the initCoord_ range. 
	bool getElementNum(int i, int j, int k, int *elementNum);

	const double* getCoordinates() const { return coordinates_.data(); }
	double* getVectorData(std::string name);
	double* getScalarData(std::string name);
	int getinitElements() { return numDataPoints_; }
	int getNumX() { return numX_; }
	int getNumY() { return numY_; }
	int getNumZ() { return numZ_; }
	double getx0() { return x0_; }
	double gety0() { return y0_; }
	double getz0() { return z0_; }
	
	void setWorldOffset(double x0, double y0, double z0);

	void initialiseVectorData(std::string name, double initValue = 0.0);
	void initialiseScalarData(std::string name, double initValue = 0.0);

private:

	std::vector<double> coordinates_;   // Coordinates of the data. Format (x0,y0,z0,x1,y1,z1...,xn,yn,zn)
	std::map<std::string, std::vector<double>> vectorData_;
	std::map<std::string, std::vector<double>> scalarData_;
	int     numDataPoints_ = -1;

	// numX_, numY_, numZ_, dx_ are only filled in if the ordered method returns true. 
	bool ordered_ = false;
	bool C_ = false;
	int numX_ = -1;
	int numY_ = -1;
	int numZ_ = -1;
	double dx_;

	double x0_ = 0; 
	double y0_ = 0;
	double z0_ = 0;

};