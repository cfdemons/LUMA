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
#pragma once
#include "stdafx.h"

// Class representing Point Cloud data
class PCpts {

public:

	PCpts(void) {};
	~PCpts(void) {};

	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
};