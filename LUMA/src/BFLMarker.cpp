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

#include "../inc/stdafx.h"
#include "../inc/BFLMarker.h"
#include "../inc/GridUtils.h"

// Implementation for the BFLMarker class //


/// Default constructor
BFLMarker::BFLMarker(void)
{
}

/// Default destructor
BFLMarker::~BFLMarker(void)
{
}

/// \brief Custom constructor with position.
/// \param x x-position of marker
/// \param y y-position of marker
/// \param z z-position of marker
BFLMarker::BFLMarker(double x, double y, double z)
{
	// Add its position
	position.push_back(x);
	position.push_back(y);
	position.push_back(z);

	// Add indices (global reference frame)
	std::vector<int> vox = GridUtils::getVoxInd(x, y, z);
	this->supp_i.push_back(vox[0]);
	this->supp_j.push_back(vox[1]);
	this->supp_k.push_back(vox[2]);

}
