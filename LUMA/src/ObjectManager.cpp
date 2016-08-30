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
#include "../inc/ObjectManager.h"


// Static declarations
ObjectManager* ObjectManager::me;

// *************************************************************************************************** //
// Instance creator
ObjectManager* ObjectManager::getInstance() {

	if (!me) me = new ObjectManager;	// Private construction
	return me;							// Return pointer to new object

}

// Instance creator (only works the first time obj
ObjectManager* ObjectManager::getInstance(GridObj* g) {

	if (!me) me = new ObjectManager(g);	// Private construction
	return me;							// Return pointer to new object

}

// Instance destuctor
void ObjectManager::destroyInstance() {

	if (me)	delete me;			// Delete pointer from static context not destructor

}

// *************************************************************************************************** //
// Constructor & destructor
ObjectManager::ObjectManager(void) { };
ObjectManager::~ObjectManager(void) {
	me = nullptr;
};
ObjectManager::ObjectManager(GridObj* g) {
	this->_Grids = g;
};


// *************************************************************************************************** //
// Voxelisation Utilities

// Return global voxel indices for a given point in global space
std::vector<int> ObjectManager::getVoxInd(double x, double y, double z) {

	std::vector<int> vox;

	if (x - (int)std::floor(x) > 0.5) vox.push_back((int)std::ceil(x));
	else vox.push_back((int)std::floor(x));

	if (y - (int)std::floor(y) > 0.5) vox.push_back((int)std::ceil(y));
	else vox.push_back((int)std::floor(y));

	if (z - (int)std::floor(z) > 0.5) vox.push_back((int)std::ceil(z));
	else vox.push_back((int)std::floor(z));

	return vox;

}

// Overload of above for a single index
int ObjectManager::getVoxInd(double p) {

	if (p - (int)std::floor(p) > 0.5) return (int)std::ceil(p);
	else return (int)std::floor(p);

}