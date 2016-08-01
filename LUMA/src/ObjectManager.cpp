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



