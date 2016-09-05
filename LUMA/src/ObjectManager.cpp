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

/*******************************************************************************/
// Lift and drag calculation (Mei's formula)
void ObjectManager::computeLiftDrag(int i, int j, int k, int M_lim, int K_lim, GridObj *g) {

	// TODO: Need abounding box for object if we have walls in the domain otherwise they will also be counted
	// TODO: Also need to be able to identify which body this site relates to so we can differentiate

	// For MPI builds, ignore if part of object is in halo region
#ifdef L_BUILD_FOR_MPI
	if (!GridUtils::isOnRecvLayer(g->XPos[i], g->YPos[j], g->ZPos[k]))
#endif
	{

		// Loop over directions from solid site
		for (int n = 0; n < L_nVels; n++) {

			int n_opp = GridUtils::getOpposite(n); // Get incoming direction

			// Compute destination coordinates
			int xdest = i + c[0][n];
			int ydest = j + c[1][n];
			int zdest = k + c[2][n];

			// Only apply if streams to a fluid site
			if (g->LatTyp(xdest, ydest, zdest, M_lim, K_lim) == eFluid)
			{

				force_on_object_x += c[0][n_opp] *
					(g->f(i, j, k, n, M_lim, K_lim, L_nVels) + g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_nVels));
				force_on_object_y += c[1][n_opp] *
					(g->f(i, j, k, n, M_lim, K_lim, L_nVels) + g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_nVels));
				force_on_object_z += c[2][n_opp] *
					(g->f(i, j, k, n, M_lim, K_lim, L_nVels) + g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_nVels));

			}
		}
	}
}