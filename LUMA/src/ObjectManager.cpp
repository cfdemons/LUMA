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
#include "../inc/GridObj.h"
#include "../inc/GridUtils.h"


// Static declarations
ObjectManager* ObjectManager::me;

// ************************************************************************* //
/// Instance creator
ObjectManager* ObjectManager::getInstance() {

	if (!me) me = new ObjectManager;	// Private construction
	return me;							// Return pointer to new object

}

/// \brief Instance creator with grid hierarchy assignment.
/// \param	g	pointer to grid hierarchy.
ObjectManager* ObjectManager::getInstance(GridObj* g) {

	if (!me) me = new ObjectManager(g);	// Private construction
	return me;							// Return pointer to new object

}

/// Instance destuctor
void ObjectManager::destroyInstance() {

	if (me)	delete me;			// Delete pointer from static context not destructor

}

// ************************************************************************* //
/// Default constructor
ObjectManager::ObjectManager(void) { };

/// Default destructor
ObjectManager::~ObjectManager(void) {
	me = nullptr;
};

/// \brief Constructor with grid hierarchy assignment.
/// \param	g	pointer to grid hierarchy.
ObjectManager::ObjectManager(GridObj* g) {
	this->_Grids = g;
};

// ************************************************************************* //

/// \brief	Compute forces on a rigid object.
///
///			Uses momentum exchange to compute forces on rigid bodies.
///			Currently working with bounce-back objects only. There is no 
///			bounding box so if we have walls in the domain they will be counted 
///			as well. Also only possible to differentiate between bodies. Lumps 
///			all bodies together.identify which body this site relates to so 
///			we can differentiate.
///
/// \param	i	local i-index of solid site.
/// \param	j	local j-index of solid site.
/// \param	k	local k-index of solid site.
/// \param	g	pointer to grid on which object resides.
void ObjectManager::computeLiftDrag(int i, int j, int k, GridObj *g) {

	// TODO: Need abounding box for object if we have walls in the domain otherwise they will also be counted
	// TODO: Also need to be able to identify which body this site relates to so we can differentiate

	int N_lim = g->N_lim;
	int M_lim = g->M_lim;
	int K_lim = g->K_lim;

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