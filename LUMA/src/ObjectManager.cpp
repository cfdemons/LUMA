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
ObjectManager::ObjectManager(void) {
	_Grids = nullptr;
};

/// Default destructor
ObjectManager::~ObjectManager(void) {
	me = nullptr;
};

/// \brief Constructor with grid hierarchy assignment.
/// \param	g	pointer to grid hierarchy.
ObjectManager::ObjectManager(GridObj* g) : _Grids(g)
{
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

	// TODO: Need a bounding box for object if we have walls in the domain otherwise they will also be counted
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
		for (int n = 0; n < L_NUM_VELS; n++) {

			// Get incoming direction
			int n_opp = GridUtils::getOpposite(n);

			// Compute destination coordinates
			int xdest = i + c[0][n];
			int ydest = j + c[1][n];
			int zdest = k + c[2][n];

			// Reject site on grid edges (like single-cell walls)
			if (GridUtils::isOffGrid(xdest, ydest, zdest, g)) return;

			// Only apply if streams to a fluid site
			if (g->LatTyp(xdest, ydest, zdest, M_lim, K_lim) == eFluid)
			{

				bbbForceOnObjectX +=
					2.0 * c[0][n_opp] * g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_NUM_VELS);
				bbbForceOnObjectY += 
					2.0 * c[1][n_opp] * g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_NUM_VELS);
				bbbForceOnObjectZ +=
					2.0 * c[2][n_opp] * g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_NUM_VELS);

			}
		}
	}
}

// ************************************************************************* //
/// \brief	Adds a bounce-back body to the grid by labelling sites.
/// \param	g		pointer to grid on which object resides.
/// \param	geom	pointer to structure containing object information read from config file.
/// \param	_PCpts	pointer to point cloud information.
void ObjectManager::addBouncebackObject(GridObj *g, GeomPacked *geom, PCpts *_PCpts)
{
	// Store information about the body in the Object Manager
	bbbOnGridLevel = geom->on_grid_lev;
	bbbOnGridReg = geom->on_grid_reg;

	// Declarations
	std::vector<int> ijk;
	eLocationOnRank loc = eNone;

	// Label the grid sites
	for (int a = 0; a < static_cast<int>(_PCpts->x.size()); a++) {

		// Get indices if on this rank
		if (GridUtils::isOnThisRank(_PCpts->x[a], _PCpts->y[a], _PCpts->z[a], &loc, g, &ijk))
		{
			// Update Typing Matrix and correct macroscopic
			if (g->LatTyp(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim) == eFluid)
			{
				// Change type
				g->LatTyp(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim) = eSolid;

				// Change macro
				g->u(ijk[0], ijk[1], ijk[2], 0, g->M_lim, g->K_lim, L_DIMS) = 0.0;
				g->u(ijk[0], ijk[1], ijk[2], 1, g->M_lim, g->K_lim, L_DIMS) = 0.0;
#if (L_DIMS == 3)
				g->u(ijk[0], ijk[1], ijk[2], 0, g->M_lim, g->K_lim, L_DIMS) = 0.0;
#endif
				g->rho(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim) = L_RHOIN;

			}
		}
	}
}


// ************************************************************************* //
/// Geometry data structure container constructor.
ObjectManager::GeomPacked::GeomPacked()
{
}

/// Geometry data structure container destructor.
ObjectManager::GeomPacked::~GeomPacked()
{
}

/// Geometry data structure container custom constructor.
ObjectManager::GeomPacked::GeomPacked(
	eObjectType objtype, int bodyID, std::string fileName,
	int on_grid_lev, int on_grid_reg,
	double body_start_x, double body_start_y, double body_centre_z,
	double body_length, eCartesianDirection scale_direction,
	eMoveableType moveProperty, bool clamped
	)
	: objtype(objtype), bodyID(bodyID), fileName(fileName),
	on_grid_lev(on_grid_lev), on_grid_reg(on_grid_reg),
	body_start_x(body_start_x), body_start_y(body_start_y), body_centre_z(body_centre_z), 
	body_length(body_length), scale_direction(scale_direction),
	moveProperty(moveProperty), clamped(clamped)
{
}