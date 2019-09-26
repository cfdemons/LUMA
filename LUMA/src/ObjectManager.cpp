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

	// Resize vector of flexible body flags
	hasIBMBodies.resize(L_NUM_LEVELS+1 ,false);
	hasFlexibleBodies.resize(L_NUM_LEVELS+1 ,false);

	// Set sub-iteration loop values
	timeav_subResidual = 0.0;
	timeav_subIterations = 0.0;
};

// ************************************************************************* //
/// \brief	Compute forces on a BB rigid object.
///
///			Uses momentum exchange to compute forces on rigid bodies.
///			Currently working with bounce-back objects only. There is no 
///			bounding box so if we have walls in the domain they will be counted 
///			as well.
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

	// For MPI builds, ignore if part of object is in halo region as represented on another rank
#ifdef L_BUILD_FOR_MPI
	if (!GridUtils::isOnRecvLayer(g->XPos[i], g->YPos[j], g->ZPos[k]))
#endif
	{

#ifdef L_MOMEX_DEBUG
		// Write position of solid site to debugging file
		if (debugstream.is_open())
			debugstream << std::endl << g->XPos[i] << "," << g->YPos[j] << "," << g->ZPos[k];
#endif
		// Loop over directions from solid site
		for (int n = 0; n < L_NUM_VELS; n++)
		{
			// Declare some local stores for this site
			double contrib_x = 0.0, contrib_y = 0.0, contrib_z = 0.0;

			// Get incoming direction
			int n_opp = GridUtils::getOpposite(n);

			// Compute destination coordinates (does not assume any periodicity)
			int xdest = i - c[eXDirection][n_opp];
			int ydest = j - c[eYDirection][n_opp];
			int zdest = k - c[eZDirection][n_opp];

			// Objects on edges will not get a periodic contribution && only applies to fluid sites
			if (!GridUtils::isOffGrid(xdest, ydest, zdest, g) &&
				g->LatTyp(xdest, ydest, zdest, M_lim, K_lim) == eFluid)
			{
				/* For HWBB:
				 *
				 *	Force =
				 *		(pre-stream population toward wall +
				 *		post-stream population away from wall)
				 *
				 * since population is simply bounced-back, we can write as:
				 *
				 *	Force =
				 *		(2 * pre-stream population toward wall)
				 *
				 * Multiplication by c unit vector resolves the result in
				 * appropriate direction.
				 */

				 // Store contribution in this direction
				contrib_x = 2.0 * c[eXDirection][n_opp] * g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_NUM_VELS);
				contrib_y = 2.0 * c[eYDirection][n_opp] * g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_NUM_VELS);
				contrib_z = 2.0 * c[eZDirection][n_opp] * g->f(xdest, ydest, zdest, n_opp, M_lim, K_lim, L_NUM_VELS);
			}

#ifdef L_MOMEX_DEBUG
			// Write contribution to file for this site for this lattice link
			if (debugstream.is_open())
				debugstream << "," << std::to_string(contrib_x) << "," << std::to_string(contrib_y) << "," << std::to_string(contrib_z);
#endif
			// Add the contribution of this link to the body forces
			bbbForceOnObjectX += contrib_x;
			bbbForceOnObjectY += contrib_y;
			bbbForceOnObjectZ += contrib_z;

		}
	}
}

// ************************************************************************* //
/// \brief	Compute forces on a BFL rigid object.
///
///			Uses momentum exchange to compute forces on a marker than makes up
///			a BFL body. Currently only works with a single BFL body but can 
///			easily be upgraded.
///
///	\param	v			lattice direction of link being considered.
///	\param	id			collapsed ijk index for site on which BFL BC is being applied.
/// \param	g			pointer to grid on which marker resides.
/// \param	markerID	id of marker on which force is to be updated.
void ObjectManager::computeLiftDrag(int v, int id, GridObj *g, int markerID)
{
	// Get opposite once
	int v_opp = GridUtils::getOpposite(v);

	// Similar to BBB but we cannot assume that bounced-back population is the same anymore
	pBody[0].markers[markerID].forceX +=
		c[eXDirection][v_opp] * (g->f[v_opp + id * L_NUM_VELS] + g->fNew[v + id * L_NUM_VELS]);
	pBody[0].markers[markerID].forceY +=
		c[eYDirection][v_opp] * (g->f[v_opp + id * L_NUM_VELS] + g->fNew[v + id * L_NUM_VELS]);
	pBody[0].markers[markerID].forceZ +=
		c[eZDirection][v_opp] * (g->f[v_opp + id * L_NUM_VELS] + g->fNew[v + id * L_NUM_VELS]);
}

// ************************************************************************* //
/// \brief	Resets the body force members prior to a new force calculation
///			using momentum exchange.
///
///	\param	grid	Grid object on which method was called
void ObjectManager::resetMomexBodyForces(GridObj * grid)
{
	if (grid->level == bbbOnGridLevel && grid->region_number == bbbOnGridReg)
	{
		bbbForceOnObjectX = 0.0;
		bbbForceOnObjectY = 0.0;
		bbbForceOnObjectZ = 0.0;

#ifdef L_MOMEX_DEBUG
		// Open file for momentum exchange information
		toggleDebugStream(grid);
#endif

	}

	// Reset the BFL body marker forces
	for (BFLBody& body : pBody)
	{
		// Only reset if body on this grid
		if (body._Owner->level == grid->level &&
			body._Owner->region_number == grid->region_number)
		{

			for (BFLMarker& marker : body.markers)
			{
				marker.forceX = 0.0;
				marker.forceY = 0.0;
				marker.forceZ = 0.0;
			}
		}
	}
}

// ************************************************************************* //
/// \brief	Adds a bounce-back body to the grid by labelling sites.
///
///			Override of the usual method which tries to place the object on the 
///			finest grid it can rather than a given grid. This will allow objects
///			to span multiple levels.
///
/// \param	geom	pointer to structure containing object information read from config file.
/// \param	_PCpts	pointer to point cloud information.
void ObjectManager::addBouncebackObject(GeomPacked *geom, PCpts *_PCpts)
{

	// Store information about the body in the Object Manager
	bbbOnGridLevel = geom->onGridLev;
	bbbOnGridReg = geom->onGridReg;

	// Declarations
	std::vector<int> ijk;
	eLocationOnRank loc = eNone;
	GridObj *g = nullptr;
	bool bPointAdded = false;
	eType localType;

	// Loop over the points
	for (int a = 0; a < static_cast<int>(_PCpts->x.size()); a++)
	{
		// Reset flag
		bPointAdded = false;

		// Loop over possible grids from bottom up
		for (int lev = L_NUM_LEVELS; lev >= 0; lev--)
		{
			for (int reg = 0; reg < L_NUM_REGIONS; reg++)
			{
				GridUtils::getGrid(lev, reg, g);

				// Skip if cannot find grid
				if (!g) continue;

				// If found grid then check in range
				if (GridUtils::isOnThisRank(_PCpts->x[a], _PCpts->y[a], _PCpts->z[a], &loc, g, &ijk))
				{

					localType = g->LatTyp(ijk[0], ijk[1], ijk[2], g->M_lim, g->K_lim);

					/* Update Typing Matrix and correct macroscopic.
					 * We must allow labelling on TL but recall that TL2C sites 
					 * which pull from a refined region may need to have BB applied
					 * so also label all the refined sites behind the fine grids 
					 * with the solid shape to make sure this is consistent. */
					if (localType != eVelocity)
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

				/* Do not break but try add the solid site on every grid behind 
				 * the finest grid -- see comment above as to why */
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
	bbbOnGridLevel = geom->onGridLev;
	bbbOnGridReg = geom->onGridReg;

	// Declarations
	std::vector<int> ijk;
	eLocationOnRank loc = eNone;

	// Label the grid sites
	for (int a = 0; a < static_cast<int>(_PCpts->x.size()); a++)
	{

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
/// Private method for opening/closing a debugging file
///	\param	g	pointer to grid toggling the stream
void ObjectManager::toggleDebugStream(GridObj *g)
{
	// Only do this if on correct time interval
	if (g->t == 0 || 
		(g->t + 1) % static_cast<int>(L_EXTRA_OUT_FREQ * (1.0 / g->refinement_ratio)) != 0) return;

	// Open file if not open, close if already open
	if (!debugstream.is_open())
	{
		debugstream.open(GridUtils::path_str + "/momex_debug_" + 
			std::to_string(static_cast<int>((g->t + 1) * g->refinement_ratio)) + 
			"_Rnk" + std::to_string(GridUtils::safeGetRank()) + ".csv", std::ios::out);

		// Add header for MomEx debug
		debugstream << "X Position,Y Position,Z Position";

		for (int v = 0; v < L_NUM_VELS; ++v)
		{
			debugstream << ",F" + std::to_string(v) + "X,F" + std::to_string(v) + "Y,F" + std::to_string(v) + "Z";
		}
	}
	else
	{
		debugstream.close();
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
	int onGridLev, int onGridReg,
	bool isCentreX, double refX,
	bool isCentreY, double refY,
	bool isCentreZ, double refZ,
	double bodyLength, eCartesianDirection scaleDirection,
	eMoveableType moveProperty, bool isClamped
	)
	: objtype(objtype), bodyID(bodyID), fileName(fileName),
	onGridLev(onGridLev), onGridReg(onGridReg),
	isRefXCentre(isCentreX), bodyRefX(refX),
	isRefYCentre(isCentreY), bodyRefY(refY),
	isRefZCentre(isCentreZ), bodyRefZ(refZ),
	bodyLength(bodyLength), scaleDirection(scaleDirection),
	moveProperty(moveProperty), isClamped(isClamped)
{
}

/// Method to interpret the reference type read in from the gerometry file
bool ObjectManager::GeomPacked::interpretRef(std::string refType)
{
	if (refType == "CENTRE")
		return true;
	else if (refType != "START")
		L_ERROR("Unknown reference type in geometry file. Exiting.", GridUtils::logfile);

	return false;
}