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
#ifndef IBBODY_H
#define IBBODY_H

#include "stdafx.h"

// Forward declarations
#include "Body.h"	// This is a templated class so include the whole file
class IBMarker;
class PCpts;
class GridObj;

/// \brief	Immersed boundary body.
class IBBody : public Body<IBMarker> {

	// Add friend classes so they can access the protected data of IBBody objects
	friend class ObjectManager;
	friend class IBInfo;

public:
	// Constructor and destructor
	IBBody(void);
	~IBBody(void);

	// Custom constructor which takes pointer to point cloud data and a pointer to the grid owner for the labelling
	IBBody(GridObj* g, int bodyID, PCpts* _PCpts, eMoveableType moveProperty, bool clamped);

	// Custom constructor for building prefab circle or sphere
	IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point, double radius, eMoveableType moveProperty);

	// Custom constructor for building prefab square or cuboid
	IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		std::vector<double> &dimensions, std::vector<double> &angles, eMoveableType moveProperty);

	// Custom constructor for building prefab filament
	IBBody(GridObj* g, int bodyID, std::vector<double> &start_position,
		double length, std::vector<double> &angles, eMoveableType moveProperty, bool clamped);

protected:

	/************** Member Data **************/

	bool isFlexible;					///< Flag to indicate flexibility: false == rigid body; true == flexible filament
	bool isMovable;						///< Flag to indicate if body is movable or not.


	/************** Member Methods **************/

public:

};

#endif
