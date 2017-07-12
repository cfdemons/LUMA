/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) The University of Manchester 2017
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * further distribution commericially or otherwise without written consent.
 *
 */
#ifndef IBBODY_H
#define IBBODY_H

#include "stdafx.h"

// Forward declarations
#include "Body.h"		// This is a templated class so include the whole file
#include "FEMBody.h"
class IBMarker;
class PCpts;
class GridObj;
class FEMBody;

/// \brief	Immersed boundary body.
class IBBody : public Body<IBMarker> {

	// Add friend classes so they can access the protected data of IBBody objects
	friend class ObjectManager;
	friend class IBInfo;
	friend class FEMBody;
	friend class MpiManager;

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
		double length, double height, double depth, std::vector<double> &angles, eMoveableType moveProperty, int nElement, bool clamped, double density, double E);

	// Custom constructor for building dummy iBody for epsilon calculation
	IBBody(IBBody &iBody, std::vector<std::vector<double>> &recvBuffer);


protected:

	/************** Member Data **************/

	bool isFlexible;					///< Flag to indicate flexibility: false == rigid body; true == flexible filament
	bool isMovable;						///< Flag to indicate if body is movable or not.

	FEMBody *fBody;						///< Pointer to FEM body object


	/************** Member Methods **************/

private:

	void initialise(eMoveableType moveProperty);		// Initialisation wrapper for setting flags
	void getValidMarkers();								// Get valid markers in iBody (only relevant for owning rank)

};

#endif
