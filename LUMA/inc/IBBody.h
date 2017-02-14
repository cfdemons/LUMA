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
	IBBody(GridObj* g, size_t id);
	IBBody(GridObj* g, size_t id, PCpts* _PCpts);
	IBBody(GridObj* g, size_t bodyID, int lev, int reg, std::vector<double> &start_position,
		double length, std::vector<double> &angles, eMoveableType moveProperty, bool clamped);

protected:

	/************** Member Data **************/

	bool isFlexible;					///< Flag to indicate flexibility: false == rigid body; true == flexible filament
	bool isMovable;						///< Flag to indicate if body is movable or not.


	/************** Member Methods **************/

public:

	//////////////////////////////////
	// Prefab body building methods //
	//////////////////////////////////

	// Method to construct sphere/circle
	void makeBody(double radius, std::vector<double> centre, bool isFlexible, bool isMovable, int group);
	// Method to construct cuboid/rectangle
	void makeBody(std::vector<double> width_length_depth, std::vector<double> angles, std::vector<double> centre,
		bool isFlexible, bool isMovable, int group);
	// Method to construct filament
	void makeBody(int numbermarkers, std::vector<double> start_point, double fil_length, std::vector<double> angles, std::vector<int> BCs,
		bool isFlexible, bool isMovable, int group);
	// Method to construct a 3D plate
	double makeBody(std::vector<double> width_length, double angle, std::vector<double> centre,
		bool isFlexible, bool isMovable, int group, bool plate);

};

#endif
