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

protected:

	/************** Member Data **************/

	bool isFlexible;					///< Flag to indicate flexibility: false == rigid body; true == flexible filament
	bool isMovable;						///< Flag to indicate if body is movable or not.
	int groupID;						///< ID of IBbody group -- position updates can be driven from a flexible body in a group

	// Flexible body properties
	double delta_rho;					///< Difference in density between fluid and solid in lattice units
	double flexural_rigidity;			///< Young's modulus E * Second moment of area I
	std::vector<double> tension;		///< Tension between the current marker and its neighbour
	std::vector<int> BCs;				///< BCs type flags (flexible bodies)


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

	// Overload base class function for adding markers with flexible flag already specified
	void addMarker(double x, double y, double z, bool isFlexible);

};

#endif
