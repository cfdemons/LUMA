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

#pragma once
#include "IBMarker.h"
#include "Body.h"
#include "PCpts.h"
#include "BFLBody.h"	// Need this because we use the Voxelising routine in makeBody()

/** Represents an IB body **/
class IBBody : public Body<IBMarker> {

	// Make ObjectManager a friend class so it can access the protected data of IBBody objects
	friend class ObjectManager;

public:
	// Constructor and destructor
	IBBody(void);
	~IBBody(void);
	IBBody(GridObj* g);

protected:

	/*
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/

	bool flex_rigid;					// Set flag for flexibility: false == rigid body; true == flexible filament
	bool deformable;					// Set flag for deformable body: false == rigid; true == deformable
	int groupID;				// ID of IBbody group -- position updates can be driven from a flexible body in a group

	// Flexible body properties
	double delta_rho;					// Difference in density between fluid and solid in lattice units
	double flexural_rigidity;			// Young's modulus E * Second moment of area I
	std::vector<double> tension;		// Tension between the current marker and its neighbour
	std::vector<int> BCs;				// BCs type flags (flexible bodies)


	/*
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

public:

	// Add marker to the body (overload parent addMarker method)
	void addMarker(double x, double y, double z, bool flex_rigid);

	//////////////////////////////////
	// Prefab body building methods //
	//////////////////////////////////

	// Method to construct sphere/circle
	void makeBody(double radius, std::vector<double> centre, bool flex_rigid, bool moving, int group);
	// Method to construct cuboid/rectangle
	void makeBody(std::vector<double> width_length_depth, std::vector<double> angles, std::vector<double> centre,
		bool flex_rigid, bool deform, int group);
	// Method to construct filament
	void makeBody(int numbermarkers, std::vector<double> start_point, double fil_length, std::vector<double> angles, std::vector<int> BCs,
		bool flex_rigid, bool deform, int group);
	// Method to construct a 3D plate
	double makeBody(std::vector<double> width_length, double angle, std::vector<double> centre,
		bool flex_rigid, bool deform, int group, bool plate);
	// Method to construct body from points read in from file
	void makeBody(PCpts* _PCpts);

};

