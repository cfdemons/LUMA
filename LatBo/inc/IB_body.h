#pragma once

#include <vector>
#include "IB_marker.h"
#include "GridUtils.h"

class IB_body {

	// Make GridObj a friend class so it can access the protected data of IB_body objects
	friend class GridObj;

public:
	// Constructor and destructor
	IB_body();
	IB_body(GridUtils&);	// Custom constructor which passes utility class
	~IB_body(void);

protected:

	/*
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/

	std::vector<IB_marker> markers;		// Array of Lagrange markers which make up the body
	double spacing;						// Spacing of the Lagrange markers in physical units
	bool flex_rigid;					// Set flag for flexibility: false == rigid body; true == flexible filament
	bool deformable;					// Set flag for deformable body: false == rigid; true == deformable
	unsigned int groupID;				// ID of IBbody group -- position updates can be driven from a flexible body in a group
	GridUtils gUtils;					// Reference to parent GridObj's utility class
	bool closed_surface;                // Flag to specify whether or not it is a closed surface (for output)

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
	// Add Lagrange marker to the body
	void addMarker(double x, double y, double z, bool flex_rigid);


	//////////////////////////////////
	// Prefab body building methods //
	//////////////////////////////////

	// Method to construct sphere/circle
	void makeBody(double radius, std::vector<double> centre, bool flex_rigid, bool moving, unsigned int group);
	// Method to construct cuboid/rectangle
	void makeBody(std::vector<double> width_length_depth, std::vector<double> angles, std::vector<double> centre,
		bool flex_rigid, bool deform, unsigned int group);
	// Method to construct filament
	void makeBody(unsigned int numbermarkers, std::vector<double> start_point, double fil_length, std::vector<double> angles, std::vector<int> BCs,
		bool flex_rigid, bool deform, unsigned int group);
	// Method to construct a 3D plate
	double makeBody(std::vector<double> width_length, double angle, std::vector<double> centre,
		bool flex_rigid, bool deform, unsigned int group, bool plate);

};

