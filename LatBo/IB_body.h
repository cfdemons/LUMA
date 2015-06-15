#pragma once

#include <vector>
#include "IB_marker.h"

class IB_body {

	// Make GridObj a friend class so it can access the protected data of IB_body objects
	friend class GridObj;

public:
	// Default constructor and destructor
	IB_body(void);
	~IB_body(void);

protected:

	/*	
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/

	std::vector<IB_marker> markers;		// Array of particles which make up the body
	double spacing;	// Spacing of the Lagrnage markers
	
	/*	
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

public:
	// Add Lagrange marker to the body
	void addMarker(double x, double y, double z, bool flex_rigid);

	// Prefab body building methods
	void makeBody(double radius, std::vector<double> centre, bool flex_rigid);		// Method to construct sphere/circle
	void makeBody(std::vector<double> width_length_depth, std::vector<double> centre, bool flex_rigid);		// Method to construct cuboid/rectangle

};

