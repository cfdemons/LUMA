#include "../inc/stdafx.h"
#include "../inc/IBMarker.h"
#include "../inc/definitions.h"

// ***************************************************************************************************
// Custom constructor
IBMarker::IBMarker(double xPos, double yPos, double zPos, bool flex_rigid) {

	// Assign position and type
	this->flex_rigid = flex_rigid;
	this->position.push_back(xPos);
	this->position.push_back(yPos);
	this->position.push_back(zPos);

	// Stationary point
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);
	this->desired_vel.push_back(0.0);

	// Resize vectors
	this->fluid_vel.resize(dims);
	this->force_xyz.resize(dims);

}


// ***************************************************************************************************
