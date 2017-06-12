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

/* This file defines the constructors and methods for the immersed boundary body object.
*/

#include "../inc/stdafx.h"
#include "../inc/IBBody.h"
#include "../inc/IBMarker.h"
#include "../inc/PCpts.h"
#include "../inc/ObjectManager.h"

// ***************************************************************************************************
/// \brief Constructor which sets group ID to zero by default.
IBBody::IBBody()
{
	this->_Owner = nullptr;
	this->id = 0;
}

/// Default destructor
IBBody::~IBBody(void)
{
}

/******************************************************************************/
/// \brief	Initialiser for a general IB body.
///
///			This initialiser performs debugging IO and sets the property flags. 
///			It is called immediately after the appropriate constructor by all 
///			derived constructors.
///
///	\param	moveProperty	movable body flag passed on by the constructor.
void IBBody::initialise(eMoveableType moveProperty)
{

	// Set movable and flexible parameters
	if (moveProperty == eFlexible) {
		this->isFlexible = true;
		this->isMovable = true;
	}
	else if (moveProperty == eMovable) {
		this->isFlexible = false;
		this->isMovable = true;
	}
	else if (moveProperty == eRigid) {
		this->isFlexible = false;
		this->isMovable = false;
	}

	// Delete any markers which are on the receiver layer as not needed for IBM
#ifdef L_BUILD_FOR_MPI
	*GridUtils::logfile << "Deleting IB markers which exist on receiver layer..." << std::endl;
	deleteRecvLayerMarkers();
#endif

}


/*********************************************/
/// \brief Custom constructor to populate body from array of points.
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param _PCpts			pointer to point cloud data
/// \param moveProperty		determines if body is moveable, flexible or rigid
/// \param clamped			boundary condition for fixed end (only relevant if body is flexible)
IBBody::IBBody(GridObj* g, int bodyID, PCpts* _PCpts, eMoveableType moveProperty, bool clamped)
	: Body(g, bodyID, _PCpts)
{
	initialise(moveProperty);	
}


/*********************************************/
/// \brief 	Custom constructor for building prefab filament
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param start_position	start position of base of filament
/// \param length			length of filament
/// \param angles			angle of filament
/// \param moveProperty		determines if body is moveable, flexible or rigid
/// \param clamped			boundary condition for fixed end (only relevant if body is flexible)
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &start_position,
		double length, std::vector<double> &angles, eMoveableType moveProperty, bool clamped)
		: Body(g, bodyID, start_position, length, angles)
{
	initialise(moveProperty);
}


/*********************************************/
/// \brief 	Custom constructor for building prefab circle/sphere
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of circle
/// \param radius			radius of circle
/// \param moveProperty		determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		double radius, eMoveableType moveProperty)
		: Body(g, bodyID, centre_point, radius)
{
	initialise(moveProperty);
}

/*********************************************/
/// \brief 	Custom constructor for building square/cuboid
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of square
/// \param dimensions		dimensions of square
/// \param angles			angle of square
/// \param moveProperty		determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		std::vector<double> &dimensions, std::vector<double> &angles, eMoveableType moveProperty)
		: Body(g, bodyID, centre_point, dimensions, angles)
{
	initialise(moveProperty);
}
