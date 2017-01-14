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

#include "../inc/GridObj.h"

/// \brief	GridUnits.
///
///			This class contains static methods for unit conversion
///         (the only ones implemented are from m to cm and velocity from lattice units to m/s)
class GridUnits
{
public:
	GridUnits(){};
	~GridUnits(){};

	///< Convert from m to cm
	template <typename T>
	static T m2cm(const T meters) { return meters * 100; }

	// *****************************************************************************
	/// \brief	Velocity in lattice units to velocity in physical units.
	///
	///			Converts velocity component from lattice units to m/s. 
	///         It uses the L_PHYSICAL_U introduced by the user, dh and dt. 
	///			You can introduce any L_PHYSICAL_U you want, but the reference lenght (usualy the width of the domain)
	///         , the Re number and the LBM parameters will remain the same. So you will be implicitly changing the physical viscosity
	///         of your fluid when you change L_PHYSICAL_U 
	///
	/// \param ulat	Lattice velocity.
	/// \param currentGrid Pointer to the current grid. 
	/// \return physical velocity
	template <typename T>
	static T ulat2uphys(T ulat, GridObj* currentGrid)
	{
		return (ulat*currentGrid->dh*L_PHYSICAL_U) / currentGrid->dt;
	}

	
};

