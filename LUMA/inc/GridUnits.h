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
	GridUnits();
	~GridUnits();

	///< Convert from m to cm
	template <typename T>
	static T m2cm(const T meters) { return meters * 100; }

	// *****************************************************************************
	/// \brief	Velocity in lattice units to velocity in physical units.
	///
	///			Converts velocity component from lattice units to m/s. 
	///         It assumes dx=dt , so vlattice = vadimensional (i.e. v for the dummy fluid with non physical viscosity)
	///         Then computes the physical velocity using the physical viscosity in definitions.h
	///
	/// \param ulat	Lattice velocity.
	/// \param currentGrid Pointer to the current grid. 
	/// \return physical velocity
	template <typename T>
	static T ulat2uphys(T ulat, GridObj* currentGrid)
	{
		T lphys = (currentGrid->nu*L_Re*currentGrid->dx) / L_u_ref;
		
		return ((L_Re*L_nu) / lphys)*ulat;
	}

	
};

