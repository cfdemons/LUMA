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

#include "../inc/stdafx.h"
#include "../inc/IBInfo.h"

// Default constructor
epsCalcMarkerClass::epsCalcMarkerClass() {

	// Default values
	bodyID = 0;
	local_area = 0.0;
	dilation = 0.0;
}


/// \brief Custom constructor for epsilon calculation class.
///	\param	ib				id of current body.
///	\param	positionIn		position of IB marker.
///	\param	areaIn			local area of IB marker.
///	\param	dilationIn		dilation parameter for IB marker.
///	\param	supp_position	positions of all support points for that IB marker.
///	\param	deltavalIn		delta value for all support points for that IB marker.
epsCalcMarkerClass::epsCalcMarkerClass(int bodyIDIn, std::vector<double> positionIn, double areaIn, double dilationIn, std::vector<std::vector<double>> supp_position, std::vector<double> deltavalIn) {

	// Assign values
	bodyID = bodyIDIn;
	local_area = areaIn;
	dilation = dilationIn;

	// Resize vectors and assign
	position.resize(positionIn.size());
	deltaval.resize(deltavalIn.size());
	position = positionIn;
	deltaval = deltavalIn;

	// Assign position of support points
	for (int i = 0; i < deltaval.size(); i++) {
		x.push_back(supp_position[i][eXDirection]);
		y.push_back(supp_position[i][eYDirection]);

#if (L_DIMS == 3)
		z.push_back(supp_position[i][eZDirection]);
#endif
	}
}
