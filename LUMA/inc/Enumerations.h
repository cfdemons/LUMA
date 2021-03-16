/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

#ifndef ENUM_H
#define ENUM_H

/// \enum eCartesianDirection
/// \brief Enumeration for directional options.
enum eCartesianDirection
{
	eXDirection,	///< X-direction
	eYDirection,	///< Y-direction
	eZDirection,	///< Z-direction
	eNoDirection	///< Null option
};

/// \enum eCartMinMax
/// \brief	Enumeration for the combination of eCartesianDirection and eMinMax 
///			as these are often used together to index arrays.
enum eCartMinMax
{
	eXMin,
	eXMax,
	eYMin,
	eYMax,
	eZMin,
	eZMax
};

/// \enum eLocationOnRank
/// \brief Enumeration indicating the location of a site when queried using isOnThisRank()
enum eLocationOnRank
{
	eNone,	///< No information provided (default).
	eCore,	///< Site on core (including send layer).
	eHalo	///< Site in halo (recv layer).

};

/// \enum eMinMax
/// \brief	Enumeration for minimum and maximum.
///
///			Some utility methods need to know whether they should be looking
///			at or for a maximum or minimum edge of a grid so we use this 
///			enumeration to specify.
enum eMinMax
{
	eMinimum,		///< Minimum
	eMaximum		///< Maximum
};


/// \enum eEdgeMinMax
/// \brief	Enumeration for the combination of Left and Right min and max edges.
enum eEdgeMinMax
{
	eLeftMin,
	eLeftMax,
	eRightMin,
	eRightMax
};

/// \enum  eType
/// \brief Lattice typing labels
enum eType
{
	eSolid,					///< Rigid, solid site with no-slip BC
	eFluid,					///< Fluid site
	eRefined,				///< Fluid site which is represented on a finer grid
	eTransitionToCoarser,	///< Fluid site coupled to a coarser grid
	eTransitionToFiner,		///< Fluid site coupled to a finer grid
	eBFL,					///< Site containing a BFL marker
	eVelocity,				///< Velocity boundary
	ePressure,				///< Pressure boundary
	eSlip,					///< Slip boundary
	eExtrapolateRight			///< Extrapolation boundary
};

/// \enum eWallLocation
/// \brief Enumeration to describe locations in terms of domain walls.
enum eWallLocation
{
	eLeftWall,		///< Left wall location
	eRightWall,		///< Right wall location
	eBottomWall,	///< Bottom wall location
	eTopWall,		///< Top wall location
	eFrontWall,		///< Front wall location
	eBackWall,		///< Back wall location

	eLeftBottomWall,	///< Left-bottom edge/corner
	eLeftTopWall,		///< Left-top edge/corner
	eRightBottomWall,	///< Right-bottom edge/corner
	eRightTopWall,		///< Right-top edge/corner

	eLeftFrontWall,		///< Left-front edge
	eLeftBackWall,		///< Left-back edge
	eRightFrontWall,	///< Right-front edge
	eRightBackWall,		///< Right-back edge

	eBottomFrontWall,	///< Bottom-front edge
	eBottomBackWall,	///< Bottom-front edge
	eTopFrontWall,		///< Top-front edge
	eTopBackWall,		///< Top-front edge

	eLeftBottomFrontWall,	///< Bottom-left-front corner
	eRightBottomFrontWall,	///< Bottom-right-front corner
	eLeftTopFrontWall,		///< Top-left-front corner
	eRightTopFrontWall,		///< Top-right-front corner
	eLeftBottomBackWall,	///< Bottom-left-back corner
	eRightBottomBackWall,	///< Bottom-right-back corner
	eLeftTopBackWall,		///< Top-left-back corner
	eRightTopBackWall,		///< Top-right-back corner
};

/// \enum  eIOFlag
/// \brief Flag for indicating write or read action for IO methods
enum eIOFlag
{
	eWrite,				///< Write to file
	eRead,				///< Read from file
};

/// \enum  eObjectType
/// \brief Specifies the type of body being processed.
enum eObjectType {
	eBBBCloud,	///< Bounce-back body
	eBFLCloud,	///< BFL body
	eIBBCloud	///< Immersed boundary body
};

/// \enum	eHdf5SlabType
///	\brief	Defines the type of storage arrangement of the variable in memory.
///
///			The write wrapper can then extract the data from memeory and write it
///			to an HDF5 file using a particular hyperslab selection.
enum eHdf5SlabType {
	eScalar,		///< 2/3D data	-- One variable per grid site
	eVector,		///< 2/3D data	-- L_DIMS variables per grid site
	eLatticeVector,	///< 2/3D data	-- L_NUM_VELS variables per grid site
	eProductVector,	///< 1D data	-- 3*L_DIMS-3 variables per grid site
	ePosX,			///< 1D data	-- Single L_dim vector per dimension
	ePosY,			///< 1D data	-- Single L_dim vector per dimension
	ePosZ			///< 1D data	-- Single L_dim vector per dimension
};

/// \enum  eMoveableType
/// \brief Specifies the whether body is movable, flexible or rigid.
enum eMoveableType {
	eFlexible,	///< Bounce-back body
	eMovable,	///< BFL body
	eRigid	///< Immersed boundary body
};

///	\enum eSDReturnType
///	\brief	Return types for smart decomposition methods.
enum eSDReturnType {
	eSDNoUnknown,
	eSDSnappedToEdges,
	eSDShrinkTarget,
	eSDProceed,
	eSDEarlyExit
};

#endif
