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

#ifndef ENUM_H
#define ENUM_H

/// \enum eCartesianDirection
/// \brief Enumeration for directional options.
enum eCartesianDirection
{
	eXDirection,	///< X-direction
	eYDirection,	///< Y-direction
	eZDirection		///< Z-direction
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

///	\enum eIBInfoType
///	\brief	Type of container required.
enum eIBInfoType {
	eIBDeltaSum,
	eIBEpsilon,
	eIBVelocityInterpolation,
	eIBVelocitySpreading,
	eIBMarkerPositions
};

/// \enum  eType
/// \brief Lattice typing labels
enum eType
{
	eSolid,					///< Rigid, solid site
	eFluid,					///< Fluid site
	eRefined,				///< Fluid site which is represented on a finer grid
	eTransitionToCoarser,	///< Fluid site coupled to a coarser grid
	eTransitionToFiner,		///< Fluid site coupled to a finer grid
	eBFL,					///< Site containing a BFL marker
	eSymmetry,				///< Symmetry boundary
	eInlet,					///< Inlet boundary
	eOutlet,				///< Outlet boundary
	eRefinedSolid,			///< Rigid, solid site represented on a finer grid
	eRefinedSymmetry,		///< Symmtery boundary represented on a finer grid
	eRefinedInlet			///< Inlet site represented on a finer grid
};

/// \enum  eBCType
/// \brief Flag for indicating which BCs to apply
enum eBCType
{
	eBCAll,				///< Apply all BCs
	eBCSolidSymmetry,	///< Apply just solid and symmetry BCs
	eBCInlet,			///< Apply just inlet BCs
	eBCOutlet,			///< Apply just outlet BCs
	eBCInletOutlet,		///< Apply inlet and outlet BCs
	eBCBFL				///< Apply just BFL BCs
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