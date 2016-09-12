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

#include <vector>
#include "IBBody.h"
#include "Body.h"
#include "BFLBody.h"
#include "IVector.h"
#include "globalvars.h"
/// \enum  eObjectType
/// \brief Specifies the type of body being processed.
enum eObjectType {
	eBBBCloud,	///< Bounce-back body
	eBFLCloud,	///< BFL body
	eIBBCloud	///< Immersed boundary body
};

/// \brief	Object Manager class.
///
///			Class to manage all objects in the domain from creation through 
///			manipulation to destruction.
class ObjectManager
{

	// Make Grid a friend so boundary conditions can access the body data
	friend class GridObj;

	/** Members **/

private:

	// Bounce-back object fields
	double force_on_object_x = 0.0;			///< Instantaneous X-direction force on BB bodies in domain
	double force_on_object_y = 0.0;			///< Instantaneous Y-direction force on BB bodies in domain
	double force_on_object_z = 0.0;			///< Instantaneous Z-direction force on BB bodies in domain

	// Objects
	std::vector<Body<Marker>> bBody;		///< Array of default bodies
	std::vector<IBBody> iBody;				///< Array of immersed boundary bodies
	std::vector<BFLBody> pBody;				///< Array of BFL bodies

	/// Pointer to grid hierarchy
	GridObj* _Grids;

	/// Pre-stream distribution functions for applying BFL BCs
	IVector<double> f_prestream;

	/// Pointer to self
	static ObjectManager* me;


	/** Methods **/

private:
	ObjectManager(void);		///< Private constructor
	~ObjectManager(void);		///< Private destructor
	/// Overloaded constructor to set pointer to grid hierarchy
	ObjectManager(GridObj* _Grids);		

public:
	// Singleton design
	static ObjectManager *getInstance();		///< Get instance method
	static void destroyInstance();				///< Destroy instance method
	static ObjectManager *getInstance(GridObj* g);	///< Overloaded get instance passing in pointer to grid hierarchy

	// IBM methods //
	void ibm_apply(GridObj& g);					// Apply interpolate, compute and spread operations for all bodies and with GridObj g
	void ibm_build_body(int body_type);			// Build a new pre-fab IBM body
	void ibm_build_body(PCpts* _PCpts, GridObj *owner);	// Overloaded build-body to build from point cloud
	void ibm_initialise(GridObj& g);			// Initialise a built immersed body with support on the supplied grid
	double ibm_deltakernel(double rad, double dilation);	// Evaluate kernel (delta function approximation)
	void ibm_interpol(int ib, GridObj& g);			// Interpolation of velocity field on GridObj g onto markers of ib-th body
	void ibm_spread(int ib, GridObj& g);			// Spreading of restoring force from ib-th body to GridObj g
	void ibm_findsupport(int ib, int m, GridObj& g);		// Populates support information for the m-th marker of ib-th body with 
															// support nodes on GridObj g.
	void ibm_computeforce(int ib, GridObj& g);		// Compute restorative force at each marker in ib-th body using 
													// charcteristic time scale of GridObj g.
	double ibm_findepsilon(int ib, GridObj& g);		// Method to find epsilon weighting parameter for ib-th body given that the 
													// support points are on GridObj g.
	void ibm_move_bodies(GridObj& g);		// Update all IBBody positions and support based on GridObj g
	double ibm_bicgstab(std::vector< std::vector<double> >& Amatrix,
		std::vector<double>& bVector, std::vector<double>& epsilon,
						   double tolerance, int maxiterations);	// Biconjugate gradient stablised method for solving asymmetric 
																	// linear system required by finding epsilon

	// Flexible body methods
	void ibm_jacowire(int ib, GridObj& g);			// Computes the tension and position of a 2D inextensible, flexible iBody 
													// filament and normailises lengths based on spacing on GridObj g.
	void ibm_position_update(int ib, GridObj& g);	// Updates the position of deformable body markers and searches for support 
													// on GridObj g.
	void ibm_position_update_grp(int group, GridObj& g );	// Updates the positions of deformable bodies in a group based on the 
															// flexible group member position and searches for support on GridObj g.
	// Methods to solve the Jacobian system associated with Jacowire
	void ibm_banbks(double **a, long n, int m1, int m2, double **al,
		unsigned long indx[], double b[]);
	void ibm_bandec(double **a, long n, int m1, int m2, double **al,
		unsigned long indx[], double *d);


	// BFL methods //
	void bfl_build_body(int body_type);		// Build a new pre-fab bounce-back body
	void bfl_build_body(PCpts* _PCpts);		// Overload to build from point cloud data

	// Voxeliser methods //
	std::vector<int> getVoxInd(double x, double y, double z);
	int getVoxInd(double p);

	// Force calculation
	void computeLiftDrag(int i, int j, int k, GridObj *g);		// Compute force for BBB or BFLB

	// IO methods //
	void io_vtk_IBwriter(double tval);				// VTK body writer
	void io_write_body_pos(int timestep);			// Write out IBBody positions at specified timestep to text files
	void io_write_lift_drag(int timestep);			// Write out IBBody lift and drag at specified timestep
	void io_restart(bool IO_flag, int level);		// Restart read and write for IBBodies given grid level
	void io_readInCloud(PCpts* _PCpts, eObjectType objtype);	// Method to read in Point Cloud data
	void io_writeForceOnObject(double tval);			// Method to write object forces to a csv file
};

