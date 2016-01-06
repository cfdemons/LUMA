#pragma once

#include <vector>
#include "IBBody.h"


class ObjectManager
{

	/** Members **/

private:

	// IBM objects
	static std::vector<IBBody> iBody;		// Array of immersed boundary bodies

	/** Methods **/

public:

	// Constructor and Destructor
	ObjectManager(void);
	~ObjectManager(void);

	// IBM methods
	static void ibm_apply(GridObj& g);					// Apply interpolate, compute and spread operations for all bodies and with GridObj g
	static void ibm_build_body(int body_type);			// Build a new pre-fab body
	void ibm_initialise(GridObj& g);					// Initialise a built immersed body with support on the supplied grid
	static double ibm_deltakernel(double rad, double dilation);		// Evaluate kernel (delta function approximation)
	static void ibm_interpol(unsigned int ib, GridObj& g);			// Interpolation of velocity field on GridObj g onto markers of ib-th body
	static void ibm_spread(unsigned int ib, GridObj& g);			// Spreading of restoring force from ib-th body to GridObj g
	static void ibm_findsupport(unsigned int ib, 
		unsigned int m, GridObj& g);					// Populates support information for the m-th marker of ib-th body with support nodes 
														// on GridObj g.
	static void ibm_computeforce(unsigned int ib, GridObj& g);		// Compute restorative force at each marker in ib-th body using 
																	// charcteristic time scale of GridObj g.
	static double ibm_findepsilon(unsigned int ib, GridObj& g);		// Method to find epsilon weighting parameter for ib-th body given that the 
																	// support points are on GridObj g.
	static void ibm_move_bodies(GridObj& g);			// Update all IBBody positions and support based on GridObj g
	// Biconjugate gradient stablised method for solving asymmetric linear system required by finding epsilon
	static double ibm_bicgstab(std::vector< std::vector<double> >& Amatrix, std::vector<double>& bVector, std::vector<double>& epsilon,
						   double tolerance, unsigned int maxiterations);

	// Flexible body methods
	static void ibm_jacowire(unsigned int ib, GridObj& g);			// Computes the tension and position of a 2D inextensible, flexible iBody 
																	// filament and normailises lengths based on spacing on GridObj g.
	static void ibm_position_update(unsigned int ib, GridObj& g);	// Updates the position of deformable body markers and searches for support 
																	// on GridObj g.
	static void ibm_position_update_grp(unsigned int group, GridObj& g );	// Updates the positions of deformable bodies in a group based on the 
																			// flexible group member position and searches for support on GridObj g.
	// Methods to solve the Jacobian system associated with Jacowire
	static void ibm_banbks(double **a, unsigned long n, unsigned int m1, unsigned int m2, double **al,
		unsigned long indx[], double b[]);
	static void ibm_bandec(double **a, unsigned long n, unsigned int m1, unsigned int m2, double **al,
		unsigned long indx[], double *d);

	// IBM IO methods
	void io_vtk_IBwriter(double tval);							// VTK body writer
	void io_write_body_pos(unsigned int timestep);				// Write out IBBody positions at specified timestep to text files
	void io_write_lift_drag(unsigned int timestep);				// Write out IBBody lift and drag at specified timestep
	static void io_restart(bool IO_flag, unsigned int level);	// Restart read and write for IBBodies given grid level


};

