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

#ifndef IVECTOR_H
#define IVECTOR_H

#include "stdafx.h"

/// \brief	Index-collapsing vector class.
///
///			This class has all the behaviour of std::vector but 
///			has a overriden operator() to allow automatic flattening of indices 
///			before returning a reference of value at indexed location.
///			Needs to be able to accept different datatypes so templated.
template <typename GenTyp>
class IVector :	public std::vector<GenTyp>		// Define IVector class which inherits from std::vector
{
	
public:

	// C++ standard does not automatically inherit constructors
	/// Default constructor
	IVector( )
	{		
	}

	/// Default destructor
	~IVector( )
	{
	}

	/* Would like to inherit non-default constructor as required but the following
	 * actually crashes the compiler which is a compiler bug.
	 * using std::vector :: vector(size_type n, const value_type& val,
     *           const allocator_type& alloc = allocator_type());
	 */
	
	// Instead, define the required constructor myself
	/// \brief	Custom constructor taking type and value.
	///
	/// \param size the desired size of vector
	/// \param val the value to fill the new vector with
	IVector(size_t size, GenTyp val) {	

		this->resize(size);

		for (size_t i = 0; i < size; i++) {
			this->operator[] (i) = val;
		}

	}
	


	/*	
	 * :::: USE OF REFERENCES ::::
	 * Return reference (not the same as a pointer but similar) of the particular element in the 1D array.
	 * Reference is an alias for an object. Anbything done to the reference is done to the object itself.
	 * A pointer is an object itself, a reference is only a reference to an actual object.
	 * Here we delcare a reference to type GenTyp as the return type.
	 */

	// Define overloaded operator() for different array index flattening
	/// \brief	4D array index flatten.
	///
	///			Override of parentheses to auto-flatten indices to a single index.
	/// \param i the i index
	/// \param j the j index
	/// \param k the k index
	/// \param v the index in the fourth dimension
	/// \param j_max the number of j elements
	/// \param k_max the number of k elements
	/// \param v_max the number of elements in the fourth dimension
	/// \return GenTyp& a reference to the value at this position in the vector
	inline GenTyp& operator() (size_t i, size_t j, size_t k, size_t v, size_t j_max, size_t k_max, size_t v_max) {

		// Loop over vel then k then j then i
		size_t idx = v + (k*v_max) + (j*v_max*k_max) + (i*v_max*k_max*j_max);

		return this->operator[] (idx);

	}

	/// \brief	3D array index flatten.
	///
	///			Override of parentheses to auto-flatten indices to a single index.
	/// \param i the i index
	/// \param j the j index
	/// \param k the k index
	/// \param j_max the number of j elements
	/// \param k_max the number of k elements
	/// \return GenTyp& a reference to the value at this position in the vector
	inline GenTyp& operator() (size_t i, size_t j, size_t k, size_t j_max, size_t k_max) {

		// Loop over k then j then i
		size_t idx = k + (j*k_max) + (i*k_max*j_max);

		return this->operator[] (idx);

	}

	/// \brief	2D array index flatten.
	///
	/// \param i the i index
	/// \param j the j index
	/// \param j_max the number of j elements
	/// \return GenTyp& a reference to the value at this position in the vector
	inline GenTyp& operator() (size_t i, size_t j, size_t j_max) {
		
		// Loop over j then i
		size_t idx = j + (i*j_max);

		return this->operator[] (idx);

	}

};

#endif