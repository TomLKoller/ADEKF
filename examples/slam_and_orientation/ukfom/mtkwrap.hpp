/*
 *  Copyright (c) 2009, 2010, Universitaet Bremen
 *  All rights reserved.
 *
 *  Author: Christoph Hertzberg <chtz@informatik.uni-bremen.de>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Universitaet Bremen nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __UKF_MTKWRAP_HPP__
#define __UKF_MTKWRAP_HPP__

#include <cassert>

#include <Eigen/Core>
//#include <Eigen/LU> 

//#include <Eigen/QR>


namespace ukfom {



// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN


/**
 * mtkwrap<M> wraps an MTK-Manifold M to a ukf-compatible manifold.
 * M has to have an enum DOF and implement the methods boxplus and boxminus.
 */
template<class M>
struct mtkwrap : public M{
	typedef mtkwrap<M> self;
public:
	typedef double scalar; // MTK only works with double
	typedef scalar scalar_type;

	enum {
		DOF = M::DOF
	};
	
	typedef Matrix<scalar_type, DOF, 1> vectorized_type;

	mtkwrap(const M &m=M()) : M(m) {}
	

	/*
	 * manifold operator (+)
	 *
	 */
	self& operator+=(const vectorized_type &delta_state)
	{
		assert(delta_state.stride() == DOF);
		M::boxplus(delta_state.data());
		return *this;
	}

	const self operator+(const vectorized_type &delta_state) const
	{
		self result = *this;
		result += delta_state;
		
		return result;
	}

	/*
	 * manifold operator (-)
	 */
	const vectorized_type operator-(const self &other) const
	{
		vectorized_type result;
		assert(result.stride()==DOF);
		M::boxminus(result.data(), other);

		return result;
	}

	bool operator==(const self &other) const
	{
		vectorized_type diff = (*this) - other;
		return diff.isZero(1e-12);
	}

	bool operator!=(const self &other) const
	{
		return !(*this == other);
	}

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace ukfom
#endif /* __UKF_MTKWRAP_HPP__ */
