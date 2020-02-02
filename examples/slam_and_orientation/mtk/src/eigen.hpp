/*
 *  Copyright (c) 2008--2011, Universitaet Bremen
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
/**
 * @file mtk/src/eigen.hpp
 * @brief Eigen2/Eigen3 compatibility layer
 */

#ifndef EIGEN_HPP_
#define EIGEN_HPP_

#define EIGEN_NO_DEBUG

//#define SLOM_EXPERIMENTAL

#include <Eigen/Core>
#include <boost/static_assert.hpp>


// TODO Find out which beta version of Eigen3 works and give errors for others
#if EIGEN_VERSION_AT_LEAST(2,91,0)

#ifdef MTK_EXPERIMENTAL
#define MTK_EIGEN 299
// This will be obsolete when Eigen::Map gets optimized:
#include "MTKMap.hpp"

#else

#define MTK_EIGEN 300

#endif /* MTK_EXPERIMENTAL */

#else

#include <Eigen/Array>

#define MTK_EIGEN 200

// For debugging/compatibility to Eigen3:
#define EIGEN_ASM_COMMENT(X)  asm("#" X)

#endif

namespace MTK {

namespace internal {

template<class Base, class T1, class T2>
struct CovBlock{
#if MTK_EIGEN >= 300
	// Eigen3 has const correct block return types:
	typedef typename Eigen::Block<Eigen::Matrix<typename Base::scalar, Base::DOF, Base::DOF>, T1::DOF, T2::DOF> Type;
	typedef typename Eigen::Block<const Eigen::Matrix<typename Base::scalar, Base::DOF, Base::DOF>, T1::DOF, T2::DOF> ConstType;
#else /* EIGEN 2 */
	typedef typename Eigen::BlockReturnType<Eigen::Matrix<typename Base::scalar, Base::DOF, Base::DOF>, T1::DOF, T2::DOF>::Type Type;
	typedef const Type ConstType;
#endif
};


template<class scalar, int dim>
struct VectviewBase {
	typedef Eigen::Matrix<scalar, dim, 1> matrix_type;
#if MTK_EIGEN >= 300
	typedef typename matrix_type::MapType Type;
	typedef typename matrix_type::ConstMapType ConstType;
#elif MTK_EIGEN >=299
	// FIXME EXPERIMENTAL!!!
	typedef Eigen::VectMap<scalar, dim> Type;
	typedef Eigen::VectMap<const scalar, dim> ConstType;
#else /* EIGEN 2 */
	typedef typename matrix_type::UnalignedMapType Type;
	typedef Type ConstType;
#endif
};




template<class T>
struct UnalignedType {
	typedef T type;
};



}  // namespace internal

}  // namespace MTK


#define MTK_ALWAYS_INLINE EIGEN_ALWAYS_INLINE
#define MTK_DEPRECATED __attribute__((deprecated))

#endif /* EIGEN_HPP_ */
