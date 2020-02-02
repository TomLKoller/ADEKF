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
 * @file mtk/src/SubManifold.hpp
 * @brief Defines the SubManifold class
 */


#ifndef SUBMANIFOLD_HPP_
#define SUBMANIFOLD_HPP_


#include "vectview.hpp"


namespace MTK {

/**
 * @ingroup SubManifolds
 * Helper class for compound manifolds. 
 * This class wraps a manifold T and provides an enum IDX refering to the 
 * index of the SubManifold within the compound manifold. 
 *  
 * Memberpointers to a submanifold can be used for @ref SubManifolds "functions accessing submanifolds".
 * 
 * @tparam T   The manifold type of the sub-type
 * @tparam idx The index of the sub-type within the compound manifold
 */
template<class T, int idx>
struct SubManifold : public T 
{
	enum {IDX = idx /*!< index of the sub-type within the compound manifold */ };
	//! manifold type
	typedef T type;
	
	//! Construct from derived type
	template<class X>
	explicit
	SubManifold(const X& t) : T(t) {}
	
	//! Construct from internal type
	SubManifold(const T& t) : T(t) {}
	
	//! inherit assignment operator
	using T::operator=;
	
};

}  // namespace MTK


#endif /* SUBMANIFOLD_HPP_ */
