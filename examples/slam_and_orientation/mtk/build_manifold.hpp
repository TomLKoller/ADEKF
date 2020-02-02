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
 * @file mtk/build_manifold.hpp
 * @brief Macro to automatically construct compound manifolds.
 * 
 */
#ifndef MTK_AUTOCONSTRUCT_HPP_
#define MTK_AUTOCONSTRUCT_HPP_


#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/cat.hpp>

#include "src/SubManifold.hpp"
#include "startIdx.hpp"

#ifndef PARSED_BY_DOXYGEN
//////// internals //////

#define MTK_APPLY_MACRO_ON_TUPLE(r, macro, tuple) macro tuple

#define MTK_TRANSFORM_COMMA(macro, entries) BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM_S(1, MTK_APPLY_MACRO_ON_TUPLE, macro, entries))

#define MTK_TRANSFORM(macro, entries) BOOST_PP_SEQ_FOR_EACH_R(1, MTK_APPLY_MACRO_ON_TUPLE, macro, entries)

#define MTK_CONSTRUCTOR_ARG(  type, id) const type& id = type()
#define MTK_CONSTRUCTOR_COPY( type, id) id(id)
#define MTK_BOXPLUS(          type, id) id.boxplus(MTK::subvector(__vec, &self::id), __scale);
#define MTK_BOXMINUS(         type, id) id.boxminus(MTK::subvector(__res, &self::id), __oth.id);
#define MTK_OSTREAM(          type, id) << __var.id << " "
#define MTK_ISTREAM(          type, id) >> __var.id

#define MTK_SUBVARLIST(seq) \
BOOST_PP_FOR_1( \
		( \
				BOOST_PP_SEQ_SIZE(seq), \
				BOOST_PP_SEQ_HEAD(seq), \
				BOOST_PP_SEQ_TAIL(seq) (~), \
				0 ),\
		MTK_ENTRIES_TEST, MTK_ENTRIES_NEXT, MTK_ENTRIES_OUTPUT)

#define MTK_PUT_TYPE(type, id, dof) \
	MTK::SubManifold<type, dof> id; 
#define MTK_PUT_TYPE_AND_ENUM(type, id, dof) \
	MTK_PUT_TYPE(type, id, dof) \
	enum {DOF = type::DOF + dof}; \
	typedef type::scalar scalar;

#define MTK_ENTRIES_OUTPUT(r, state) MTK_ENTRIES_OUTPUT_I state
#define MTK_ENTRIES_OUTPUT_I(s, head, seq, dof) \
	MTK_APPLY_MACRO_ON_TUPLE(~, \
		BOOST_PP_IF(BOOST_PP_DEC(s), MTK_PUT_TYPE, MTK_PUT_TYPE_AND_ENUM), \
		( BOOST_PP_TUPLE_REM_2 head, dof)) 

#define MTK_ENTRIES_TEST(r, state) MTK_TUPLE_ELEM_4_0 state

//! this used to be BOOST_PP_TUPLE_ELEM_4_0:
#define MTK_TUPLE_ELEM_4_0(a,b,c,d) a

#define MTK_ENTRIES_NEXT(r, state) MTK_ENTRIES_NEXT_I state
#define MTK_ENTRIES_NEXT_I(len, head, seq, dof) ( \
		BOOST_PP_DEC(len), \
		BOOST_PP_SEQ_HEAD(seq), \
		BOOST_PP_SEQ_TAIL(seq), \
		dof + BOOST_PP_TUPLE_ELEM_2_0 head::DOF)

#endif /* not PARSED_BY_DOXYGEN */


/**
 * Construct a manifold.
 * @param name is the class-name of the manifold, 
 * @param entries is the list of sub manifolds 
 * 
 * Entries must be given in a list like this:
 * @code
 * typedef MTK::trafo<MTK::SO3<double> > Pose;
 * typedef MTK::vect<double, 3> Vec3;
 * MTK_BUILD_MANIFOLD(imu_state,
 *    ((Pose, pose))
 *    ((Vec3, vel))
 *    ((Vec3, acc_bias))
 * )
 * @endcode
 * Whitespace is optional, but the double parentheses are necessary.
 * Construction is done entirely in preprocessor.
 * After construction @a name is also a manifold. Its members can be 
 * accessed by names given in @a entries.
 * 
 * @note Variable types are not allowed to have commas, thus types like
 *       @c vect<double, 3> need to be typedef'ed ahead.
 */
#define MTK_BUILD_MANIFOLD(name, entries) \
struct name { \
	typedef name self; \
	MTK_SUBVARLIST(entries) \
	name ( \
		MTK_TRANSFORM_COMMA(MTK_CONSTRUCTOR_ARG, entries) \
		) : \
		MTK_TRANSFORM_COMMA(MTK_CONSTRUCTOR_COPY, entries) {}\
	int getDOF() const { return DOF; } \
	void boxplus(const MTK::vectview<const scalar, DOF> & __vec, scalar __scale = 1 ) { \
		MTK_TRANSFORM(MTK_BOXPLUS, entries) \
	} \
	void boxminus(MTK::vectview<scalar,DOF> __res, const name& __oth) const { \
		MTK_TRANSFORM(MTK_BOXMINUS, entries) \
	} \
	friend std::ostream& operator<<(std::ostream& __os, const name& __var){ \
		return __os MTK_TRANSFORM(MTK_OSTREAM, entries); \
	} \
	friend std::istream& operator>>(std::istream& __is, name& __var){ \
		return __is MTK_TRANSFORM(MTK_ISTREAM, entries); \
	} \
};



#endif /*MTK_AUTOCONSTRUCT_HPP_*/
