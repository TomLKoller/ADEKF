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
 * @file mtk/types/pose.hpp
 * @brief Definition of poses (or transformations between them).
 * 
 * This file replaces the old MakePose.h file. Now poses just consist of 
 * orientation and translation, but they can be components of larger manifolds.
 * 
 * @todo Maybe @c trafo should be renamed to @c SEn?
 */

#ifndef POSE_HPP_
#define POSE_HPP_

#include "../startIdx.hpp"

namespace MTK {

/**
 * Representation of Standard Euclidean Group as manifold.
 * 
 * Transformations can be combined using operators @c *, @c / and 
 * @c % [corresponding to MATLAB's @c \\ ], as if they where 
 * @f$ (n+1) \times (n+1) @f$ matrices. In the same manner they can be 
 * multiplied by vectors (interpreted as projective vectors).
 * 
 * 
 * @tparam RotationType rotation type which must be a manifold and must
 *                      implement operators @c *, @c / and @c % [corresponding to MATLAB's @c \\] 
 *                      for @c RotationType and operators @c * and @c % for @c VectType.
 * @tparam Trans1st     decides whether the translational or rotational part shall be used 
 *                      first when applying @c boxplus or @c boxminus.
 * @tparam VectType     type of the vector part, defaults to @c RotationType::vect_type.
 */
template<class RotationType, bool Trans1st = true, class VectType = typename RotationType::vect_type>
class trafo
{
	enum {RotationDOF = RotationType::DOF,             TranslationDOF = VectType::DOF};
	enum {RotationIdx = Trans1st ? TranslationDOF : 0, TranslationIdx = Trans1st ? 0 : RotationDOF};
	
public:
	typedef RotationType rotation_type;
	typedef VectType vect_type;
	typedef typename rotation_type::scalar scalar;
	enum {DOF = RotationDOF + TranslationDOF};
	
	//! Rotational part
	MTK::SubManifold<rotation_type, RotationIdx> orient;
	//! Translational part
	MTK::SubManifold<vect_type, TranslationIdx> pos;
	
	//! construct from a rotation and translation
	trafo(const rotation_type &rot = rotation_type(), const vect_type &vec = vect_type()) : orient(rot), pos(vec) {}
	//! construct from a translation and rotation
	trafo(const vect_type &vec, const rotation_type &rot = rotation_type()) : orient(rot), pos(vec) {}
	
	//! construct form a different @c trafo
	template<class RotType2, bool T, class VectType2>
	trafo(const trafo<RotType2, T, VectType2> &tr) : orient(tr.orient), pos(tr.pos) {}
	
	void boxplus(MTK::vectview<const scalar, DOF> vec, scalar _scale = 1)
	{
		orient.boxplus(MTK::subvector(vec, &trafo::orient), _scale);
		pos.boxplus(MTK::subvector(vec, &trafo::pos), _scale);
	}
	void boxminus(MTK::vectview<scalar, DOF> ret, const trafo &other) const
	{
		orient.boxminus(MTK::subvector(ret, &trafo::orient), other.orient);
		pos.boxminus(MTK::subvector(ret, &trafo::pos), other.pos);
	}
	
	friend std::ostream& operator<<(std::ostream &os, const trafo<RotationType, Trans1st, VectType> &tr){
		if(Trans1st)
			return os << tr.pos << " " << tr.orient;
		else
			return os << tr.orient << " " << tr.pos;
	}
	friend std::istream& operator>>(std::istream &is,  trafo<RotationType, Trans1st, VectType> &tr){
		if(Trans1st)
			return is >> tr.pos >> tr.orient;
		else
			return is >> tr.orient >> tr.pos;
	}

	
	vect_type operator*(const vect_type& vec) const
	{
		return local2World(vec);
	}
	vect_type local2World(const vect_type& vec) const
	{
		return orient * vec + pos;
	}
	
	vect_type operator%(const vect_type& vec) const
	{
		return world2Local(vec);
	}
	vect_type world2Local(const vect_type& vec) const
	{
		return orient % (vec - pos);
	}
	
	trafo operator*(const trafo &traf) const
	{
		return local2World(traf);
	}
	trafo local2World(const trafo &traf) const
	{
		return trafo(orient * traf.orient, operator*(traf.pos));
	}

	trafo operator%(const trafo &traf) const
	{
		return world2Local(traf);
	}
	trafo world2Local(const trafo &traf) const
	{
		return trafo(orient % traf.orient, operator%(traf.pos));
	}
	
	trafo operator/(const trafo &traf) const
	{
		rotation_type rot = orient / traf.orient;
		return trafo(rot, pos - rot * traf.pos);
	}
	
	trafo inverse() const
	{
		return world2Local(trafo());
	}
	
	
	// FIXME fatal also implement /, % and provide unit-tests
	trafo operator*(const rotation_type& rot) const {
		return trafo(orient * rot, pos);
	}
	friend 
	trafo operator*(const rotation_type& rot, const trafo& traf) {
		return trafo(rot * traf.orient, rot*traf.pos);
	}
};



namespace internal {
template<class RotationType, bool Trans1st, class VectType>
struct UnalignedType<trafo<RotationType, Trans1st, VectType > >{
	typedef trafo<typename UnalignedType<RotationType>::type, Trans1st, typename UnalignedType<VectType>::type > type;
};


}  // namespace internal



}  // namespace MTK


#endif /* POSE_HPP_ */
