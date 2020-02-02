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
 * @file mtk/src/vectview.hpp
 * @brief Wrapper class around a pointer used as interface for plain vectors.
 */

#ifndef VECTVIEW_HPP_
#define VECTVIEW_HPP_

#include "eigen.hpp"

namespace MTK {

/**
 * A view to a vector.
 * Essentially, @c vectview is only a pointer to @c scalar but can be used directly in @c Eigen expressions.
 * The dimension of the vector is given as template parameter and type-checked when used in expressions.
 * Data has to be modifiable.
 * 
 * @tparam scalar Scalar type of the vector.
 * @tparam dim    Dimension of the vector.
 * 
 * @todo @c vectview can be replaced by simple inheritance of @c Eigen::Map, as soon as they get const-correct
 */
template<class scalar, int dim>
class vectview : public internal::VectviewBase<scalar, dim>::Type {
	typedef internal::VectviewBase<scalar, dim> VectviewBase;
public:
	//! plain matrix type
	typedef typename VectviewBase::matrix_type matrix_type;
	//! base type
	typedef typename VectviewBase::Type base;
	//! construct from pointer
	explicit
	vectview(scalar* data, int dim_=dim) : base(data, dim_) {}
	//! construct from plain matrix
	vectview(matrix_type& m) : base(m.data(), m.size()) {}
	//! construct from another @c vectview
	vectview(const vectview &v) : base(v) {}
	//! construct from Eigen::Block:
#if MTK_EIGEN >= 300
	template<class Base>
	vectview(Eigen::VectorBlock<Base, dim> block) : base(&block.coeffRef(0), block.size()) {}
	template<class Base, bool PacketAccess>
	vectview(Eigen::Block<Base, dim, 1, PacketAccess> block) : base(&block.coeffRef(0), block.size()) {}
#else
	template<class Base, int PacketAccess, int DirectAccessStatus>
	vectview(Eigen::Block<Base, dim, 1, PacketAccess, DirectAccessStatus> block) : base(&block.coeffRef(0)) {}
#endif
	//! inherit assignment operator
	using base::operator=;
	//! data pointer
	scalar* data() {return const_cast<scalar*>(base::data());}
};

/**
 * @c const version of @c vectview.
 * Compared to @c Eigen::Map this implementation is const correct, i.e.,
 * data will not be modifiable using this view.
 * 
 * @tparam scalar Scalar type of the vector.
 * @tparam dim    Dimension of the vector.
 * 
 * @sa vectview
 */
template<class scalar, int dim>
class vectview<const scalar, dim> : public internal::VectviewBase<scalar, dim>::ConstType {
	typedef internal::VectviewBase<scalar, dim> VectviewBase;
public:
	//! plain matrix type
	typedef typename VectviewBase::matrix_type matrix_type;
	//! base type
	typedef typename VectviewBase::ConstType base;
	//! construct from const pointer
	explicit
	vectview(const scalar* data, int dim_ = dim) : base(data, dim_) {}
	//! construct from column vector
	template<int options>
	vectview(const Eigen::Matrix<scalar, dim, 1, options>& m) : base(m.data()) {}
	//! construct from row vector
	template<int options, int phony>
	vectview(const Eigen::Matrix<scalar, 1, dim, options, phony>& m) : base(m.data()) {}
	//! construct from another @c vectview
	vectview(vectview<scalar, dim> x) : base(x.data()) {}
	//! construct from base
	vectview(const base &x) : base(x) {}
	/**
	 * Construct from Block
	 * @todo adapt this, when Block gets const-correct
	 */
#if MTK_EIGEN >= 300
	template<class Base>
	vectview(Eigen::VectorBlock<Base, dim> block) : base(&block.coeffRef(0)) {}
	template<class Base, bool PacketAccess>
	vectview(Eigen::Block<Base, dim, 1, PacketAccess> block) : base(&block.coeffRef(0)) {}
#else
	template<class Base, int PacketAccess, int DirectAccessStatus>
	vectview(Eigen::Block<Base, dim, 1, PacketAccess, DirectAccessStatus> block) : base(&block.coeffRef(0)) {}
#endif
//	//FIXME construct from temporary expression, this is not save!
//	template<typename OtherDerived>
//	vectview(const Eigen::MatrixBase<OtherDerived>& other) : base(other.eval().data()) {
//		std::cout << __PRETTY_FUNCTION__ << ", data: " << base::data() << std::endl;
//	}
	//! constant index operator
//	const scalar operator[](int idx) const {return base::operator[](idx);}
	//! constant index operator
//	const scalar operator()(int idx) const {return base::operator()(idx);}
private:
	void operator=(const vectview& DONT_USE) const;
	template<class DONT_USE> void operator= (const DONT_USE&) const;
	template<class DONT_USE> void operator+=(const DONT_USE&) const;
	template<class DONT_USE> void operator-=(const DONT_USE&) const;
	template<class DONT_USE> void operator*=(const DONT_USE&) const;
	template<class DONT_USE> void operator/=(const DONT_USE&) const;
	const scalar& operator()(int, int); // DONT USE!
	const scalar& coeffRef(int, int = 0);
};


} // namespace MTK

#endif /* VECTVIEW_HPP_ */
