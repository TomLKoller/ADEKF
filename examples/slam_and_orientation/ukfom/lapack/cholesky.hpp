/*
 *  Copyright (c) 2009, Rene Wagner
 *  All rights reserved.
 *
 *  Author: Rene Wagner <rw@nelianur.org>
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
 *   * Neither the name of Rene Wagner nor the names of any
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

#ifndef __UKFOM_LAPACK_CHOLESKY_HPP__
#define __UKFOM_LAPACK_CHOLESKY_HPP__

#include "lapack.h"

#include <Eigen/Core>

namespace ukfom {
namespace lapack {
using namespace Eigen;

template<size_t M>
class cholesky
{
public:
	cholesky(const Matrix<double, M, M> &m)
	{
		L_ = m;
		
		char UPLO = 'L';
		int N = L_.cols();
		int LDA = L_.stride();
		int INFO;

		dpotrf_(&UPLO, &N, L_.data(), &LDA, &INFO);

		spd_ = INFO == 0;

		// clear everything but the lower triangular matrix
		for (int j = 1; j < L_.cols(); ++j)
			for (int i = 0; i < j; ++i)
				L_(i,j) = 0;
	}

	
	const Matrix<double, M, M> &getL() const
	{
		if (!spd_)
			throw "not SPD";
		return L_;
	}

	bool isSPD() const
	{
		return spd_;
	}
	
private:
	Matrix<double, M, M> L_;
	bool spd_;
};

} // namespace lapack
} // namespace ukfom

#endif // __UKFOM_LAPACK_CHOLESKY_HPP__
