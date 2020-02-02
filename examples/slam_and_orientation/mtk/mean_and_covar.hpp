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
 * @file mtk/mean_and_covar.hpp
 * @brief Functions to estimate mean value and covariance on manifolds.
 */

#ifndef MEAN_AND_COVAR_HPP_
#define MEAN_AND_COVAR_HPP_

#include <Eigen/Core>

namespace MTK {

/** 
 * \defgroup MeanCov Mean and Covariance Calculation
 * @todo provide functions, which only calculate either mean or covariance
 */
//@{


/**
 * Estimate mean value and covariance of a set of manifold values.
 * 
 * @tparam M    Manifold type. Must implement boxminus, boxplus and have 
 *              @c typedef scalar and @c enum DOF.
 * @tparam Cont Container Type. Elements must be convertible to @c M.
 * 
 * @param mean   reference to mean value (output)
 * @param cov    reference to covariance matrix (output)
 * @param values const reference to input container
 * @param max_it maximum number of iterations (optional).
 * 
 * Mean value and covariance are estimated using algorithm described in 
 * @cite{Hertzberg2011}
 */
template<class M, class Cont>
double mean_and_covariance(M& mean, Eigen::Matrix<typename M::scalar, M::DOF, M::DOF> &cov, 
                           const Cont &values, int max_it = 16)
{
	enum {DOF = M::DOF};
	typedef typename M::scalar scalar;
	mean = values[0];
	double res;
	int i=0;
	do {
		Eigen::Matrix<scalar, DOF, 1> mean_delta, delta;
		mean_delta.setZero();
		for (typename Cont::const_iterator Xi = values.begin(); Xi != values.end(); ++Xi)
		{
			Xi->boxminus(delta.data(), mean);
			mean_delta += delta;
		}
		mean_delta /= values.size();
		res = mean_delta.norm();
		mean.boxplus(mean_delta.data());
	} while (res > 1e-6 && ++i < max_it);
	
	
	cov.setZero();
	for (typename Cont::const_iterator Xi = values.begin(); Xi != values.end(); ++Xi)
	{
		Eigen::Matrix<scalar, DOF, 1> delta;
		Xi->boxminus(delta.data(), mean);
		cov += delta * delta.transpose();
	}
	cov *= 0.5;

	return res;
}

//@}

}  // namespace MTK


#endif /* MEAN_AND_COVAR_HPP_ */
