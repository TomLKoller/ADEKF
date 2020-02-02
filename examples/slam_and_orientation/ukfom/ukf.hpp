/*
 *  Copyright (c) 2009, Rene Wagner
 *  Copyright (c) 2010, 2011 DFKI GmbH
 *  All rights reserved.
 *
 *  Author: Rene Wagner <rene.wagner@dfki.de>
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
 *   * Neither the name of the DFKI GmbH nor the names of its
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

#ifndef __UKFOM_UKF_HPP__
#define __UKFOM_UKF_HPP__

#include <vector>
#include <algorithm>
#include <numeric>

#include <boost/bind.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

#include <Eigen/QR>

#include "lapack/cholesky.hpp"
#include "traits/dof.hpp"
#include "util.hpp"

namespace ukfom {
	
// import most common Eigen types 
using namespace Eigen;

template <typename state>
class ukf {
	typedef ukf self;

	enum {
		n = state::DOF
	};
	
public:
	
	typedef typename state::scalar_type scalar_type;
	typedef typename state::vectorized_type vectorized_state;
	typedef Matrix<scalar_type, int(state::DOF), int(state::DOF)> cov;
	typedef std::vector<state> state_vector;
	
	ukf(const state &mu,
		const cov &sigma)
		: mu_(mu),
		  sigma_(sigma)
	{
	}

	template<typename ProcessModel>
	void predict(ProcessModel g, const cov &R)
	{
		predict(g, boost::bind(id<cov>, R));
	}

	template<typename ProcessModel, typename ProcessNoiseCovariance>
	void predict(ProcessModel g, ProcessNoiseCovariance R)
	{
		state_vector X(2 * n + 1);
		generate_sigma_points(mu_, sigma_, X);
		
		std::transform(X.begin(), X.end(), X.begin(), g);

		mu_ = sigma_points_mean(X);

		//std::cout << "mu':" << std::endl << mu_ << std::endl;
		
		sigma_ = sigma_points_cov<state::DOF>(mu_, X) + R();
	}

	template<typename Measurement,
			 typename MeasurementModel,
			 typename MeasurementNoiseCovariance>
	void update(const Measurement &z,
				MeasurementModel h,
				MeasurementNoiseCovariance Q)
	{
		update(z, h, Q,
			   accept_any_mahalanobis_distance<scalar_type>);
	}
	
	template<typename Measurement,
			 typename MeasurementModel>
	void update(const Measurement &z,
				MeasurementModel h,
				const Eigen::Matrix<scalar_type, dof<Measurement>::value, dof<Measurement>::value> &Q)
	{
		typedef Eigen::Matrix<scalar_type, dof<Measurement>::value, dof<Measurement>::value> measurement_cov;
		update(z, h,
			   boost::bind(id<measurement_cov>, Q),
			   accept_any_mahalanobis_distance<scalar_type>);
	}
	
	template<typename Measurement,
			 typename MeasurementModel,
			 typename MeasurementNoiseCovariance,
			 typename MahalanobisTest>
	void update(const Measurement &z,
				MeasurementModel h,
				MeasurementNoiseCovariance Q,
				MahalanobisTest mt)
	{
		const static int measurement_rows = dof<Measurement>::value;
		typedef Measurement measurement;
		typedef Eigen::Matrix<scalar_type, measurement_rows, 1> vectorized_measurement;
		typedef Matrix<scalar_type, measurement_rows, measurement_rows> measurement_cov;
		typedef Matrix<scalar_type, state::DOF, measurement_rows> cross_cov;

		state_vector X(2 * n + 1);
		generate_sigma_points(mu_, sigma_, X);

		std::vector<measurement> Z(X.size());
		std::transform(X.begin(), X.end(), Z.begin(), h);
		
		const measurement meanZ = sigma_points_mean(Z);
		const measurement_cov S = sigma_points_cov<measurement_rows>(meanZ, Z) + Q();
		const cross_cov covXZ = sigma_points_cross_cov<measurement_rows>(mu_, meanZ, X, Z);

		measurement_cov S_inverse(S.inverse());

		const cross_cov K = covXZ * S_inverse;
		
		const vectorized_measurement innovation = z - meanZ;

		const scalar_type mahalanobis2 = (innovation.transpose() * S_inverse * innovation)(0);

		if (mt(mahalanobis2))
		{
			sigma_ -= K * S * K.transpose();
			apply_delta(K * innovation);
		}
	}

	const state &mu() const
	{
		return mu_;
	}

	const cov &sigma() const
	{
		return sigma_;
	}
	
private:

	void generate_sigma_points(const state &mu,
							   const vectorized_state &delta,
							   const cov &sigma,
							   state_vector &X) const
	{
		assert(X.size() == 2 * n + 1);

		MatrixXd L(sigma.llt().matrixL());

		/*if (!L.isSPD())
		{
			std::cerr << std::endl << "sigma is not SPD:" << std::endl
					  << sigma << std::endl
					  << "---" << std::endl;
			Eigen::EigenSolver<cov> eig(sigma);
			std::cerr << "eigen values: " << eig.eigenvalues().transpose() << std::endl;
		}
		
		assert(L.isSPD());*/

		/*
		std::cout << ">> L" << std::endl
				  << L.getL() << std::endl
				  << "<< L" << std::endl;
		*/
		
		X[0] = mu + delta;
		for (std::size_t i = 1, j = 0; j < n; ++j)
		{
			//std::cout << "L.col(" << j << "): " << L.getL().col(j).transpose() << std::endl;
			X[i++] = mu + (delta + L.col(j));
			X[i++] = mu + (delta - L.col(j));
		}
		//print_sigma_points(X);
	}

	void generate_sigma_points(const state &mu,
							   const cov &sigma,
							   state_vector &X) const
	{
		generate_sigma_points(mu, vectorized_state::Zero(), sigma, X);
	}

	// manifold mean
	template<typename manifold>
	manifold
	sigma_points_mean(const std::vector<manifold> &X) const
	{
		manifold reference = X[0];
		typename manifold::vectorized_type mean_delta;
		const static std::size_t max_it = 10000;

		std::size_t i = 0;
		do {
			mean_delta.setZero();
			for (typename std::vector<manifold>::const_iterator Xi = X.begin(); Xi != X.end(); ++Xi)
			{
				mean_delta += *Xi - reference;
			}
			mean_delta /= X.size();
			reference += mean_delta;
		} while (mean_delta.norm() > 1e-6
				 && ++i < max_it);

		if (i >= max_it)
		{
			std::cerr << "ERROR: sigma_points_mean() did not converge. norm(mean_delta)=" << mean_delta.norm() << std::endl;
			assert(false);
		}
		
		return reference;
	}
	
	// vector mean
	template<int measurement_rows>
	Matrix<scalar_type, measurement_rows, 1>
	sigma_points_mean(const std::vector<Matrix<scalar_type, measurement_rows, 1> > &Z) const
	{
		typedef Matrix<scalar_type, measurement_rows, 1> measurement;
		
		return std::accumulate(Z.begin(), Z.end(), measurement(measurement::Zero())) / Z.size();
	}

#ifdef VECT_H_
	// MTK vector mean
	template<int measurement_rows>
	MTK::vect<measurement_rows, scalar_type>
	sigma_points_mean(const std::vector<MTK::vect<measurement_rows, scalar_type> > &Z) const
	{
		typedef MTK::vect<measurement_rows, scalar_type> measurement;
		
		return std::accumulate(Z.begin(), Z.end(), measurement(measurement::Zero())) / Z.size();
	}
#endif // VECT_H_
	
	template<int cov_size, typename T>
	Matrix<scalar_type, -1, -1>
	sigma_points_cov(const T &mean, const std::vector<T> &V) const
	{
		typedef Matrix<scalar_type, -1, -1> cov_mat;
		typedef Matrix<scalar_type, cov_size, 1> cov_col;
		
		cov_mat c(cov_mat::Zero(cov_size,cov_size));
		
		for (typename std::vector<T>::const_iterator Vi = V.begin(); Vi != V.end(); ++Vi)
		{
			cov_col d = *Vi - mean;
			c += d * d.transpose();
		}

		return 0.5 * c;
	}

	template<int measurement_rows, typename Measurement>
	Matrix<scalar_type, state::DOF, measurement_rows>
	sigma_points_cross_cov(const state &meanX,
						   const Measurement &meanZ,
						   const state_vector &X,
						   const std::vector<Measurement> &Z) const
	{
		assert(X.size() == Z.size());

		typedef Matrix<scalar_type, state::DOF, measurement_rows> cross_cov;

		cross_cov c(cross_cov::Zero());

		{
			typename state_vector::const_iterator Xi = X.begin();
			typename std::vector<Measurement>::const_iterator Zi = Z.begin();
			for (;Zi != Z.end(); ++Xi, ++Zi)
			{
				c += (*Xi - meanX) * (*Zi - meanZ).transpose();
			}
		}
		
		return 0.5 * c;
	}
	
	void apply_delta(const vectorized_state &delta)
	{
		state_vector X(2 * n + 1);
		generate_sigma_points(mu_, delta, sigma_, X);

		mu_ = sigma_points_mean(X);
		sigma_ = sigma_points_cov<state::DOF>(mu_, X);
	}

public: state mu_;
public: cov sigma_;

	// for debugging only

	void print_sigma_points(const state_vector &X) const
	{
		std::cout << "generated sigma points:" << std::endl;
		for (typename state_vector::const_iterator Xi = X.begin(); Xi != X.end(); ++Xi)
		{
			std::cout << *Xi << std::endl << "***" << std::endl;
		}
	}

public:
	void check_sigma_points()
    {
		state_vector X(2 * n + 1);
		generate_sigma_points(mu_, sigma_, X);

		state muX = sigma_points_mean(X);
		
		cov sigma_test = sigma_points_cov<state::DOF>(muX, X);
		if((sigma_test - sigma_).cwise().abs().maxCoeff()>1e-6){
			std::cerr << sigma_test << "\n\n" << sigma_;
			assert(false);
		}

		if (mu_ != muX)
		{
//			std::cout << "mu_:" << mu_ << std::endl;
//			std::cout << "muX:" << muX << std::endl;
			std::cout << "norm:" << ((mu_ - muX).norm() > 0. ? ">" : "=") << std::endl;
		}
		assert (mu_ == muX);
	}

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
	
} // namespace ukfom
#endif // __UKFOM_UKF_HPP__
