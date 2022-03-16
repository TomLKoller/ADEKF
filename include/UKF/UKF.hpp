/*
 * EKF.hpp
 *
 *  Created on: 27.07.2018
 *      Author: tomlucas
 */

#ifndef ESTIMATORS_UKF_HPP_
#define ESTIMATORS_UKF_HPP_

#define FILL_LATER 0

#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include "../ADEKF/ADEKFUtils.h"
#include <memory>

namespace adekf {
	/**
	 * @brief Generic UKF implementation
	 * 
	 * @tparam State state description (either Eigen::Vector or adekf::Manifold)
	 */
template<typename State>
class UKF {
protected:

public:


	UKF(bool prepareSmoothing = false) :
			state_count(0), alignment(alignment.Identity()), prepareSmoothing(prepareSmoothing) {
	}
	typedef State STATE_TYPE;     //model type
	static constexpr int DOF = adekf::DOFOf<State>;
	using ScalarType = typename adekf::StateInfo<State>::ScalarType;
	typedef Eigen::Matrix<ScalarType, DOF, 1> DELTA_TYPE;
	typedef adekf::CovarianceOf<STATE_TYPE> STATE_COV_TYPE;
	typedef Eigen::Matrix<ScalarType, DOF, DOF * 2 + 1> SIGMA_SIGMA_POINTS_TYPE;




	UKF(const STATE_TYPE & start_state, const STATE_COV_TYPE &start_cov, bool prepareSmoothing=false):UKF<State>(prepareSmoothing) {
			setStart(start_state,start_cov);
		}

	virtual ~UKF() {
		//smoothAllEstimates();

	}

	/**
	 * Returns the sigma points of the state
	 *
	 *
	 * @param state the state to get the sigma points off
	 * @param cov the covariance of the state
	 * @return a eigen matrix with the sigma points
	 */

	std::shared_ptr<STATE_TYPE[]> getSigmaPoints(const STATE_TYPE & state, const STATE_COV_TYPE & cov) {
		STATE_COV_TYPE cholesky = cov.llt().matrixL();
		std::shared_ptr<STATE_TYPE[]> sigma_points(new STATE_TYPE[2 * DOF + 1]);
		sigma_points[0] = state;
		for (int i = 0; i < DOF; i++) {
			sigma_points[2 * i + 1] = state +(cholesky.col(i)).eval();
			sigma_points[2 * i + 2] = state + (-cholesky.col(i)).eval();

		}
		return sigma_points;
	}

	/**
	 * Simple wrapper to call getSigmaPoints without arguments
	 * @return sigma points of state_vector and sigma
	 */
	virtual std::shared_ptr<STATE_TYPE[]> getSigmaPoints() {
		return getSigmaPoints(mu, sigma);
	}

	/**
	 * Calculates the mean of sigma points
	 * @param sigma_points  the matrix with all sigma points
	 * @param epsilon the stopping criteria
	 * @param max_iterations max iterations of convergence
	 * @return the mean of sigma points
	 */
	template<typename MEASURE_TYPE>
	MEASURE_TYPE meanOfSigmaPoints(const std::shared_ptr<MEASURE_TYPE[]> sigma_points, double epsilon = 1e-6,
			int max_iterations = 30) {

		MEASURE_TYPE mean = sigma_points[0];
		MEASURE_TYPE old_mean = sigma_points[0];
		Eigen::Matrix<ScalarType,adekf::DOFOf<MEASURE_TYPE>,1> diff_sum = diff_sum.Zero();
		int iterations = 0;
		do {
			iterations++;
			old_mean = mean;
			diff_sum = diff_sum.Zero();
			for (int i = 0; i <= DOF * 2; ++i) {
				diff_sum += sigma_points[i] - mean;
			}
			diff_sum /= 2. * DOF + 1.;
			mean = mean + diff_sum;
		} while (iterations <= max_iterations && diff_sum.norm() > epsilon);
		if (iterations > max_iterations)
			printf("Warning: stopped due to excess of iterations");
		return mean;
	}
	/**
	 * Saves all releveant states for smoothing
	 * @param input the input u
	 * @param time_diff time since last call
	 */
	void saveStatesForSmoothing(double time_diff) {
		state_count++;
		past_states.push_back(mu);
		past_states_smoothed.push_back(mu);
		past_covs.push_back(sigma);
	}
	/**
	 * Does  dynamic step in EKF
	 *
	 * @param input the input u
	 * @param time_diff time since last call
	 */
	template<typename DynamicModel, typename ... Controls>
	void predict(DynamicModel dynamicModel, const STATE_COV_TYPE &processNoise, const Controls ...u) {
		auto f = std::bind(dynamicModel, std::placeholders::_1, u...);
		if (prepareSmoothing) {
			dynamic_functions.push_back(FunctionWrapperForDynamics<STATE_TYPE>(f));
			past_process_noises.push_back(processNoise);
		}
		auto sigma_points = getSigmaPoints();
		for (int i = 0; i <= DOF * 2; i++) {
			f(sigma_points[i]);
		}
		mu = meanOfSigmaPoints(sigma_points);
		//printf(state_vector);
		//printf(" ");
		SIGMA_SIGMA_POINTS_TYPE result = SIGMA_SIGMA_POINTS_TYPE::Zero();
		for (int i = 0; i <= DOF * 2; i++) {
			result.col(i) = sigma_points[i] - mu;
		}
		sigma = 0.5 * (result * result.transpose()) + processNoise;
		//printf("dynamic");
		//printf(state_vector.block(3, 0, 3, 1).norm());
	}
	/**
	 * Does  a  measurement update in UKF
	 * @param measurement the measurement
	 * @param time_diff time since last call
	 * @param measure_function function to map a state to a predicted measurement
	 * @param noise_function gives the measurement noise with time_diff
	 * @param boxplus_m boxplus for measurement
	 * @param boxminus_m boxminus for measurement
	 * measure dim is the dimension of the measurement vector
	 */
	template<typename Measurement, typename MeasurementModel, typename ... Variables>
	void update(MeasurementModel measurementModel, const adekf::CovarianceOf<Measurement> &measurementNoise,
			const Measurement &measurement, const Variables ...variables) {
		auto sigma_points = getSigmaPoints();
		constexpr int measure_dim = adekf::DOFOf<Measurement>;
		auto h = std::bind(measurementModel, std::placeholders::_1, variables...);
		std::shared_ptr<Measurement[]> expected_zs(new Measurement[2*DOF+1]);
		for (int i = 0; i <= 2 * DOF; ++i) {
			expected_zs[i] = h(sigma_points[i]);
		}
		auto mean_z = meanOfSigmaPoints(expected_zs);     //< expected measurement mean
		Eigen::Matrix<double, measure_dim, DOF * 2 + 1> diff_z;
		for (int i = 0; i <= DOF * 2; ++i) {
			diff_z.col(i) = expected_zs[i] - mean_z;
		}
		Eigen::Matrix<double, measure_dim, measure_dim> sigma_z = 0.5 * (diff_z * diff_z.transpose())
				+ measurementNoise;
		SIGMA_SIGMA_POINTS_TYPE result = SIGMA_SIGMA_POINTS_TYPE::Zero();
		for (int i = 0; i <= DOF * 2; ++i) {
			result.col(i) = sigma_points[i] - mu;
		}
		Eigen::Matrix<double, DOF, measure_dim> sigma_xz = 0.5 * (result * diff_z.transpose());
		Eigen::Matrix<double, DOF, measure_dim> kalman_gain = sigma_xz * sigma_z.inverse();
		Eigen::Matrix<double, DOF, 1> delta = kalman_gain * (measurement - mean_z);
		STATE_COV_TYPE sigma_t = sigma - (kalman_gain * sigma_z * kalman_gain.transpose());

		STATE_COV_TYPE cholesky = sigma_t.llt().matrixL();
		std::shared_ptr<STATE_TYPE[]> sigma_points_second(new STATE_TYPE[2 * DOF + 1]);
		sigma_points_second[0] = mu + delta;
		for (int i = 0; i < DOF; i++) {
			sigma_points_second[2 * i + 1] = mu + (delta + cholesky.col(i)).eval();
			sigma_points_second[2 * i + 2] = mu + (delta - cholesky.col(i)).eval();
		}

		mu = meanOfSigmaPoints(sigma_points_second);
		for (int i = 0; i <= DOF * 2; ++i) {
			result.col(i) = sigma_points_second[i] - mu;
		}
		sigma = 0.5 * (result * result.transpose());
	}
	/**
	 * Does  a  measurement update in UKF
	 * @param measurement the measurement
	 * @param time_diff time since last call
	 * @param measure_function function to map a state to a predicted measurement
	 * @param noise_function gives the measurement noise with time_diff
	 * measure dim is the dimension of the measurement vector
	 */
	template<typename MeasurementDerived, typename MeasurementModel, typename ... Variables>
	void update(MeasurementModel measurementModel,
			const adekf::CovarianceOf<Eigen::MatrixBase<MeasurementDerived>> &measurementNoise,
			const Eigen::MatrixBase<MeasurementDerived> &measurement, const Variables ...variables) {
		auto h = std::bind(measurementModel, std::placeholders::_1, variables...);
		constexpr int measure_dim = adekf::DOFOf<Eigen::MatrixBase<MeasurementDerived>>;
		auto sigma_points = getSigmaPoints();
		Eigen::Matrix<double, measure_dim, DOF * 2 + 1> expected_zs;
		for (int i = DOF * 2; i >= 0; --i) {
			expected_zs.col(i) = h(sigma_points[i]);
		}
		Eigen::Matrix<double, measure_dim, 1> mean_z = expected_zs.rowwise().mean();     //< expected measurement mean
		Eigen::Matrix<double, measure_dim, DOF * 2 + 1> diff_z = expected_zs.colwise() - mean_z;
		Eigen::Matrix<double, measure_dim, measure_dim> sigma_z = 0.5 * (diff_z * diff_z.transpose())
				+ measurementNoise;

		SIGMA_SIGMA_POINTS_TYPE result = SIGMA_SIGMA_POINTS_TYPE::Zero();
		for (int i = 0; i <= DOF * 2; ++i) {
			result.col(i) = sigma_points[i] - mu;
		}
		Eigen::Matrix<double, DOF, measure_dim> sigma_xz = 0.5 * (result * diff_z.transpose());

		Eigen::Matrix<double, DOF, measure_dim> kalman_gain = sigma_xz * sigma_z.inverse();

		mu = mu + (kalman_gain * (measurement - mean_z));

		sigma = sigma - kalman_gain * sigma_xz.transpose();
	}
	/**
	 * wrapper to call measurementstep with different arguments
	 */
	template<int measure_dim>
	struct MeasurementWrapper {
		Eigen::Matrix<double, measure_dim, 1> (*function)(const STATE_TYPE &, void *);
		MeasurementWrapper(Eigen::Matrix<double, measure_dim, 1> (*function)(const STATE_TYPE &, void *)) :
				function(function) {

		}
		Eigen::Matrix<double, measure_dim, 1> operator()(const STATE_TYPE & state,
				const Eigen::Matrix<double, 4, 4> & alignment, void *prior) const {
			return function(state, prior);
		}

	};
	template<int measure_dim>
	void measurementStep(const Eigen::Matrix<double, measure_dim, 1> & measurement, double time_diff,
			Eigen::Matrix<double, measure_dim, 1> (*measure_function)(const STATE_TYPE &, void * prior),
			const Eigen::Matrix<double, measure_dim, measure_dim> & noise, void * prior = NULL) {
		measurementStep(measurement, time_diff, MeasurementWrapper<measure_dim>(measure_function), noise, prior);
	}

	/**
	 * Smooth previous estimates with new knowledge
	 *
	 * This is taken from:
	 *
	 * Unscented Rauch–Tung–Striebel Smoother by Simo S"arkk"a
	 *
	 * @param k_length the amount of steps to smooth back
	 * @param k_start the starting index for the smoothing
	 */
	void smoothEstimates(unsigned int k_length, unsigned int k_start) {
		if (k_start > state_count - 1) {
			ERR_STREAM<<"Trying to smooth from a non existent state" LOG_END;
			return;
		}

		if (k_length > k_start) {
			ERR_STREAM<< "Trying to smooth more states than are before k_start" LOG_END;
			return;
		}
		for (unsigned int k = k_start; k > k_start - k_length; k--) {
			auto sigma_points = getSigmaPoints(past_states_smoothed[k], past_covs[k]);
			decltype(sigma_points) sigma_points_plus;
			for (int i = 0; i < DOF * 2 + 1; i++) {
				sigma_points_plus[i] = dynamic_functions[k](sigma_points.col(i));
			}

			STATE_TYPE mean = meanOfSigmaPoints(sigma_points_plus);
			SIGMA_SIGMA_POINTS_TYPE result_plus;
			for (int i = 0; i < DOF * 2 + 1; i++) {
				result_plus.col(i) = sigma_points_plus[i]-mean;
			}

			SIGMA_SIGMA_POINTS_TYPE result;
			STATE_COV_TYPE cov = 0.5*(result_plus * result_plus.transpose())+ past_process_noises[k];
			for (int i = 0; i < DOF * 2 + 1; i++) {
				result.col(i) = boxMinus(past_states_smoothed[k], sigma_points.col(i));
			}
			STATE_COV_TYPE c_k_plus =0.5* result * result_plus.transpose();
			STATE_COV_TYPE d_k = c_k_plus * cov.inverse();

			//from here the second sigma propagation applies
			DELTA_TYPE delta=d_k * (past_states_smoothed[k + 1]-mean);
			STATE_COV_TYPE pst_cov_k=past_covs[k] + d_k * (past_covs[k + 1] - cov) * d_k.transpose();

			STATE_COV_TYPE cholesky = pst_cov_k.llt().matrixL();
			std::shared_ptr<STATE_TYPE[]> sigma_points_second(new STATE_TYPE[2 * DOF + 1]);
			sigma_points_second[0] = past_states_smoothed[k] + delta;
			for (int i = 0; i < DOF; i++) {
				sigma_points_second[2 * i + 1] = past_states_smoothed[k] + (delta + cholesky.col(i));
				sigma_points_second[2 * i + 2] = past_states_smoothed[k] + (delta - cholesky.col(i));
			}

			past_states_smoothed[k] = meanOfSigmaPoints(sigma_points_second);
			for (int i = 0; i <= DOF * 2; ++i) {
				result.col(i) = boxMinus(past_states_smoothed[k], sigma_points_second.col(i));
			}
			past_covs[k] = 0.5 * (result * result.transpose());
		}
	}
	/**
	 * Perform smoothing from the newest estimate
	 * @param k_length the amount of estimates to smooth
	 */
	void smoothEstimates(unsigned int k_length) {

		smoothEstimates(k_length, state_count - 2);

	}
	/**
	 * Perform smoothing from the newest estimate for all estimates
	 */
	void smoothAllEstimates() {
		smoothEstimates(state_count - 2);
	}

	/**
	 * Sets the start estimat \hat(x) (0)
	 * @param start_state the starting state vector
	 * @param start_cov  the starting covariance
	 */
	virtual inline void setStart(const STATE_TYPE & start_state, const STATE_COV_TYPE &start_cov) {
		mu = start_state;
		sigma = start_cov;
	}

	inline STATE_TYPE & getStateVector() {
		return mu;
	}
	inline STATE_COV_TYPE getCov() {
		return sigma;
	}

	std::vector<STATE_TYPE> & getSmoothedStates() {
		return past_states_smoothed;
	}

	void setAlignment(const Eigen::Matrix4d & alignment) {
		this->alignment = alignment;
	}

	template<typename pass_type>
	class FunctionWrapperForDynamics {
	private:

		struct TypeErasure {
			virtual ~TypeErasure() {
			}
			virtual void operator()(pass_type & argument) const =0;

		};
		template<typename functor>
		struct FunctionWrapper: public TypeErasure {
			FunctionWrapper(const functor &function) :
			function(function) {
			}

			virtual void operator()(pass_type & argument) const {

				function(argument);
			}
			functor function;
		};
	public:
		/**
		 * Calls the stored functor with parameters
		 * @param time time_point on wich the function es evaluated
		 * @param expand_noise whether to add a noisy increment
		 * @return a Eigen 3D Vector as the function result
		 */
		void operator()(pass_type & argument) const {
			(*function)(argument);
		}
		/**
		 * Create the funcCaller from an arbitrary object which implements the operator()((double * time, double * result, bool expand_noise)
		 * @param func function object
		 */
		template<typename functor>
		FunctionWrapperForDynamics(const functor & func) :
		function(new FunctionWrapper<functor>(func)) {
		}
		/**
		 * Copy constructor
		 * @param func another FuncCaller
		 */
		FunctionWrapperForDynamics(FunctionWrapperForDynamics & func) :
		function(func.function) {
		}

	protected:
		std::shared_ptr<TypeErasure> function;
	};
	STATE_TYPE mu;     //< the current estimated state x
	STATE_COV_TYPE sigma;//< the estimated stiffness(x)
	
protected:
	std::vector<STATE_TYPE> past_states;//< all past states
	std::vector<STATE_TYPE> past_states_smoothed;//< all past states
	//std::vector<INPUT_TYPE> past_inputs;//< all past inputs
	std::vector<FunctionWrapperForDynamics<STATE_TYPE>> dynamic_functions;
	std::vector<STATE_COV_TYPE> past_process_noises;
	std::vector<STATE_COV_TYPE> past_covs;// < all past covariance matrices
	unsigned int state_count;//< the current state index
	Eigen::Matrix<double, 4, 4> alignment;
	bool prepareSmoothing;// tells the ukf whether it has to store dynamics
};

template<typename Derived,typename deri, typename WHATEVER>
	UKF(const Eigen::CwiseNullaryOp<Derived,deri> & start_state, const WHATEVER &start_cov, bool prepareSmoothing=false) ->UKF< Eigen::Vector3d >;

}
/* namespace adekf */

#endif /* ESTIMATORS_UKF_HPP_ */
