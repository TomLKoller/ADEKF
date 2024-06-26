#pragma once

#include <Eigen/LU>

#include "ceres/jet.h"
#include "ADEKFUtils.h"

#include <iostream>


namespace adekf {


    using namespace Eigen;
    using namespace std::placeholders;

    /**
     * An EKF Implementation for automatic differentiation of Jacobian Matrices
     * @tparam State The State to be used for estimation
     */
    template<typename State>
    class ADEKF {
        /**
         * The DOF of the State
         */
        static constexpr int DOF = DOFOf<State>;
        protected:
        /**
         * The Type of the Scalars used in the State
         */
        using ScalarType = typename StateInfo<State>::ScalarType;

        /**
         * A Matrix with the Scalar Type of the State
         * @tparam N Number of Rows
         * @tparam M Number of Columns
         */
        template<int N, int M>
        using MatrixType = Matrix<ScalarType, N, M>;

        /**
         * A Square Matrix with the Scalar Type of the State
         * @tparam N Number of Rows and Columns
         */
        template<int N>
        using SquareMatrixType = MatrixType<N, N>;




        /**
         * The Jacobian of a given State or Measurement Class
         * @tparam The class to be used
         */
        template<typename T>
        using JacobianOf = typename std::conditional<dynamicMatrix<DOFOf<T>, DOF>, MatrixType<-1, -1>, MatrixType<DOFOf<T>, DOF>>::type;

        /**
         * A Vector of dual components
         * @tparam N The dimension of the Vector and dual components
         */
        template<int N>
        using Derivator = Matrix<ceres::Jet<ScalarType, N>, N, 1>;

        static_assert(DOF > 0, "Only Fixed Size States and Manifolds are supported");

    public:
        /**
         * The Covariance type of the State
         */
        using Covariance = typename std::conditional<dynamicMatrix<DOF, DOF>, MatrixType<-1, -1>, SquareMatrixType<DOF>>::type;

        /**
         * The Expected Value of the State
         */
        State mu;

        /**
         * The Covariance of the State
         */
        Covariance sigma;

        /**
         * Constructor of the ADEKF
         * @param _mu Initial Expected Value of the State
         * @param _sigma Initial Covariance of the State
         */
        ADEKF(const State &_mu, const Covariance &_sigma) : mu(_mu), sigma(_sigma) {}


        /**
         * Predict the State Estimate with automatically differentiated Jacobian Matrices
         * @tparam DynamicModel Type of the Dynamic Model Functor
         * @tparam Controls Types of the Control Vectors
         * @param dynamicModel The Dynamic Model f(x,u)
         * @param Q Additive Process Noise Covariance
         * @param u Control Vectors
         */
        template<typename DynamicModel, typename... Controls>
        void predict(DynamicModel dynamicModel, const Covariance &Q, const Controls &...u) {
            //The Jacobian to be calculated from the dynamic Model
            JacobianOf<State> F(DOF, DOF);
            //Bind the control vectors to the dynamic Model
            auto f = std::bind(dynamicModel, _1, u...);
            //Add a dual component vector to the state
            auto input = eval(mu + getDerivator<DOF>());
            //Evaluate the dynamic model
            f(input);
            //Calculate the Jacobian Matrix and set the new State Estimate
            predict_impl(input, f, input, F);
            //The dynamic model has to be differentiable
            assert(!F.hasNaN() && "Differentiation resulted in an indeterminate form");
            //Calculate the new Covariance
            sigma = F * sigma * F.transpose() + Q;
        }


        
        /**
         * Predict the State Estimate with automatically differentiated Jacobian Matrices and non additive Noise
         * @tparam NoiseDim The Dimension of the Noise Vector w
         * @tparam DynamicModel Type of the Dynamic Model Functor
         * @tparam Controls Types of the Control Vectors
         * @param dynamicModel The Dynamic Model f(x,w,u)
         * @param Q non-Additive Process Noise Covariance
         * @param u Control Vectors
         */
        template<int NoiseDim, typename DynamicModel, typename... Controls>
        void predictWithNonAdditiveNoise(DynamicModel dynamicModel, const SquareMatrixType<NoiseDim> &Q,
                                         const Controls &...u) {
            //The Jacobian to be calculated from the dynamic Model
            MatrixType<DOF, DOF + NoiseDim> F(DOF, DOF + NoiseDim);
            //Bind the control vectors to the dynamic Model
            auto f = std::bind(dynamicModel, _1, _2, u...);
            //Generate a Vector of dual components for the state and noise vector
            auto derivator = getDerivator<DOF + NoiseDim>();
            //Add the first DOF cells of the dual component vector to the State
            auto input = eval(mu + derivator.template head<DOF>());
            //Evaluate the dynamic model with the NoiseDim last cells of the dual component vector
            f(input, derivator.template tail<NoiseDim>());
            //Calculate the Jacobian Matrix and set the new State Estimate
            //The noise vector gets set to zero for this
            predict_impl(input, std::bind(f, _1, MatrixType<NoiseDim, 1>::Zero()), input, F);
            //The dynamic model has to be differentiable
            assert(!F.hasNaN() && "Differentiation resulted in an indeterminate form");
            //Calculate the new Covariance
            sigma = F.template leftCols<DOF>() * sigma * F.template leftCols<DOF>().transpose() +
                    F.template rightCols<NoiseDim>() * Q * F.template rightCols<NoiseDim>().transpose();
         }



        /**
         * Predict the State Estimate
         * @tparam DynamicModel Type of the Dynamic Model Functor
         * @tparam JacobianFunc Type of the Jacobian of the Dynamic Model
         * @tparam Controls Types of the Control Vectors
         * @param f The Dynamic Model f(x,u)
         * @param jacobianFunc The Jacobian of the Dynamic Model df/dx(x,u)
         * @param Q Additive Process Noise Covariance
         * @param u Control Vectors
         */
        template<typename DynamicModel, typename JacobianFunc, typename... Controls>
        void predictWithJacobian(DynamicModel f, JacobianFunc jacobianFunc, const Covariance &Q, const Controls &...u) {
            //The Jacobian, calculated from the given function
            auto F = jacobianFunc(mu, u...);
            //Evaluate the dynamic model and set the new state estimate
            f(mu, u...);
            //Calculate the new covariance
            sigma = F * sigma * F.transpose() + Q;
        }

        /**
         * Update the State Estimate with automatically differentiated Jacobian Matrices
         * @tparam Measurement Type of the Measurement
         * @tparam MeasurementModel Type of the Measurement Model Functor
         * @tparam Variables Types of Auxiliary Variables for the Measurement Model
         * @param measurementModel The Measurement Model h(x,variables)
         * @param R Additive Measurement Noise Covariance
         * @param z Measurement
         * @param variables Auxiliary Variables for the Measurement Model
         */
        template<typename Measurement, typename MeasurementModel, typename Derived, typename... Variables>
        void update(MeasurementModel measurementModel, const MatrixBase<Derived> &R, const Measurement &z,
                    const Variables &...variables) {
            //Bind the auxiliary variables to the measurement model
            auto h = [&measurementModel, &variables...](const auto &state) {
                return eval(measurementModel(state, variables ...));
            };
            //The jacobian matrix to be calculated from the measurement model
            JacobianOf<Measurement> H(DOFOf<Measurement>, DOF);
            //The result of the measurement model
            typename StateInfo<Measurement>::type hx;
            //The result of the measurement model with a dual component vector added to the state
            auto input = h(eval(mu + getDerivator<DOF>()));
            //Calculate the Jacobian and the result of the measurement model
            update_impl(hx, input, h, H);
            //The measurement model has to be differentiable
            assert(!H.hasNaN() && "Differentiation resulted in an indeterminate form");
            //Calculate the Innovation covariance
            auto S = H * sigma * H.transpose() + R;
            //Calcualte the Kalman Gain
            auto K = (sigma * H.transpose() * S.inverse()).eval();
            //Calculate the updated state estimate
            auto delta=eval(z-hx);
            //Calcualte the updated covariance estimate
            add_diff(mu, K * delta, K * H);

        }
        /**
         * Transpose overload to handle likelihood of scalar updates
         */
        template<typename Type>
        auto transpose(const Type & object){
            return object.transpose();
        }
        double transpose(double value){
            return value;
        }

        template<typename Derived>
        double exp(const Eigen::MatrixBase<Derived> & matrix){
            assert(matrix.rows()==1 && matrix.cols()==1);
            return std::exp(matrix(0,0));
        }

        /**
        * Update the State Estimate with automatically differentiated Jacobian Matrices
        * @tparam Measurement Type of the Measurement
        * @tparam MeasurementModel Type of the Measurement Model Functor
        * @tparam Variables Types of Auxiliary Variables for the Measurement Model
        * @param log_likelihood output: the likelihood of the given measurement
        * @param measurementModel The Measurement Model h(x,variables)
        * @param R Additive Measurement Noise Covariance
        * @param z Measurement
        * @param variables Auxiliary Variables for the Measurement Model
        */
        template<typename Measurement, typename MeasurementModel, typename Derived, typename... Variables>
        void update(double & log_likelihood, MeasurementModel measurementModel, const MatrixBase<Derived> &R, const Measurement &z,
                    const Variables &...variables) {
            //Bind the auxiliary variables to the measurement model
            auto h = [&measurementModel, &variables...](const auto &state) {
                return eval(measurementModel(state, variables ...));
            };
            //The jacobian matrix to be calculated from the measurement model
            JacobianOf<Measurement> H(DOFOf<Measurement>, DOF);
            //The result of the measurement model
            typename StateInfo<Measurement>::type hx;
            //The result of the measurement model with a dual component vector added to the state
            auto input = h(eval(mu + getDerivator<DOF>()));
            //Calculate the Jacobian and the result of the measurement model
            update_impl(hx, input, h, H);
            //The measurement model has to be differentiable
            assert(!H.hasNaN() && "Differentiation resulted in an indeterminate form");
            //Calculate the Innovation covariance
            auto S = H * sigma * H.transpose() + R;
            //Calculate the Kalman Gain
            auto K = (sigma * H.transpose() * S.inverse()).eval();
            //Calculate the updated state estimate
            auto delta=eval(z-hx);
            //std::cout << delta << std::endl;
            log_likelihood= -0.5 * (transpose(delta) * S.inverse() * delta)(0)+log(1/sqrt(S.determinant() * pow((2 * M_PI), S.rows())));
            //Calcualte the updated covariance estimate
            add_diff(mu, K * delta, K * H);

        }



        /**
         * Update the State Estimate with automatically differentiated Jacobian Matrices and non-Additive Noise
         * @tparam NoiseDim The Dimension of the Noise Vector v
         * @tparam Measurement Type of the Measurement
         * @tparam MeasurementModel Type of the Measurement Model Functor
         * @tparam Variables Types of Auxiliary Variables for the Measurement Model
         * @param measurementModel The Measurement Model h(x,v,variables)
         * @param R non-Additive Measurement Noise Covariance
         * @param z Measurement
         * @param variables Auxiliary Variables for the Measurement Model
         */
        template<int NoiseDim, typename Measurement, typename MeasurementModel, typename... Variables>
        void updateWithNonAdditiveNoise(MeasurementModel measurementModel, const SquareMatrixType<NoiseDim> &R,
                                        const Measurement &z, const Variables &...variables) {
            //The DOF of the Measurement
            constexpr int MDOF = DOFOf<Measurement>;
            //Bind the auxiliary variables to the measurement model
            auto h = [&measurementModel, &variables...](const auto &state, const auto &noise) {
                return eval(measurementModel(state,noise , variables ...));
            };
            //The jacobian matrix to be calculated from the measurement model
            MatrixType<MDOF, DOF + NoiseDim> H(MDOF, DOF);
            //The result of the measurement model
            typename StateInfo<Measurement>::type hx;
            //Generate a Vector of dual components for the state and noise vector
            auto derivator = getDerivator<DOF + NoiseDim>();
            //The result of the measurement model with a dual component vector added to the state and the noise
            auto input = h(eval(mu + derivator.template head<DOF>()), derivator.template tail<NoiseDim>());
            //Calculate the Jacobian and the result of the measurement model
            //The noise vector gets set to zero for this
            update_impl(hx, input, std::bind(h, _1, MatrixType<NoiseDim, 1>::Zero()), H);
            //The measurement model has to be differentiable
            assert(!H.hasNaN() && "Differentiation resulted in an indeterminate form");
            //Calculate the Innovation covariance
            auto S = H.template leftCols<DOF>() * sigma * H.template leftCols<DOF>().transpose() +
                     H.template rightCols<NoiseDim>() * R * H.template rightCols<NoiseDim>().transpose();
            //Calcualte the Kalman Gain
            auto K = (sigma * H.template leftCols<DOF>().transpose() * S.inverse()).eval();
            //Calculate the updated state estimate
            //Calcualte the updated covariance estimate
            add_diff(mu, K * (z - hx), K * H.template leftCols<DOF>());
        }

        /**
         * Update the State Estimate
         * @tparam Measurement Type of the Measurement
         * @tparam MeasurementModel Type of the Measurement Model Functor
         * @tparam JacobianFunc Type of the Jacobian of the Measurement Model
         * @tparam Variables Types of Auxiliary Variables for the Measurement Model
         * @param h The Measurement Model h(x,variables)
         * @param jacobianFunc The Jacobian of the Measurement Model dh/dx(x,variables)
         * @param R Additive Measurement Noise Covariance
         * @param z Measurement
         * @param variables Auxiliary Variables for the Measurement Model
         */
        template<typename Measurement, typename MeasurementModel, typename JacobianFunc, typename Derived, typename... Variables>
        void updateWithJacobian(MeasurementModel h, JacobianFunc jacobianFunc, const MatrixBase<Derived> &R,
                                const Measurement &z, const Variables &...variables) {
            //The jacobian matrix, calculated from the given function
            auto H = jacobianFunc(mu, variables...);
            //Calculate the Innovation covariance
            auto S = H * sigma * H.transpose() + R;
            //Calcualte the Kalman Gain
            auto K = (sigma * H.transpose() * S.inverse()).eval();
            //Calculate the updated state estimate
            //Calcualte the updated covariance estimate
            add_diff(mu, K * (z - h(mu, variables...)), K * H);
        }


/**
         * Update the State Estimate for manifolds and given jacobians
         * @tparam Measurement Type of the Measurement
         * @tparam MeasurementModel Type of the Measurement Model Functor
         * @tparam JacobianFunc Type of the Jacobian of the Measurement Model
         * @tparam JacobianFuncBoxPlus Type of the Jacobian of the boxplus operation
         * @tparam Variables Types of Auxiliary Variables for the Measurement Model
         * @param h The Measurement Model h(x,variables)
         * @param jacobianFunc The Jacobian of the Measurement Model dh/dx(x,variables)
         * @param jacobianFuncBoxPlus The Jacobian of the Boxplus operation ((mu boxplus delta) boxplus diff) boxminus (mu boxplus diff)
         * @param R Additive Measurement Noise Covariance
         * @param z Measurement
         * @param variables Auxiliary Variables for the Measurement Model
         */
        template<typename Measurement, typename MeasurementModel, typename JacobianFunc, typename JacobianFuncBoxPlus, typename Derived, typename... Variables>
        void updateManifoldWithJacobian(MeasurementModel h, JacobianFunc jacobianFunc,
                                        JacobianFuncBoxPlus jacobianFuncBoxPlus, const MatrixBase<Derived> &R,
                                        const Measurement &z, const Variables &...variables) {
            //The jacobian matrix, calculated from the given function
            auto H = jacobianFunc(mu, variables...);
            //Calculate the Innovation covariance
            auto S = H * sigma * H.transpose() + R;
            //Calcualte the Kalman Gain
            auto K = (sigma * H.transpose() * S.inverse()).eval();
            auto y = (K * (z - h(mu, variables...))).eval();
            //Calculate the updated state estimate
            State newMu = mu + y;
            //Definition of the Jacobian of the Boxplus Function
            JacobianOf<State> D = jacobianFuncBoxPlus(mu, y, variables...);
            //Set the new Estimated Value
            mu = newMu;
            //Calculate the new Covariance Matrix
            sigma = D * (sigma - (K * H * sigma)) * D.transpose();

        }


    private:

       


        /**
         * Calculation of the new Expeted Value and Jacobian during Prediction with a Manifold as State
         * @tparam DOFandNoise The DOF of the State plus the dimension of a possible Noise Vector
         * @tparam DynamicModel The Type of the Dynamic Model Functor
         * @tparam Manifold The Type of Manifold used as State
         * @param f The Dynamic Model f(x)
         * @param input Result of an Addition of the State and a dual component
         * @param F The resulting Jacobian. Calculated with dual numbers
         */
        template<typename Derived, typename DynamicModel, typename ManifoldType>
        void predict_impl(const Manifold &, DynamicModel f, const ManifoldType &input, MatrixBase<Derived> &F) {
            //Set the new state estimate
            f(mu);
            //The difference of a differentiated manifold with it's identity results in the jacobian
            extractJacobian(input,mu,F);
        }


        /**
       * Calculation of the new Expeted Value and Jacobian during Prediction with a CompoundManifold as State
       * @tparam Derived The MatrixType of the Covariance
       * @tparam DynamicModel The Type of the Dynamic Model Functor
       * @tparam ManifoldType The Type of Manifold used as State
       * @param f The Dynamic Model f(x)
       * @param input Result of an Addition of the State and a dual component
       * @param F The resulting Jacobian. Calculated with dual numbers
       */
        template<typename Derived, typename DynamicModel, typename ManifoldType>
        void predict_impl(const CompoundManifold &, DynamicModel f, const ManifoldType &input, MatrixBase<Derived> &F) {
            //check if state is a vector compound manifold
            if constexpr(mu.MAN_DOF==0){
                for (int i = 0; i < DOF; ++i) {
                    F.row(i) = input.vector_part(i).v;
                    mu.vector_part(i) = input.vector_part(i).a;
                }
                return;
            }
            //Set the new state estimate
            f(mu);
            //calculate the Jacobian
            extractJacobian(input, mu, F);
        }


        /**
         * Calculation of the new Expeted Value and Jacobian during Prediction with a Matrix as State
         * @tparam DOFandNoise The DOF of the State plus the dimension of a possible Noise Vector
         * @tparam DynamicModel The Type of the Dynamic Model Functor (for interface purposes)
         * @param input Result of an Addition of the State and a dual component
         * @param F The resulting Jacobian. Calculated with dual numbers
         */
        template<typename Derived, typename Derived2, typename DynamicModel>
        void predict_impl(const MatrixBase<Derived> &, DynamicModel, const MatrixBase<Derived> &input,
                          MatrixBase<Derived2> &F) {
            //The real component of the dual numbers are the result of the dynamic model
            //The dual component vectors represent the rows of the jacobian matrix
            for (int i = 0; i < DOF; ++i) {
                F.row(i) = input[i].v;
                mu[i] = input[i].a;
            }
        }


        /**
         * Pass to real update_impl_
          * @tparam Measurement The Type of Manifold used as Measurement
         * @tparam ModelReturn The Type of the Addition between a Measurement and a dual component
         * @tparam MeasurementModel The Type of the Measurement Model Functor
         * @tparam Derived Type of the Covaraince Matrix
         * @param modelResult Result of the Measurement Model
         * @param input Result of the Measurement Model with a dual compnent added to the state
         * @param h The Measurement Model h(x)
         * @param H The resulting Jacobian. Calculated with dual numbers
         */template<typename Measurement, typename ModelReturn, typename MeasurementModel, typename Derived>
        void
        update_impl(Measurement &modelResult, const ModelReturn &input, MeasurementModel h, MatrixBase<Derived> &H) {
            update_impl_(modelResult, modelResult, input, h, H);
        }

        /**
         * Calculation of the Observation and Jacobian during Update with a CompoundManifold as Measurement
         * @tparam Measurement The Type of Manifold used as Measurement
         * @tparam ModelReturn The Type of the Addition between a Measurement and a dual component
         * @tparam MeasurementModel The Type of the Measurement Model Functor
         * @tparam Derived Type of the Covaraince Matrix
         * @param modelResult Result of the Measurement Model
         * @param input Result of the Measurement Model with a dual compnent added to the state
         * @param h The Measurement Model h(x)
         * @param H The resulting Jacobian. Calculated with dual numbers
         */
        template<typename Measurement, typename ModelReturn, typename MeasurementModel, typename Derived>
        void
        update_impl_(const CompoundManifold &, Measurement &modelResult, const ModelReturn &input, MeasurementModel h,
                     MatrixBase<Derived> &H) {
            //Set the observation to the result of the measurement model
            modelResult = h(mu);
            //calculate the Jacobian
            extractJacobian(input, modelResult, H);
        }

        /**
         * Calculation of the Observation and Jacobian during Update with a Manifold as Measurement
         * @tparam Measurement The Type of Manifold used as Measurement
         * @tparam ModelReturn The Type of the Addition between a Measurement and a dual component
         * @tparam MeasurementModel The Type of the Measurement Model Functor
         * @tparam Derived Type of the Covaraince Matrix
         * @param modelResult Result of the Measurement Model
         * @param input Result of the Measurement Model with a dual compnent added to the state
         * @param h The Measurement Model h(x)
         * @param H The resulting Jacobian. Calculated with dual numbers
         */
        template<typename Measurement, typename ModelReturn, typename MeasurementModel, typename Derived>
        void
        update_impl_(const Manifold &, Measurement &modelResult, const ModelReturn &input, MeasurementModel h,
                     MatrixBase<Derived> &H) {
            //Set the observation to the result of the measurement model
            modelResult = h(mu);
            //calculate the Jacobian
            extractJacobian(input,modelResult,H);

        }


        /**
         * Calculation of the Observation and Jacobian during Update with a Matrix as Measurement
         * @tparam MDOF The DOF of the Measurement
         * @tparam DOFandNoise The DOF of the State plus the dimension of a possible Noise Vector
         * @tparam ModelReturn The Type of the Addition between a Measurement and a dual component
         * @tparam MeasurementModel The Type of the Measurement Model Functor (for interface purposes)
         * @param modelResult Result of the Measurement Model
         * @param input Result of the Measurement Model with a dual compnent added to the state
         * @param H The resulting Jacobian. Calculated with dual numbers
         */
        template<typename DerivedMeasurement, typename ModelReturn, typename MeasurementModel, typename Derived>
        void update_impl_(const MatrixBase<DerivedMeasurement> &, MatrixBase<DerivedMeasurement> &modelResult, const ModelReturn &input,
                          MeasurementModel,
                          MatrixBase<Derived> &H) {
            //The real component of the dual numbers are the result of the measurement model
            //The dual component vectors represent the rows of the jacobian matrix
                for (int i = 0; i < DOFOf<MatrixBase<DerivedMeasurement>>; ++i) {
                    H.row(i) = input[i].v;
                    modelResult[i] = input[i].a;
                }
            }




        /** update implementation for scalar measurements
        * @tparam ModelReturn The Type of the Addition between a Measurement and a dual component
        * @tparam MeasurementModel The Type of the Measurement Model Functor
        * @tparam Derived Type of the Covaraince Matrix
        * @param modelResult Result of the Measurement Model
        * @param input Result of the Measurement Model with a dual compnent added to the state
        * @param h The Measurement Model h(x)
        * @param H The resulting Jacobian. Calculated with dual numbers
        */
        template<typename MeasurementModel, typename Derived>
        void update_impl_(const ScalarType &, ScalarType &modelResult, const ceres::Jet<ScalarType,DOF> &input, MeasurementModel,
                          MatrixBase<Derived> &H) {
            //The real component of the dual numbers are the result of the measurement model
            //The dual component vectors represent the rows of the jacobian matrix
            H.row(0) = input.v;
            modelResult = input.a;
        }

        /**
         * Add an Offset to the Estimated State, if the State is a Manifold
         * @tparam Manifold The Type of Manifold used as State
         * @tparam Derived The Type of Matrix used as the Result of K*H
         * @param diff The Difference to be added to the state
         * @param KH The Multiplication of Kalman Gain and Measurement-Jacobian
         */
        template<typename Derived>
        void add_diff(const Manifold &, const MatrixType<DOF, 1> &diff, const MatrixBase<Derived> &KH) {
            //Add the Difference on the Estimated State
            State newMu = mu + diff;
            //Calculate the Jacobian of the Boxplus Function
            JacobianOf<State> D=transformReferenceJacobian(mu,newMu,diff);
            //Set the new Estimated Value
            mu = newMu;
            //Calculate the new Covariance Matrix
            sigma = D * (sigma - (KH * sigma)) * D.transpose();
        }


        /**
        * Add an Offset to the Estimated State, if the State is a CompoundManifold
         *
         * For CompoundManifolds we can optimise the Jacobian D, since each derivative is only dependent on one substate
        * @tparam Manifold The Type of Manifold used as State
        * @tparam Derived The Type of Matrix used as the Result of K*H
        * @param diff The Difference to be added to the state
        * @param KH The Multiplication of Kalman Gain and Measurement-Jacobian
        */
        template<typename Derived>
        void add_diff(const CompoundManifold &, const MatrixType<DOF, 1> &diff, const MatrixBase<Derived> &KH) {
            //check if compoundManifold is simple vector this may look a bit dirty but it allows to use the ADEKF_MANIFOLD for vector parts only without significant speed loss
            if constexpr (mu.MAN_DOF==0){
                add_diff<Derived>(diff,diff,KH);
                return;
            }

            //Add the Difference on the Estimated State
            State newMu = mu + diff;
            //Calculate  the Jacobian of the Boxplus Function
            JacobianOf<State> D = transformReferenceJacobian(mu,newMu,diff);
            //Set the new Estimated Value
            mu = newMu;
            //Calculate the new Covariance Matrix
            sigma = D * (sigma - (KH * sigma)) * D.transpose();
        }

        /**
      * Add an Offset to the Estimated State, if the State is a CompoundManifold
       *
       * For CompoundManifolds we can optimise the Jacobian D, since each derivative is only dependent on one substate
      * @tparam Manifold The Type of Manifold used as State
      * @tparam Derived The Type of Matrix used as the Result of K*H
      * @param diff The Difference to be added to the state
      * @param KH The Multiplication of Kalman Gain and Measurement-Jacobian
      */
        template<typename Derived, typename Nullspace>
        void add_diff(const CompoundManifold &, const MatrixType<DOF, 1> &diff, const MatrixBase<Derived> &KH, const MatrixBase<Nullspace> & N )   {
            //check if compoundManifold is simple vector this may look a bit dirty but it allows to use the ADEKF_MANIFOLD for vector parts only without significant speed loss
            if(mu.MAN_DOF==0){
                add_diff<Derived>(diff,diff,KH);
                return;
            }

            //Add the Difference on the Estimated State
            State newMu = mu + diff;
            //Calculate  the Jacobian of the Boxplus Function
            JacobianOf<State> D = transformReferenceJacobian(mu,newMu,diff);

            //Set the new Estimated Value
            mu = newMu;
            //nullspace constraint
            D=D-(D*N-N)*(N.transpose()*N).inverse()*N.transpose();
            //Calculate the new Covariance Matrix
            sigma = D * (sigma - (KH * sigma)) * D.transpose();
        }


        /**
         * Add an Offset to the Estimated State, if the Statesigma is a Matrix
         * @tparam Derived The Type of Matrix used as the Result of K*H
         * @param diff The Difference to be added to the state
         * @param KH The Multiplication of Kalman Gain and Measurement-Jacobian
         */
        template<typename Derived>
        void add_diff(const MatrixType<DOF, 1>&, const MatrixType<DOF, 1>& diff, const MatrixBase<Derived>& KH) {
            //Add the Difference on the State
            mu = mu + diff;
            //Calculate the new Covariance
            sigma = sigma - (KH * sigma);
        }
    };


    /**
     * General Deduction Template for the ADEKF based on StateRetriever.
     * This is needed so you can type ADEKF ekf(State,COV) without template arguments
     */
    template<typename DERIVED, typename COV_TYPE>
    ADEKF(const DERIVED &, const COV_TYPE &) -> ADEKF<typename StateInfo<DERIVED>::type>;

}
