#pragma once

#include<Eigen/Geometry>
#include <ADEKF/ADEKFUtils.h>
#include <ADEKF/ceres/jet.h>


namespace adekf {


//Todo test extensively
/**
 * This class represents a direction vector manifold with unit norm.
 *
 * With this, you can keep a vector in the state which has always unit norm and is only rotated when changing.
 * This is intended to track observable directions of the orientation e.g. the gravity in most cases.
 * It behaves in calculations like a normal Eigen::Vector3, except for the operator+ with size 2 matrices (boxplus rotation)
 * and the operator- which returns a 2 DOF difference to the other vector.
 *
 *
 */
    template<typename Scalar>
    class DirectionVector : public Manifold, public Eigen::Matrix<Scalar, 3, 1> {
    public:
        using ScalarType = Scalar;
        static constexpr unsigned DOF = 2;
        /**
         * Construct from 3 scalars.
         *
         * The vector is normalized to ensure unit norm.
         * @param x  first entry
         * @param y  second entry
         * @param z  third entry
         */
        DirectionVector(const Scalar &x, const Scalar &y, const Scalar &z) :
                Eigen::Matrix<Scalar, 3, 1>(x, y, z) {
            this->normalize();
        };

        /**
         * Constructs the DirectionVector from another vector.
         *
         * THe vector is normalized to ensure unit norm.
         * @tparam Derived The type of the other vector. Must be a 3x1 vector of any scalar type.
         * @param src The vector to copy
         */
        template<class Derived>
        DirectionVector(const Eigen::MatrixBase<Derived> &src) :
                Eigen::Matrix<Scalar, 3, 1>(src) {
            this->normalize();
        }

        /**
         * Default constructor which initializes the vector with 0,0,1.
         *
         *  Default is required to be used in the ADEKF.
         */
        DirectionVector():Eigen::Matrix<Scalar,3,1>(0.,0.,1.)  {}


        /**
         * Calculates the euclidean norm of 2 double values.
         * @param a first value
         * @param b second value
         * @return sqrt(a^2+b^2)
         */
        static double norm(const double a, const double b) {
            return sqrt(pow(a, 2) + pow(b, 2));
        }


        /**
         * Calculates the euclidean norm of 2 Jet values.
         *
         * This is required to calculate the limit of the derivative of the euclidean norm at norm=0.
         * The limit is calculated using the path where both values simultaneously approach 0.
         *
         * @tparam OtherScalar the scalar type of the Jet
         * @tparam N the DOF of the Jet
         * @param a first value
         * @param b second value
         * @return r.a=sqrt(a.a^2+b.a^2), r.v = limit(a->0,b->0)  (a.a*a.v+b.a*b.v)/sqrt(a.a^2+b.a^2)
         */
        template<typename OtherScalar, int N>
        inline
        static ceres::Jet<OtherScalar, N>
        norm(const ceres::Jet<OtherScalar, N> &a, const ceres::Jet<OtherScalar, N> &b) {
            ceres::Jet<OtherScalar, N> out;

            OtherScalar const temp1 = sqrt(pow(a.a, 2) + pow(b.a, 2));
            OtherScalar const multiplier = 1. / (temp1);
            //Check limit condition the limit is 1/sqrt(2)
            OtherScalar const temp2 = temp1 == OtherScalar(0.) ? OtherScalar(1. / sqrt(2.)) : multiplier * (a.a);
            OtherScalar const temp3 = temp1 == OtherScalar(0.) ? OtherScalar(1. / sqrt(2.)) : multiplier * (b.a);

            out.a = temp1;
            //Chain rule for derivative
            out.v = temp2 * a.v + temp3 * b.v;
            return out;
        }

        /**
         *   Calculates the Matrix that turns [1 0 0] to vector.
         *
         *  Basically, this is the rotation encoded by the direction vector. If the vector is the gravity vector, this matrix is the rotation without yaw.
         * @tparam OtherScalar The scalar type of the vector
         * @param vector  the direction vector
         * @return The 3x3 rotation matrix that turns [1 0 0 ] to vector.
         */
        template<typename OtherScalar>
        static Eigen::Matrix<OtherScalar, 3, 3> getRx(const Eigen::Matrix<OtherScalar, 3, 1> &vector) {
            OtherScalar r = norm(vector(1), vector(2));

            //limit case
            OtherScalar a = asin(vector(2));
            //usual case
            if (vector(1) != OtherScalar(0) or vector(2) != OtherScalar(0)) {
                a = atan2(vector(2), vector(1));
            }
            Eigen::Matrix<OtherScalar, 3, 3> Rx;
            Rx.col(0) = vector;
            Rx(0, 1) = -r;
            Rx(1, 1) = vector(0) * cos(a);
            Rx(2, 1) = vector(0) * sin(a);
            Rx(0, 2) = OtherScalar(0);
            Rx(1, 2) = -sin(a);
            Rx(2, 2) = cos(a);
            assert_finite(Rx);
            return Rx;
        }

        /**
         * Exponential form of delta.
         * Analog to the euler rodrigues exponential function for SO3 rotations in
         * Integrating Generic Sensor Fusion Algorithms with Sound State Representations through Encapsulation of Manifolds. Hertzberg 2013
         *
         * Construts a 3x1 vector out of the 2 angles encoded in delta.
         * e.g. [0 0] -> [1 0 0]
         *
         * @param delta an 2d angle offset  of a 3d vector
         * @return The exponential form of delta.
         */
        template<typename OtherScalar>
        static Eigen::Matrix<OtherScalar, 3, 1> getExp(const Eigen::Matrix<OtherScalar, 2, 1> &delta) {
            OtherScalar n = norm(delta(0), delta(1));
            Eigen::Matrix<OtherScalar, 3, 1> exp;
            exp(0) = cos(n);
            //Norm check to apply the limit on jets
            OtherScalar sinc = n == OtherScalar(0.) ? cos(n) : sin(n) / n;
            exp(1) = sinc * delta(0);
            exp(2) = sinc * delta(1);
            return exp;
        }

        /**
         * The vector logarithm which derives the 2d angles of a vector.
         *
         * This is the inverse operation of getExp. Analog to the Matrix logarithm in
         * Integrating Generic Sensor Fusion Algorithms with Sound State Representations through Encapsulation of Manifolds. Hertzberg 2013.
         *
         * @tparam Derived The type of the passed vector. Must be a 3x1 Vector
         * @tparam OtherScalar The scalar type of the vector
         * @param vector The vector to retrieve the 2d angles from
         * @return The 2d angles of vector
         */
        template<typename Derived, typename OtherScalar=typename Derived::Scalar>
        static Eigen::Matrix<OtherScalar, 2, 1> getLog(const Eigen::MatrixBase<Derived> &vector) {
            Eigen::Matrix<OtherScalar, 2, 1> lg;
            OtherScalar delta_n = norm(vector(1), vector(2));
            //limit case where vector is [1 0 0]
            if (delta_n == OtherScalar(0.)) {
                lg =vector.template segment<2>(1) /vector(0);
            } else {
                lg = atan2(delta_n, vector(0)) / delta_n * vector.template segment<2>(1);
            }
            return lg;
        }


        template<typename Derived, typename OtherScalar = typename Derived::Scalar>
        auto operator+(const Eigen::MatrixBase<Derived> &delta) const {
            return  DirectionVector<ADEKF_PLUSRESULT(Scalar, OtherScalar)>(
                    getRx<Scalar>(*this) * getExp<OtherScalar>(delta));;
        }

        template<typename OtherScalar>
        Eigen::Matrix<ADEKF_MINUSRESULT(Scalar, OtherScalar), DOF, 1>
        operator-(const DirectionVector<OtherScalar> &other) const {
            return getLog(getRx<OtherScalar>(other).transpose() * *this);
        }

        friend std::ostream &operator<<(std::ostream &stream, const DirectionVector<Scalar> &state) {
            return stream << state.transpose() << " ";
        }

        friend std::istream &operator>>(std::istream &is, DirectionVector<Scalar> &state) {
            is >> state;
            return is;
        }
    };

    template<typename Derived>
    DirectionVector(const Eigen::MatrixBase<Derived> &) ->DirectionVector<typename Derived::Scalar>;

    using DVf = DirectionVector<float>;
    using DVd = DirectionVector<double>;
}





