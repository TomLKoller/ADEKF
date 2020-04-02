#pragma once

#include<Eigen/Geometry>
#include <ADEKF/ADEKFUtils.h>


namespace adekf {


//Todo test extensively
    template<typename Scalar>
    class DirectionVector : public Manifold, public Eigen::Matrix<Scalar, 3, 1> {
    public:
        using ScalarType = Scalar;
        static constexpr unsigned DOF = 2;

        DirectionVector(const Scalar &x, const Scalar &y, const Scalar &z) :
                Eigen::Matrix<Scalar, 3, 1>(x, y, z) {
            this->normalize();
        };

        template<class Derived>
        DirectionVector(const Eigen::MatrixBase<Derived> &src) :
                Eigen::Matrix<Scalar, 3, 1>(src) {
            this->normalize();
        }

        DirectionVector()  {
            this->setZero();
        }

        static double norm(const double a, const double b) {
            return sqrt(pow(a, 2) + pow(b, 2));
        }


        template<typename OtherScalar, int N>
        inline
        static ceres::Jet<OtherScalar, N>
        norm(const ceres::Jet<OtherScalar, N> &a, const ceres::Jet<OtherScalar, N> &b) {
            ceres::Jet<OtherScalar, N> out;

            OtherScalar const temp1 = sqrt(pow(a.a, 2) + pow(b.a, 2));
            OtherScalar const multiplier = 1. / (temp1);
            OtherScalar const temp2 = temp1 == OtherScalar(0.) ? OtherScalar(1. / sqrt(2.)) : multiplier * (a.a);
            OtherScalar const temp3 = temp1 == OtherScalar(0.) ? OtherScalar(1. / sqrt(2.)) : multiplier * (b.a);

            out.a = temp1;
            out.v = temp2 * a.v + temp3 * b.v;
            return out;
        }

        /**
         *  Matrix that turns [1 0 0 ] to vector
         *
         */
        template<typename OtherScalar>
        static Eigen::Matrix<OtherScalar, 3, 3> getRx(const Eigen::Matrix<OtherScalar, 3, 1> &vector) {
            OtherScalar r = norm(vector(1), vector(2));

            OtherScalar a = asin(vector(2));
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
         * Exponential form of delta
         * @param delta an 2d offset of a 3d vector
         * @return
         */
        template<typename OtherScalar>
        static Eigen::Matrix<OtherScalar, 3, 1> getExp(const Eigen::Matrix<OtherScalar, 2, 1> &delta) {
            OtherScalar n = norm(delta(0), delta(1));
            Eigen::Matrix<OtherScalar, 3, 1> exp;
            exp(0) = cos(n);
            OtherScalar sinc = n == OtherScalar(0.) ? cos(n) : sin(n) / n;
            exp(1) = sinc * delta(0);
            exp(2) = sinc * delta(1);
            assert_finite(exp);
            return exp;
        }

        template<typename Derived, typename OtherScalar=typename Derived::Scalar>
        static Eigen::Matrix<OtherScalar, 2, 1> getLog(const Eigen::MatrixBase<Derived> &vector) {
            Eigen::Matrix<OtherScalar, 2, 1> lg;
            OtherScalar n = norm(vector(1), vector(2));
            if (n == OtherScalar(0)) {

                lg(0, 0) = atan2(n, vector(0));
                lg(1, 0) = n;
            } else {
                lg = atan2(n, vector(0)) / n * vector.template block<2, 1>(1, 0);
            }
            assert_finite(lg);
            return lg;
        }


        template<typename Derived, typename OtherScalar = typename Derived::Scalar>
        auto operator+(const Eigen::MatrixBase<Derived> &delta) const {
            DirectionVector<ADEKF_PLUSRESULT(Scalar, OtherScalar)> temp(
                    getRx<Scalar>(*this) * getExp<OtherScalar>(delta));
            assert_finite(temp);
            return temp;
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

    using SO3f = SO3<float>;
    using SO3d = SO3<double>;
}





