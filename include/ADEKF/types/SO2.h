#pragma once

#include<Eigen/Geometry>
#include "../ceres/rotation.h"
#include "../ceres/jet.h"
#include "type_utils.h"

using namespace Eigen;

namespace ADEKF {
    /**
     * A Representation of the SO(2) Group of Rotations in 2D Space
     * @tparam Scalar The used Scalar Type
     */
    template<typename Scalar = float>
    class SO2 : public Rotation2D<Scalar>, public Manifold {
    public:
        //Necessities of the Manifold Interface
        using ScalarType = Scalar;
        static constexpr unsigned DOF = 1;

        /**
         * Construct a rotation from an angle
         * @param angle The angle
         */
        SO2(const Scalar &angle = 0) : Rotation2D<Scalar>(anglemod(angle)) {};

        /**
         * Construct a rotation from an Eigen 2D Rotation
         * @param src The Eigen 2D Rotation
         */
        SO2(const Rotation2D<Scalar> &src) : Rotation2D<Scalar>(src) {};

        /**
         * Construct a rotation from a vector
         * @param vector The vector
         */
        SO2(const Matrix<Scalar, 2, 1> &vector) : Rotation2D<Scalar>(atan2(vector[1], vector[0])) {};

        template<int jetsize>
        static SO2<Scalar> extractFromJet(const SO2<ceres::Jet<Scalar,jetsize> >& other){
            return SO2{other.angle().a};
        }

        /**
         * Modulo for angles. Changes the result, but not the derivative
         * @tparam T The Scalar Type of the Jet
         * @tparam N The Size of the Jet
         * @param f The Jet to be normalized
         * @return A normalized Jet
         */
        template<typename T, int N> inline
        ceres::Jet<T, N> const &anglemod(const ceres::Jet<T, N> &f) {
            T angle = fmod(f.a + M_PI, M_PI*2);
            if(angle < T(0))
                angle += M_PI*2;
            return ceres::Jet<T, N>(angle - M_PI, f.v);
        }

        /**
         * Modulo for angles. Changes the result, but not the derivative
         * @tparam T The Scalar Type
         * @param f The Number to be normalized
         * @return The normalized number
         */
        template<typename T> inline
        T const     anglemod(const T f) {
            T angle = fmod(f + M_PI, M_PI*2);
            if(angle < T(0))
                angle += M_PI*2;
            return angle - M_PI;
        }

        template<typename OtherScalar, typename ResultScalar = decltype(std::declval<Scalar>() +
                                                                        std::declval<OtherScalar>())>
        SO2<ResultScalar> operator*(const SO2<OtherScalar> &other) const {
            return SO2<ResultScalar>(Rotation2D<Scalar>::angle() + other.angle());
        }

        SO2<Scalar> inverse() const {
            return SO2<Scalar>(Rotation2D<Scalar>::inverse());
        }

        template<typename Derived, typename OtherScalar = typename internal::traits<Derived>::Scalar, typename ResultScalar = decltype(
        std::declval<Scalar>() + std::declval<OtherScalar>())>
        SO2<ResultScalar>
        operator+(const MatrixBase<Derived> &delta) const {
            return SO2<ResultScalar>(delta[0]) * *this;
        }

        template<typename OtherScalar, typename ResultScalar = decltype(std::declval<Scalar>() -
                                                                        std::declval<OtherScalar>())>
        Matrix<ResultScalar, DOF, 1>
        operator-(const SO2<OtherScalar> &other) const {
            SO2<ResultScalar> delta = other.inverse() * *this;
            return Matrix<ResultScalar, DOF, 1>(anglemod(delta.angle()));
        }
        LOCAL_PARAMETRISATION(SO2,3);

        friend std::ostream &operator<<(std::ostream &stream, const SO2<Scalar> &state) {
            return stream << state.angle();
        }

        friend std::istream &operator>>(std::istream &is, SO2<Scalar> &state) {
            return is >> state.angle();
        }
    };
	
	using SO2f = SO2<float>;
    using SO2d = SO2<double>;
}
