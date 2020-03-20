#pragma once

#include<Eigen/Geometry>
#include "../ceres/rotation.h"

namespace adekf {
template<typename Scalar = float>
class SO3: public Eigen::Quaternion<Scalar> ,public Manifold{
public:
	using ScalarType = Scalar;
	static constexpr unsigned DOF = 3;

	SO3(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) :
			Eigen::Quaternion<Scalar>(w, x, y, z) {
		Eigen::Quaternion<Scalar>::normalize();
};

	SO3(const Eigen::Quaternion<Scalar> &src = Eigen::Quaternion<Scalar>::Identity()) :
			Eigen::Quaternion<Scalar>(src) {};

	SO3(const Eigen::Matrix<Scalar, DOF, DOF> &rotationMatrix) :
			Eigen::Quaternion<Scalar>(rotationMatrix) {};

	template<typename OtherScalar>
	auto operator*(const SO3<OtherScalar> &other) const {
		return adekf::SO3(this->w() * other.w() - this->x() * other.x() - this->y() * other.y() - this->z() * other.z(),
				this->w() * other.x() + this->x() * other.w() + this->y() * other.z() - this->z() * other.y(),
				this->w() * other.y() + this->y() * other.w() + this->z() * other.x() - this->x() * other.z(),
				this->w() * other.z() + this->z() * other.w() + this->x() * other.y() - this->y() * other.x());
	}

	template<typename Derived, typename OtherScalar = typename Eigen::internal::traits<Derived>::Scalar,
			typename ResultScalar = decltype(
					std::declval<Scalar>() * std::declval<OtherScalar>())>
	EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Eigen::Matrix<ResultScalar, DOF, 1> operator*(const Eigen::MatrixBase<Derived> &other) const {
        Eigen::Matrix<ResultScalar,3,1> uv = this->vec().cross(other);
        uv+=uv;
        return other + this->w() * uv + this->vec().cross(uv);
	}

	SO3<Scalar> inverse() const {
		return SO3<Scalar>(Eigen::Quaternion<Scalar>::inverse());
	}

	template<typename Derived, typename OtherScalar = typename Derived::Scalar>
	auto operator+(const Eigen::MatrixBase<Derived> &delta) const {
		OtherScalar expData[4];
		Eigen::Ref<const Eigen::Matrix<OtherScalar, DOF, 1>> buffer = delta;

		ceres::AngleAxisToQuaternion(buffer.data(), expData);

		SO3<OtherScalar> exp(expData[0], expData[1], expData[2], expData[3]);
		return *this *exp;
	}

	template<typename OtherScalar>
	auto operator-(const SO3<OtherScalar> &other) const {
		adekf::SO3 delta(other.inverse() * *this);
		using ResultScalar=typename decltype(delta)::ScalarType;
		Eigen::Matrix<ResultScalar, DOF, 1> result;
		ResultScalar deltaData[4] = { delta.coeffs().data()[3], delta.coeffs().data()[0], delta.coeffs().data()[1],
				delta.coeffs().data()[2] };
		ceres::QuaternionToAngleAxis(deltaData, result.data());
		return result;
	}

	friend std::ostream &operator<<(std::ostream &stream, const SO3<Scalar> &state) {
		return stream << state.coeffs().transpose() << " ";
	}

	friend std::istream &operator>>(std::istream &is, SO3<Scalar> &state) {
		Eigen::Matrix<Scalar, 4, 1> coeffs;
		is >> coeffs;
		state.coeffs() = coeffs.normalized();
		return is;
	}
};

using SO3f = SO3<float>;
using SO3d = SO3<double>;
}

