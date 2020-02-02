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
	Eigen::Matrix<ResultScalar, DOF, 1> operator*(const Eigen::MatrixBase<Derived> &other) const {
		const Scalar scale = Scalar(1)
				/ sqrt(this->w() * this->w() + this->x() * this->x() + this->y() * this->y() + this->z() * this->z());

		// Make unit-norm version of q.
		const Scalar unit[4] = { scale * this->w(), scale * this->x(), scale * this->y(), scale * this->z(), };

		const Scalar t2 = unit[0] * unit[1];
		const Scalar t3 = unit[0] * unit[2];
		const Scalar t4 = unit[0] * unit[3];
		const Scalar t5 = -unit[1] * unit[1];
		const Scalar t6 = unit[1] * unit[2];
		const Scalar t7 = unit[1] * unit[3];
		const Scalar t8 = -unit[2] * unit[2];
		const Scalar t9 = unit[2] * unit[3];
		const Scalar t1 = -unit[3] * unit[3];

		return Eigen::Matrix<ResultScalar, 3, 1>(
				Scalar(2) * ((t8 + t1) * other.x() + (t6 - t4) * other.y() + (t3 + t7) * other.z()) + other.x(),
				Scalar(2) * ((t4 + t6) * other.x() + (t5 + t1) * other.y() + (t9 - t2) * other.z()) + other.y(),
				Scalar(2) * ((t7 - t3) * other.x() + (t2 + t9) * other.y() + (t5 + t8) * other.z()) + other.z());
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

