#pragma once

#include<Eigen/Geometry>
#include "../ceres/rotation.h"
#include "../ADEKFUtils.h"
#include "SO3.h"
#include "type_utils.h"

namespace adekf {
    


template<typename Scalar>
class SO3RightInvariant: public Manifold, public Eigen::Quaternion<Scalar> {
public:
	using ScalarType = Scalar;
	static constexpr unsigned DOF = 3;
	
	SO3RightInvariant(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) :
			Eigen::Quaternion<Scalar>(w, x, y, z) {
		Eigen::Quaternion<Scalar>::normalize();
};

	SO3RightInvariant(const Eigen::Quaternion<Scalar> &src = Eigen::Quaternion<Scalar>::Identity()) :
			Eigen::Quaternion<Scalar>(src) {};


    SO3RightInvariant(const Eigen::Matrix<Scalar, DOF, DOF>  &rotationMatrix) :
            Eigen::Quaternion<Scalar>(rotationMatrix) {};
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	template<typename Derived>
	SO3RightInvariant(const Eigen::MatrixBase<Derived> & omega){
        Scalar expData[4];
        Eigen::Ref<const Eigen::Matrix<Scalar, DOF, 1>> buffer = omega;
        ceres::AngleAxisToQuaternion(buffer.data(), expData);
        *this=SO3RightInvariant<Scalar>(expData[0], expData[1], expData[2], expData[3]);
	}
	/**
	 * Constructor to assign from other scalars (if OtherScalar is convertible non explicitly to Scalar)
	 * @tparam OtherScalar  The scalar of the other SO3
	 * @param other the other SO3
	 */
    template<typename OtherScalar>
	SO3RightInvariant(const SO3RightInvariant<OtherScalar> & other):Eigen::Quaternion<Scalar>(other){

	}

    SO3RightInvariant(const Scalar * src): Eigen::Quaternion<Scalar>(src){

    }

    void toPointer(Scalar *dest){
        (Eigen::Map<Eigen::Matrix<Scalar,4,1> >(dest))=this->coeffs();
    }



    template<typename OtherScalar>
	auto operator*(const SO3RightInvariant<OtherScalar> &other) const {
		return adekf::SO3RightInvariant(this->w() * other.w() - this->x() * other.x() - this->y() * other.y() - this->z() * other.z(),
				this->w() * other.x() + this->x() * other.w() + this->y() * other.z() - this->z() * other.y(),
				this->w() * other.y() + this->y() * other.w() + this->z() * other.x() - this->x() * other.z(),
				this->w() * other.z() + this->z() * other.w() + this->x() * other.y() - this->y() * other.x());
	}

	template<typename Derived, typename OtherScalar = typename Eigen::internal::traits<Derived>::Scalar,
			typename ResultScalar = decltype(
					std::declval<Scalar>() * std::declval<OtherScalar>())>
	EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE QuatRotation<adekf::SO3RightInvariant<Scalar>,Eigen::Matrix<OtherScalar,3,1>> operator*(const Eigen::MatrixBase<Derived> &other) const {
        return QuatRotation<adekf::SO3RightInvariant<Scalar>,Eigen::Matrix<OtherScalar,3,1>>(*this,other);
	}

	SO3RightInvariant<Scalar> inverse() const {
		return SO3RightInvariant<Scalar>(Eigen::Quaternion<Scalar>::inverse());
	}
	SO3RightInvariant<Scalar> conjugate() const{
	    return SO3RightInvariant<Scalar>(Eigen::Quaternion<Scalar>::conjugate());
	}


	template<typename Derived, typename OtherScalar = typename Derived::Scalar, typename ResultScalar=ADEKF_PLUSRESULT(Scalar, OtherScalar)>
	SO3RightInvariant<ResultScalar> operator+(const Eigen::MatrixBase<Derived> &delta) const {
        return SO3RightInvariant<OtherScalar>(delta)* *this;
	}

    /**
     * Operator for ceres local parameterization
     */
    template<typename T>
    bool operator()(const T* x,const T* delta,T* x_plus_delta) const {
        (SO3RightInvariant<T>(x)+Eigen::Map<const Eigen::Matrix<T,DOF,1>>(delta)).toPointer(x_plus_delta);
        return true;
    }


	template<typename OtherScalar>
	auto operator-(const SO3RightInvariant<OtherScalar> &other) const {
		adekf::SO3RightInvariant delta(*this *other.conjugate());
		using ResultScalar=typename decltype(delta)::ScalarType;
		Eigen::Matrix<ResultScalar, DOF, 1> result;
		ResultScalar deltaData[4] = { delta.coeffs().data()[3], delta.coeffs().data()[0], delta.coeffs().data()[1],
				delta.coeffs().data()[2] };//Eigen storage order is x y z w while ceres is w x y z
		ceres::QuaternionToAngleAxis(deltaData, result.data());
		return result;
	}

	LOCAL_PARAMETRISATION(SO3RightInvariant,4);




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

template<typename Derived>
SO3RightInvariant(const Eigen::MatrixBase<Derived> & ) -> SO3<typename Derived::Scalar>;

using SO3RightInvariantf = SO3RightInvariant<float>;
using SO3RightInvariantd = SO3RightInvariant<double>;
}








