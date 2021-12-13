#pragma once

#include<Eigen/Geometry>
#include "../ceres/rotation.h"
#include "../ADEKFUtils.h"
#include "type_utils.h"


namespace adekf {
    /**
     * Expression Class for Quaternion Rotation on non primitive scalar types (as ceres::Jet)
     */
    template <class QuaternionType, class ArgType,typename ResultScalar = decltype(
    std::declval<typename QuaternionType::ScalarType >() * std::declval<typename Eigen::internal::traits<ArgType>::Scalar>())>
    class QuatRotation;

    template <class QuaternionType, class ArgType,typename ResultScalar >
    class QuatRotation : public Eigen::MatrixBase<QuatRotation<QuaternionType,ArgType,ResultScalar> >
    {

    public:
        typedef ResultScalar CoeffReturnType;
        typedef Eigen::Matrix<ResultScalar,3,1> RESULT_TYPE;
        QuatRotation(const QuaternionType & quat,const ArgType& arg)
        {
            RESULT_TYPE uv = quat.vec().cross(arg);
            uv+=uv;
            result= arg + quat.w() * uv + quat.vec().cross(uv);
        }
        typedef typename Eigen::internal::ref_selector<QuatRotation>::type Nested;
        typedef Eigen::Index Index;
        Index rows() const { return 3; }
        Index cols() const { return 1; }
        RESULT_TYPE result;
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };



template<typename Scalar>
class SO3: public Manifold, public Eigen::Quaternion<Scalar> {
public:
	using ScalarType = Scalar;
	static constexpr unsigned DOF = 3;

	SO3(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) :
			Eigen::Quaternion<Scalar>(w, x, y, z) {
		Eigen::Quaternion<Scalar>::normalize();
};

	SO3(const Eigen::Quaternion<Scalar> &src = Eigen::Quaternion<Scalar>::Identity()) :
			Eigen::Quaternion<Scalar>(src) {};


    SO3(const Eigen::Matrix<Scalar, DOF, DOF>  &rotationMatrix) :
            Eigen::Quaternion<Scalar>(rotationMatrix) {};
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	template<typename Derived>
	SO3(const Eigen::MatrixBase<Derived> & omega){
        Scalar expData[4];
        Eigen::Ref<const Eigen::Matrix<Scalar, DOF, 1>> buffer = omega;
        ceres::AngleAxisToQuaternion(buffer.data(), expData);
        *this=SO3<Scalar>(expData[0], expData[1], expData[2], expData[3]);
	}
	/**
	 * Constructor to assign from other scalars (if OtherScalar is convertible non explicitly to Scalar)
	 * @tparam OtherScalar  The scalar of the other SO3
	 * @param other the other SO3
	 */
    template<typename OtherScalar>
	SO3(const SO3<OtherScalar> & other):Eigen::Quaternion<Scalar>(other){

	}

    SO3(const Scalar * src): Eigen::Quaternion<Scalar>(src){

    }

    void toPointer(Scalar *dest){
        (Eigen::Map<Eigen::Matrix<Scalar,4,1> >(dest))=this->coeffs();
    }



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
	EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE QuatRotation<adekf::SO3<Scalar>,Eigen::Matrix<OtherScalar,3,1>> operator*(const Eigen::MatrixBase<Derived> &other) const {
        return QuatRotation<adekf::SO3<Scalar>,Eigen::Matrix<OtherScalar,3,1>>(*this,other);
	}

	SO3<Scalar> inverse() const {
		return SO3<Scalar>(Eigen::Quaternion<Scalar>::inverse());
	}
	SO3<Scalar> conjugate() const{
	    return SO3<Scalar>(Eigen::Quaternion<Scalar>::conjugate());
	}


	template<typename Derived, typename OtherScalar = typename Derived::Scalar, typename ResultScalar=ADEKF_PLUSRESULT(Scalar, OtherScalar)>
	SO3<ResultScalar> operator+(const Eigen::MatrixBase<Derived> &delta) const {
        return *this * SO3<OtherScalar>(delta);
	}

   


	template<typename OtherScalar>
	auto operator-(const SO3<OtherScalar> &other) const {
		adekf::SO3 delta(other.conjugate() * *this);
		using ResultScalar=typename decltype(delta)::ScalarType;
		Eigen::Matrix<ResultScalar, DOF, 1> result;
		ResultScalar deltaData[4] = { delta.coeffs().data()[3], delta.coeffs().data()[0], delta.coeffs().data()[1],
				delta.coeffs().data()[2] };//Eigen storage order is x y z w while ceres is w x y z
		ceres::QuaternionToAngleAxis(deltaData, result.data());
		return result;
	}

    LOCAL_PARAMETRISATION(SO3,4);




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
SO3(const Eigen::MatrixBase<Derived> & ) -> SO3<typename Derived::Scalar>;

using SO3f = SO3<float>;
using SO3d = SO3<double>;
}






namespace Eigen {
    namespace internal {
        /*
         * Traits class for the quaternion rotation expression
         */
        template <class QuaternionType,class ArgType>
        struct traits<adekf::QuatRotation<QuaternionType,ArgType>>
        {
            typedef Eigen::Dense StorageKind;
            typedef Eigen::MatrixXpr XprKind;
            typedef typename ArgType::StorageIndex StorageIndex;
            typedef decltype(std::declval<typename QuaternionType::ScalarType >() * std::declval<typename Eigen::internal::traits<ArgType>::Scalar>()) Scalar;
            enum {
                Flags = Eigen::ColMajor,
                RowsAtCompileTime = 3,
                ColsAtCompileTime = 1,
                MaxRowsAtCompileTime = 3,
                MaxColsAtCompileTime = 1
            };
        };

        /**
         * Evaluator for quaternion rotation. It simply reads the calculated vector
         */
        template <class QuaternionType, class ArgType>
        struct evaluator<adekf::QuatRotation<QuaternionType,ArgType> >
        : evaluator_base<adekf::QuatRotation<QuaternionType,ArgType> >
    {
        typedef adekf::QuatRotation<QuaternionType,ArgType> XprType;
        typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
        typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
        typedef typename XprType::CoeffReturnType CoeffReturnType;
        enum {
            CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
            Flags = Eigen::ColMajor | Eigen::LinearAccessBit
        };

        evaluator(const XprType& xpr)
                : m_argImpl(xpr.result), m_rows(xpr.rows())
        { }
        CoeffReturnType coeff(Index row, Index col=0) const
        {
           return m_argImpl.coeff(row,col);
        }
        evaluator<Eigen::Matrix<CoeffReturnType,3,1>> m_argImpl;
        const Index m_rows;
    };
}
}


