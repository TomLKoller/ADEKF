//
// Created by tomlucas on 20.05.20.
//

#ifndef CYLINDEREXAMPLE_JETNORM_H
#define CYLINDEREXAMPLE_JETNORM_H
#include <Eigen/Core>
#include <ADEKF/ceres/jet.h>
#include "ADEKFUtils.h"

/**
 * Calculates the euclidean norm of 2 Jet values.
 *
 * This is required to calculate the limit of the derivative of the euclidean norm at norm=0.
 * The limit is calculated using the path where all values simultaneously approach 0.
 * @param vector the vector to compute the norm of
 * @return the norm of vector correctly differentiated if requested
 *
 */
template<typename OtherScalar, int N>
inline
static ceres::Jet<OtherScalar, N>
normDeduced(const Eigen::Matrix<ceres::Jet<OtherScalar, N>,3,1> &vector) {
    ceres::Jet<OtherScalar, N> out;

    OtherScalar const temp1 = sqrt(pow(vector(0).a, 2) + pow(vector(1).a, 2)+pow(vector(2).a, 2));
    OtherScalar const multiplier = 1. / (temp1);
    //Check limit condition the limit is 1/sqrt(2)
    OtherScalar const temp2 = temp1 == OtherScalar(0.) ? OtherScalar(1. / sqrt(3.)) : multiplier *  (vector(0).a);
    OtherScalar const temp3 = temp1 == OtherScalar(0.) ? OtherScalar(1. / sqrt(3.)) : multiplier * (vector(1).a);
    OtherScalar const temp4 = temp1 == OtherScalar(0.) ? OtherScalar(1. / sqrt(3.)) : multiplier * (vector(2).a);

    out.a = temp1;
    //Chain rule for derivative
    out.v = temp2 * vector(0).v + temp3 * vector(1).v+temp4*vector(2).v;
    return out;
}


//use type deduction so we can write norm(vector)

static double normDeduced(const Eigen::Matrix<double,3,1> & vector) {
    return vector.norm();
}
template<typename T>
static T normMatrix(const Eigen::Matrix<T,3,1> & vector){
    return normDeduced(vector);
}

template <typename Derived>
static typename adekf::StateInfo<Derived>::ScalarType norm (const Eigen::MatrixBase<Derived> & vector){
    return normMatrix<typename adekf::StateInfo<Derived>::ScalarType>(vector);
}





#endif //CYLINDEREXAMPLE_JETNORM_H
