//
// Created by tomlucas on 06.04.20.
//

#include <gtest/gtest.h>
#include "DirectionVector.h"
#include <ADEKF/ceres/jet.h>

TEST (ManifoldTests , DirectionVectorExpLogTest) {
    Eigen::Matrix<ceres::Jet<double,2>,2,1> test_vector=test_vector.setZero();
    test_vector(0)=ceres::Jet<double,2>(1,0);
    test_vector(1)=ceres::Jet<double,2>(0,1);
    auto exp=adekf::DVd::getExp(test_vector).eval();
    auto log=adekf::DVd::getLog(exp);
    for(int i=0;i < 2; i++) {
        //Test real part
        ASSERT_EQ(test_vector(i), log(i));
        //test dual part
        ASSERT_EQ(test_vector(i).v, log(i).v);
    }
}


TEST (ManifoldTests , DirectionVectorExpLogLimitTest ) {
    Eigen::Matrix<ceres::Jet<double,2>,2,1> test_vector=test_vector.setZero();
    test_vector(0)=ceres::Jet<double,2>(0,0);
    test_vector(1)=ceres::Jet<double,2>(0,1);
    auto exp=adekf::DVd::getExp(test_vector).eval();
    auto log=adekf::DVd::getLog(exp);
    for(int i=0;i < 2; i++) {
        //Test real part
        ASSERT_EQ(test_vector(i), log(i));
        //test dual part
        ASSERT_EQ(test_vector(i).v, log(i).v);
    }
}