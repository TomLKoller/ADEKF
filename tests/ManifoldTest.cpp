//
// Created by tomlucas on 06.04.20.
//

#include <gtest/gtest.h>
#include "types/DirectionVector.h"
#include "ceres/jet.h"

/**
 * Tests whether log(exp([1 0 ]) yields the correct result for derivatives
 */
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



/**
 * Tests whether log(exp([0 0 ])  (limit case) yields the correct result for derivatives
 */
TEST (ManifoldTests , DirectionVectorExpLogLimitTest ) {
    Eigen::Matrix<ceres::Jet<double,2>,2,1> test_vector=test_vector.setZero();
    test_vector(0)=ceres::Jet<double,2>(0,0);
    test_vector(0).v=Eigen::Vector2d(0.2,0.3);
    test_vector(1)=ceres::Jet<double,2>(0,1);
    test_vector(1).v=Eigen::Vector2d(0.4,0.5);
    auto exp=adekf::DVd::getExp(test_vector).eval();
    auto log=adekf::DVd::getLog(exp);
    for(int i=0;i < 2; i++) {
        //Test real part
        ASSERT_EQ(test_vector(i), log(i));
        //test dual part
        ASSERT_EQ(test_vector(i).v, log(i).v);
    }
}

TEST (ManifoldTests , DirectionVectorBoxPlusMinusTest ) {
    Eigen::Matrix<ceres::Jet<double,3>,3,1> test_vector=test_vector.setZero();
    test_vector(0)=ceres::Jet<double,3>(0,0);
    test_vector(0).v=Eigen::Vector3d(0.2,0.3,0.7);
    test_vector(1)=ceres::Jet<double,3>(0,1);
    test_vector(1).v=Eigen::Vector3d(0.4,0.5,0.6);
    test_vector(2)=ceres::Jet<double,3>(0,2);
    test_vector(2).v=Eigen::Vector3d(0.4,0.5,0.6);


    Eigen::Matrix<ceres::Jet<double,3>,3,1> test_vector2=test_vector2.setZero();
    test_vector2(0)=ceres::Jet<double,3>(1,0);
    test_vector2(0).v=Eigen::Vector3d(0.2,0.3,0.7);
    test_vector2(1)=ceres::Jet<double,3>(1,1);
    test_vector2(1).v=Eigen::Vector3d(0.4,0.5,0.6);
    test_vector2(2)=ceres::Jet<double,3>(1,2);
    test_vector2(2).v=Eigen::Vector3d(0.4,0.5,0.6);




    auto minus=test_vector2-test_vector;
    auto result=test_vector+minus;
    for(int i=0;i < 3; i++) {
        //Test real part
        ASSERT_EQ(test_vector2(i), result(i));
        //test dual part
        ASSERT_EQ(test_vector2(i).v, result(i).v);
    }
}