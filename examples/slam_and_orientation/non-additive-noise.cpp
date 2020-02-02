//
// Created by tomlucas on 30.01.20.
//

#include "ADEKF.h"
#include <Eigen/Core>
#include "types/SO3.h"
//Position estimate in m
void test(){
    // Set the start orientation to identity
    adekf::ADEKF ekf(adekf::SO3d(),Eigen::Matrix3d::Identity());
    //Create dynamic model which adds the angular velocity times time_diff to the position
    auto dynamic_model=[](auto & state, Eigen::Vector3d angular_rate, double time_diff){
        state= state + angular_rate * time_diff;
    };
    //create a measurement_model which gets a orientation measurement in mm
    auto measurement_model=[](auto state){
        return state;
    };
    //time step length
    double time_diff=0.01;
    // input measurement -> random velocity
    Eigen::Vector3d angular_rate=angular_rate.Random();
    //A orientation measurement -> random orientation;
    adekf::SO3d orientation(Eigen::Quaterniond::UnitRandom());
    //dynamic noise
    Eigen::Matrix3d process_noise=process_noise.Identity()*0.01;
    // measurement noise
    Eigen::Matrix3d measurement_noise=measurement_noise.Identity();
    //prediction:pass model, noise followed by the parameters you want in your dynamic _model (variadic acceptance)
    ekf.predict(dynamic_model, process_noise, angular_rate, time_diff);
    //measurement update: pass model , noise and a measurement followed by parameters you need for your model  (variadic acceptance)
    ekf.update(measurement_model, measurement_noise, orientation);
}