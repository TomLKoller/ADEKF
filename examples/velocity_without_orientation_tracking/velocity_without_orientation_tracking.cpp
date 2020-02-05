//
//

#include <ADEKF/ADEKF.h>
int main(){
    //Initialize adekf with 6 states
    adekf::ADEKF ekf(Eigen::Vector3d::Zero(),Eigen::Matrix3d::Identity());
    //dynamic model which propagates the acceleration on to the velocity
    auto dynamic_model=[](auto &state, const Eigen::Vector3d & accel, double time_diff){
        state+=accel*time_diff;
    };

    //external velocity measurement
    auto measurement_model=[](auto state){
        return state;
    };

    //Initialise measurements
    Eigen::Vector3d __IN_acceleration=__IN_acceleration.Random();
    Eigen::Vector3d __IN_velocity=__IN_velocity.Random();
    double __IN_time_diff=0.01;
    // Initialise noise
    Eigen::Matrix3d __IN_processNoise=__IN_processNoise.Identity();
    Eigen::Matrix3d __IN_measurement_noise=__IN_measurement_noise.Identity();
    //Call predict to apply acceleration on velocity
    ekf.predict(dynamic_model,__IN_processNoise,__IN_acceleration,__IN_time_diff);

    //Correct estimation with a measurement update
    ekf.update(measurement_model,__IN_measurement_noise,__IN_velocity);


    //print out ekf state
    std::cout << ekf.mu << std::endl;


    return 0;
}
