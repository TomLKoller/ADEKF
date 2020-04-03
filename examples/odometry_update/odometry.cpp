#include <ADEKF/ADEKF.h>
#include <ADEKF/ManifoldCreator.h>
#include <ADEKF/types/SO3.h>
#include <thread>
#include "adekf_viz.h"
#include "DirectionVector.h"
#include "GaussianNoiseVector.h"

ADEKF_MANIFOLD(Pose3d, ((adekf::SO3, orientation)), (3, position))

ADEKF_MANIFOLD(Position, ((adekf::SO3, orientation)), (3, position), (3, velocity))

ADEKF_MANIFOLD(Old_Position,((adekf::SO3, orientation)) , (3, position), (3, velocity), (3, old_position))
ADEKF_MANIFOLD(GravityState,((adekf::SO3, orientation))((adekf::DirectionVector,gravity)) , (3, position), (3, velocity))

#define GRAVITY_CONSTANT 9.81

int main(int argc, char * argv[]) {


    //Simple Example with no orientation
    adekf::viz::initGuis(argc,argv);
    Eigen::Vector3d acc{1, 0, 0};
    Eigen::Vector3d omega{0,0.1,0.};
    Eigen::Vector3d base_vel{1, 0, 0};
    Eigen::Vector3d vel{base_vel};
    Eigen::Vector3d pos{20,0,0};

    Eigen::Vector3d gravity{0,0,GRAVITY_CONSTANT};

    adekf::SO3d orient;
    adekf::StateInfo<adekf::SO3d> test;
    adekf::ADEKF ekf_body(Position(orient, pos, vel), Eigen::Matrix<double, 9, 9>::Identity());
    adekf::ADEKF ekf_true(Position(orient,pos, vel), Eigen::Matrix<double, 9, 9>::Identity());
    adekf::ADEKF ekf_world(Position(orient, pos, vel),
                           Eigen::Matrix<double, 9, 9>::Identity());
    adekf::ADEKF ekf_gravity(GravityState(orient,Eigen::Vector3d(0,0,1),pos,vel),Eigen::Matrix<double,11,11>::Identity());
    adekf::viz::displayPose(&ekf_body, "Salmon");
    adekf::viz::displayPose(&ekf_true, "Green");
    adekf::viz::displayPose(&ekf_world, "Blue");
    adekf::viz::displayPose(&ekf_gravity,"Yellow");

    //adekf::viz::HeatMap::displayCovariance(&ekf_true, "Ground truth Filter");
    adekf::viz::displayCovariance(&ekf_body, "Body velocity Filter");
    adekf::viz::displayCovariance(&ekf_world, "World velocity Filter");
    adekf::viz::displayCovariance(&ekf_gravity,"Gravity vector Filter");
    double deltaT = 0.02;
    auto dynamic_model_body = [&gravity](auto &state, auto noise, auto acc, auto omega, double deltaT) {
        state.orientation=state.orientation*adekf::SO3((omega+NOISE(0,3))*deltaT);
        auto world_acc=state.orientation*(acc+NOISE(3,3))-gravity;
        state.position += state.orientation*state.velocity * deltaT +0.5* world_acc * deltaT*deltaT;
        state.velocity +=state.orientation.conjugate()*world_acc * deltaT+omega.cross(state.velocity)*deltaT;
    };
    auto dynamic_model_world = [&gravity](auto &state, auto noise, auto acc, auto omega, double deltaT) {
        state.orientation=state.orientation*adekf::SO3((omega+NOISE(0,3))*deltaT);
        auto world_acc=state.orientation*(acc+NOISE(3,3))-gravity;
        state.position += state.velocity * deltaT +0.5* world_acc * deltaT*deltaT;
        state.velocity +=world_acc * deltaT;
    };

    auto dynamic_model_gravity = [](auto &state, auto noise, auto acc, auto omega, double deltaT) {
        auto noiseRotation=adekf::SO3((omega+NOISE(0,3))*deltaT);
        state.orientation=state.orientation*noiseRotation;
        state.gravity=noiseRotation.conjugate()*state.gravity;
        auto body_acc= acc + NOISE(3, 3) - state.gravity * GRAVITY_CONSTANT;
        state.position += state.orientation*(state.velocity * deltaT + 0.5 * body_acc * deltaT * deltaT);
        state.velocity +=  (body_acc  + omega.cross(state.velocity) )* deltaT;
    };
    auto measurement_model_body = [](auto state) { return state.velocity; };
    auto measurement_model_world = [](auto state) { return state.orientation.conjugate() * (state.velocity); };

    auto consistency_model =[&gravity](auto state){return adekf::DirectionVector(state.orientation.conjugate()*gravity);};
    Eigen::Matrix<double, 6, 6> cov=cov.Identity() *0.1;
    cov(3,3)=cov(4,4)=cov(5,5)=0.1;
    adekf::GaussianNoiseVector<double,3> random_vel{0,0.1}, random_acc{0,0.1},random_omega{0,0.001};
    Eigen::Vector2d consistency_target=Eigen::Vector2d::Zero();
    //adekf::viz::LinePlot::plotVector(ekf_true.mu.position, "acc",1000,"xyz");
    //LOG_STREAM << ekf_gravity.sigma  <<std::endl LOG_END
    Eigen::Matrix<double, 11,11> erase_sigma=erase_sigma.Zero();
    erase_sigma(3,3)=erase_sigma(4,4)=erase_sigma(8,8)=erase_sigma(9,9)=erase_sigma(10,10)=1;
    Eigen::Matrix<double, 11,11> start_cov=(start_cov.Identity()-erase_sigma)*0.1;
    Eigen::Matrix<double,11,1> nullspace=nullspace.Zero();
    nullspace(2)=1;
    std::thread loop([&]() {
                        int cnt=0;
                         while (!adekf::viz::PoseRenderer::isDone()) {
                             random_acc.poll();
                             random_omega.poll();
                             //random_acc=random_omega=random_omega.Zero();
                             acc=-omega.cross(base_vel)+(ekf_true.mu.orientation*adekf::SO3(omega*deltaT)).conjugate()*gravity;
                             //adekf::viz::LinePlot::plotVector(ekf_true.mu.position, "True velocity",1000,"xyz");
                             //adekf::viz::LinePlot::plotVector(ekf_true.mu.position, "True velocity stride",1000,"xyz",3);
                             //adekf::viz::LinePlot::plotVector(ekf_world.mu.position, "Position World EKF",3000,"xyz");
                             //adekf::viz::LinePlot::plotVector(ekf_body.mu.velocity, "Body World EKF",3000,"xyz");
                             ekf_true.predictWithNonAdditiveNoise(dynamic_model_body, cov, acc, omega, deltaT);
                             //ekf_body.predictWithNonAdditiveNoise(dynamic_model_body, cov, acc + random_acc, omega + random_omega, deltaT);
                             //ekf_world.predictWithNonAdditiveNoise(dynamic_model_world, cov, acc + random_acc, omega + random_omega, deltaT);
                             ekf_gravity.predictWithNonAdditiveNoise(dynamic_model_gravity, cov, acc + random_acc, omega + random_omega, deltaT);
                             //LOG_STREAM << ekf_gravity.mu.gravity LOG_END
                             random_vel.poll();
                             //ekf_body.update(measurement_model_body, Eigen::Matrix3d::Identity() * 0.1, (ekf_true.mu.velocity + random_vel).eval());
                             //ekf_world.update(measurement_model_world, Eigen::Matrix3d::Identity() * 0.1 , (ekf_true.mu.velocity + random_vel).eval());
                             ekf_gravity.updateWithNullSpaceConstraint(measurement_model_body, Eigen::Matrix3d::Identity() * 0.1, (ekf_true.mu.velocity + random_vel).eval(),nullspace);
                             ekf_gravity.updateWithNullSpaceConstraint(consistency_model,Eigen::Matrix2d::Identity()*0.1,ekf_gravity.mu.gravity,nullspace);
                             //if(cnt++ % 100 ==0)
                             //ekf_gravity.sigma=erase_sigma*ekf_gravity.sigma*erase_sigma+start_cov;
                             adekf::viz::LinePlot::plotVector(ekf_gravity.sigma.diagonal().segment<3>(0),"Covariance of Orientation GRavity Filter",2000,"xyz",200);
                         }
                     }
    );

    adekf::viz::runGuis();
    std::cout << "Waiting for thread to stop ..." << std::endl;
    loop.join();
    return EXIT_SUCCESS;
}

