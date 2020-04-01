#include <ADEKF/ADEKF.h>
#include <ADEKF/ManifoldCreator.h>
#include <ADEKF/types/SO3.h>
#include <thread>
#include "adekf_viz.h"

ADEKF_MANIFOLD(Pose3d, ((adekf::SO3, orientation)), (3, position))

ADEKF_MANIFOLD(Position, ((adekf::SO3, orientation)), (3, position), (3, velocity))

ADEKF_MANIFOLD(Old_Position,((adekf::SO3, orientation)) , (3, position), (3, velocity), (3, old_position))

int main(int argc, char * argv[]) {


    //Simple Example with no orientation
    adekf::viz::initGuis(argc,argv);
    Eigen::Vector3d acc{1, 0, 0};
    Eigen::Vector3d omega{0,0,0.1};
    Eigen::Vector3d base_vel{1, 0, 0};
    Eigen::Vector3d vel{base_vel};
    Eigen::Vector3d pos{20,0,0};
    Eigen::Vector3d gravity{0,0,0};

    adekf::SO3d orient;
    adekf::StateInfo<adekf::SO3d> test;
    adekf::ADEKF ekf(Position(orient,pos, vel), Eigen::Matrix<double, 9, 9>::Identity());
    adekf::ADEKF ekf_true(Position(orient,pos, vel), Eigen::Matrix<double, 9, 9>::Identity());
    adekf::ADEKF ekf_old(Old_Position(orient,pos, vel, pos),
                         Eigen::Matrix<double, 12, 12>::Identity());
    adekf::viz::PoseRenderer::displayPose(&ekf, "Salmon");
    adekf::viz::PoseRenderer::displayPose(&ekf_true, "Green");
    adekf::viz::PoseRenderer::displayPose(&ekf_old, "Blue");

    adekf::viz::HeatMap::displayCovariance(&ekf_true, "Ground truth Filter");
    adekf::viz::HeatMap::displayCovariance(&ekf, "Normal Filter");
    adekf::viz::HeatMap::displayCovariance(&ekf_old, "Old state Filter");
    double deltaT = 0.02;
    auto dynamic_model = [&gravity](auto &state, auto noise, auto acc, auto omega, double deltaT) {
        state.orientation=state.orientation+(omega+NOISE(0,3))*deltaT;
        auto world_acc=state.orientation*(acc+NOISE(3,3))-gravity;
        state.position += state.orientation*state.velocity * deltaT +0.5* world_acc * deltaT*deltaT;
        state.velocity +=state.orientation.conjugate()*world_acc * deltaT-omega.cross(state.velocity)*deltaT;
    };
    auto dynamic_model_old = [dynamic_model](auto &state, auto noise, auto acc, auto omega, double deltaT) {
        state.old_position = state.position;
        dynamic_model(state,noise, acc, omega, deltaT);
    };
    auto measurement_model = [](auto state) { return state.velocity; };
    auto measurement_model_old = [](auto state) { return state.orientation.conjugate()*(state.position - state.old_position); };
    Eigen::Matrix<double, 6, 6> cov=cov.Identity() *0.1;
    Eigen::Vector3d random_acc,random_omega,random_vel;
    //adekf::viz::LinePlot::plotVector(ekf_true.mu.position, "acc",1000,"xyz");

    std::thread loop([&]() {
                         while (!adekf::viz::PoseRenderer::isDone()) {
                             random_acc= random_acc.Random()*0.1;
                             random_omega=random_omega.Random()*0.001;
                             acc=omega.cross(base_vel)+ekf_true.mu.orientation.conjugate()*gravity;
                             adekf::viz::LinePlot::plotVector(ekf_true.mu.position, "acc",1000,"xyz");
                             Eigen::Vector3d old_position=ekf_true.mu.position;
                             ekf_true.predictWithNonAdditiveNoise(dynamic_model, cov, acc,omega, deltaT);
                             ekf.predictWithNonAdditiveNoise(dynamic_model, cov, acc + random_acc, omega+random_omega, deltaT);
                             ekf_old.predictWithNonAdditiveNoise(dynamic_model_old, cov, acc + random_acc, omega+random_omega, deltaT);
                             random_vel= random_vel.Random() * 0.1;
                             ekf.update(measurement_model, Eigen::Matrix3d::Identity() * 0.1, (ekf_true.mu.velocity+random_vel).eval());
                             ekf_old.update(measurement_model_old, Eigen::Matrix3d::Identity() * 0.1*deltaT*deltaT, (ekf_true.mu.orientation.conjugate()*(ekf_true.mu.position-old_position)+random_vel*deltaT).eval());
                         }
                     }
    );

    adekf::viz::runGuis();
    std::cout << "Waiting for thread to stop ..." << std::endl;
    loop.join();
    return EXIT_SUCCESS;
}

