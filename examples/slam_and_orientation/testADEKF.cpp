#define EIGEN_STACK_ALLOCATION_LIMIT 0

#include <iostream>
#include <fstream>
#include "ADEKF.h"
#include "ManifoldCreator.h"

#include "types/SO3.h"

#include "ukfom/ukf.hpp"
#include "mtk/types/SOn.hpp"

#include <chrono>
#include <vector>

using namespace Eigen;

using namespace adekf;

class CSVRow
{
public:
    std::string const& operator[](std::size_t index) const
    {
        return m_data[index];
    }
    std::size_t size() const
    {
        return m_data.size();
    }
    void readNextRow(std::istream& str)
    {
        std::string         line;
        std::getline(str, line);

        std::stringstream   lineStream(line);
        std::string         cell;

        m_data.clear();
        while(std::getline(lineStream, cell, ';'))
        {
            m_data.emplace_back(cell);
        }
        // This checks for a trailing comma with no data after it.
        if (!lineStream && cell.empty())
        {
            // If there was a trailing comma then add an empty element.
            m_data.emplace_back("");
        }
    }
private:
    std::vector<std::string> m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class IMUData
{
public:
    unsigned long timestamp;
    Matrix<double, 3, 1> acc;
    Matrix<double, 3, 1> gyro;

    explicit IMUData(CSVRow& data)
    {
        timestamp = std::stoul(data[1]);
        acc << std::stod(data[2]), std::stod(data[3]), std::stod(data[4]);
        gyro << std::stod(data[5]), std::stod(data[6]), std::stod(data[7]);
    }
};

class TrackerData
{
public:
    unsigned long timestamp;
    Matrix<double, 3, 3> rot;

    explicit TrackerData(CSVRow& data)
    {
        timestamp = std::stoul(data[1]);
        rot.col(0) << std::stod(data[2]), std::stod(data[3]), std::stod(data[4]);
        rot.col(1) << std::stod(data[5]), std::stod(data[6]), std::stod(data[7]);
        rot.col(2) << std::stod(data[8]), std::stod(data[9]), std::stod(data[10]);
    }
};

auto imuDynMTK(const MTK::SO3<double> &state, const MTK::vect<3> &gyro) {
    return MTK::SO3<double>::exp(gyro * 0.01) * state;
}

auto trackMeasMTK(const MTK::SO3<double> &state) {
    return state;
}

void testHanddataset(bool withPrecision, unsigned repetitions) {
    std::cout << "----------BEGIN TEST HAND----------" << std::endl;
    std::ifstream imufile("0_IMU.csv");
    std::ifstream trackerfile("0_Tracker.csv");

    CSVRow row;

    std::vector<IMUData> imuDataVector;
    std::vector<TrackerData> trackerDataVector;

    bool first = true;

    while(imufile >> row)
    {
        if(!first) {
            imuDataVector.emplace_back(IMUData(row));
        }
        first = false;
    }

    first = true;

    while(trackerfile >> row)
    {
        if(!first) {
            trackerDataVector.emplace_back(TrackerData(row));
        }
        first = false;
    }

    ADEKF ekfwj(SO3d(trackerDataVector[0].rot), Matrix3d::Identity());

    ADEKF ekf(SO3d(trackerDataVector[0].rot), Matrix3d::Identity());

    ukfom::ukf<MTK::SO3<double>> ukf(MTK::SO3<double>(trackerDataVector[0].rot), Matrix3d::Identity());


    auto epoch = std::chrono::high_resolution_clock::now();
    for (unsigned j = 0; j < (withPrecision ? 1 : repetitions); j++) {
        unsigned trackerIdx = 0;
        for (const IMUData& data : imuDataVector) {
            TrackerData tData = trackerDataVector[trackerIdx];

            ekfwj.predictWithJacobian([](auto state, auto u){state = state + (0.01 * u);}, [](auto, auto u){Vector3d rot = -u*0.01;return AngleAxisd(rot.norm(), rot.normalized()).toRotationMatrix();}, Matrix3d::Identity(), data.gyro);

            if (tData.timestamp <= data.timestamp) {
                ekfwj.updateManifoldWithJacobian([](auto state){return state;}, [](auto){return Matrix3d::Identity();}, [](auto, auto &diff){Vector3d rot = -diff;return AngleAxisd(rot.norm(), rot.normalized()).toRotationMatrix();}, Matrix3d::Identity() * 0.1, SO3d(tData.rot));

                if (trackerIdx + 1 < trackerDataVector.size())
                    trackerIdx++;
            }
        }
    }
    auto now = std::chrono::high_resolution_clock::now();
    auto mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - epoch).count();
    std::cout << "ekfwj: " << mseconds << " ms" << std::endl;

    epoch = std::chrono::high_resolution_clock::now();
    for (unsigned j = 0; j < (withPrecision ? 1 : repetitions); j++) {
        unsigned trackerIdx = 0;
        for (const IMUData& data : imuDataVector) {
            TrackerData tData = trackerDataVector[trackerIdx];

            ekf.predict([](auto state, auto u){state = state + (0.01 * u);}, Matrix3d::Identity(), data.gyro);

            if (tData.timestamp <= data.timestamp) {
                ekf.update([](auto state){return state;}, Matrix3d::Identity() * 0.1, SO3d(tData.rot));

                if (trackerIdx + 1 < trackerDataVector.size())
                    trackerIdx++;
            }
        }
    }
    now = std::chrono::high_resolution_clock::now();
    mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - epoch).count();
    std::cout << "ekf: " << mseconds << " ms" << std::endl;

    epoch = std::chrono::high_resolution_clock::now();
    for (unsigned j = 0; j < (withPrecision ? 1 : repetitions); j++) {
        unsigned trackerIdx = 0;
        for (const IMUData& data : imuDataVector) {
            TrackerData tData = trackerDataVector[trackerIdx];

            ukf.predict(std::bind(imuDynMTK, std::placeholders::_1, data.gyro), [&] { return Matrix3d::Identity(); });

            if (tData.timestamp <= data.timestamp) {
                ukf.update(MTK::SO3<double>(tData.rot), trackMeasMTK, [&] { return Matrix3d::Identity() * 0.1; });

                if (trackerIdx + 1 < trackerDataVector.size())
                    trackerIdx++;
            }
        }
    }
    now = std::chrono::high_resolution_clock::now();
    mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - epoch).count();
    std::cout << "ukf: " << mseconds << " ms" << std::endl;
    std::cout << "----------END TEST HAND----------" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Meas
{
public:
    unsigned frameNumber;
    Vector2d pos;

    Matrix2d cov;

    unsigned id;

    explicit Meas(CSVRow& data)
    {
        frameNumber = std::stoul(data[0]);

        pos << std::stod(data[1]), std::stod(data[2]);

        cov << std::stod(data[3]), std::stod(data[4]),
                std::stod(data[4]), std::stod(data[5]);

        id = std::stoul(data[6]);
    }
};

class Step
{
public:
    Vector2d trans;
    double rot;
    double offset_array[3];

    Matrix3d cov;

    std::vector<Meas> landmarks;

    explicit Step(CSVRow& data)
    {
        trans << std::stod(data[1]), std::stod(data[2]);
        rot = std::stod(data[3]);
        offset_array[0] = std::stod(data[1]);
        offset_array[1] = std::stod(data[2]);
        offset_array[2] = std::stod(data[3]);

        cov << std::stod(data[4]), std::stod(data[5]), std::stod(data[7]),
                std::stod(data[5]), std::stod(data[6]), std::stod(data[8]),
                std::stod(data[7]), std::stod(data[8]), std::stod(data[9]);
    };

    void addMeas(const Meas& meas) {
        landmarks.emplace_back(meas);
    }
};

//constexpr unsigned MaxLandmarks = 62;
constexpr unsigned MaxLandmarks = 62;
//constexpr unsigned MaxLandmarks = 589;
constexpr unsigned StateSize =  3+(MaxLandmarks*2);

//constexpr unsigned MaxSteps = 195;
constexpr unsigned MaxSteps = 195;
//constexpr unsigned MaxSteps = 3297;

template<typename T = double>
using State = Matrix<T, StateSize,1>;
using Cov = Matrix<double, -1,-1>;

using StateMTK = MTK::vect<StateSize>;

struct InitLandmark {
    template<typename T>
    void operator()(State<T> &state, unsigned idx) {
        Rotation2D<T> rot(state(2));
        state.template segment<2>(idx) = (rot * state.template segment<2>(idx)) + state.template head<2>();
    }
} initLand;

struct InitLandmarkJacobian {
    template<typename T>
    Matrix<double, -1,-1> operator()(const State<T> &state, unsigned idx) {
        Matrix<double, -1, -1> jacobian = Matrix<double, -1, -1>::Identity(StateSize, StateSize);
        auto _cos = cos(state(2));
        auto _sin = sin(state(2));
        jacobian(idx,0) = 1;
        jacobian(idx+1,1) = 1;
        jacobian(idx,2) = (-state(idx+1)*_cos)-(state(idx)*_sin);
        jacobian(idx+1,2) = (state(idx)*_cos)-(state(idx+1)*_sin);
        jacobian(idx,idx) = _cos;
        jacobian(idx,idx+1) = -_sin;
        jacobian(idx+1,idx) = _sin;
        jacobian(idx+1,idx+1) = _cos;
        return jacobian;
    }
} initLandJac;

struct InitLandmarkMTK {
    auto operator()(const StateMTK &state, Meas &landmark) {
        Rotation2Dd rot(state(2));
        unsigned idx = 3 + (2*(landmark.id - 1));
        StateMTK result = state;
        result.template segment<2>(idx) = (rot * state.template segment<2>(idx)) + state.template head<2>();
        return result;
    }
} initLandMTK;

struct StepDynamics {
    template<typename T>
    void operator()(State<T> &state, const Step &step) {
        Rotation2D<T> rot(state(2));
        state.template head<2>() += rot * step.trans;
        state(2) += step.rot;
    }
} stepDyn;

struct StepDynamicsJacobian {
    template<typename T>
    Matrix<double, -1,-1> operator()(const State<T> &state, const Step &step) {
        Matrix<double, -1, -1> jacobian = Matrix<double, -1, -1>::Identity(StateSize, StateSize);
        T _sin = sin(state(2));
        T _cos = cos(state(2));
        jacobian(0,2) = (-step.trans(1)*_cos)-(step.trans(0)*_sin);
        jacobian(1,2) = (step.trans(0)*_cos)-(step.trans(1)*_sin);
        return jacobian;
    }
} stepDynJac;

struct StepDynamicsMTK {
    auto operator()(const StateMTK &state, const Step &step) {
        StateMTK result = state;
        Rotation2Dd rot(state(2));
        result.template head<2>() += rot * step.trans;
        result(2) += step.rot;
        return result;
    }
} stepDynMTK;

struct MeasureLandmark {
    template<typename T>
    auto operator()(const State<T> &state, const unsigned idx) {
        Rotation2D<T> rot(state(2));
        return rot.inverse() * (state.template segment<2>(idx) - state.template head<2>());
    }
} measLand;

struct MeasureLandmarkJacobian {
    template<typename T>
    Matrix<double, -1,-1> operator()(const State<T> &state, unsigned idx) {
        Matrix<double, -1, -1> jacobian = Matrix<double, -1, -1>::Zero(2,StateSize);
        auto _cos = cos(state(2));
        auto _sin = sin(state(2));
        jacobian(0,0) = -_cos;
        jacobian(1,0) = _sin;
        jacobian(0,1) = -_sin;
        jacobian(1,1) = -_cos;

        jacobian(0,2) = (state(idx+1)*_cos) - (state(1)*_cos) + (state(0)*_sin) - (state(idx)*_sin);
        jacobian(1,2) = (state(0)*_cos) - (state(idx)*_cos) + (state(1)*_sin) - (state(idx+1)*_sin);

        jacobian(0,idx) = _cos;
        jacobian(1,idx) = -_sin;
        jacobian(0,idx+1) = _sin;
        jacobian(1,idx+1) = _cos;
        return jacobian;
    }
} measLandJac;

struct MeasureLandmarkMTK {
    auto operator()(const StateMTK &state, const unsigned id) {
        Rotation2Dd rot(state(2));
        unsigned idx = 3 + ((id-1) *2);
        return rot.inverse() * (state.segment<2>(idx) - state.template head<2>());
    }
} measLandMTK;

void testSLAM(bool enableLandmarks, bool log) {
    std::cout << "----------BEGIN TEST SLAM----------" << std::endl;
    std::ifstream stepfile("dlr_odometry.csv");
    std::ifstream measfile("dlr_landmarks.csv");

    std::ofstream slam_ekf_pos;
    std::ofstream slam_ekf_map;

    std::ofstream slam_adekf_pos;
    std::ofstream slam_adekf_map;

    std::ofstream slam_ukf_pos;
    std::ofstream slam_ukf_map;

    std::ofstream slam_mat_pos;
    std::ofstream slam_mat_map;

    if(log) {
        if(enableLandmarks) {
            slam_ekf_pos.open("SLAM/slam_ekf_pos");
            slam_adekf_pos.open("SLAM/slam_adekf_pos");
            slam_ukf_pos.open("SLAM/slam_ukf_pos");

            slam_ekf_map.open("SLAM/slam_ekf_map");
            slam_adekf_map.open("SLAM/slam_adekf_map");
            slam_ukf_map.open("SLAM/slam_ukf_map");
        } else {
            slam_ekf_pos.open("SLAM/slam_ekf_odometry");
            slam_adekf_pos.open("SLAM/slam_adekf_odometry");
            slam_ukf_pos.open("SLAM/slam_ukf_odometry");
        }
    }

    CSVRow row;

    std::vector<Step> steps;
    std::vector<Meas> measurements;

    bool first = true;

    while(stepfile >> row)
    {
        if(!first) {
            steps.emplace_back(Step(row));
        }
        first = false;
    }

    first = true;

    while(measfile >> row)
    {
        if(!first) {
            measurements.emplace_back(Meas(row));
        }
        first = false;
    }

    for(const Meas& m : measurements) {
        if(enableLandmarks && m.id < 700)
            steps[m.frameNumber - 1].addMeas(m);
    }

    bool seen_landmark[MaxLandmarks];
    std::fill(seen_landmark,seen_landmark+MaxLandmarks,false);

    ADEKF ekf(State<double>::Zero(), Cov::Zero(StateSize, StateSize));

    ADEKF ekfwj(State<double>::Zero(), Cov::Zero(StateSize, StateSize));

    ukfom::ukf<StateMTK> ukf(StateMTK::Zero(), Cov::Zero(StateSize,StateSize));

    auto epoch = std::chrono::high_resolution_clock::now();

    for(unsigned i = 0; i < MaxSteps; i++) {
        const Step& s = steps[i];

        //std::cout << i+1 << "/" << MaxSteps << std::endl;

        Cov cov = Cov::Zero(StateSize, StateSize);
        cov.topLeftCorner<3,3>() = s.cov;
        ekfwj.predictWithJacobian(stepDyn,stepDynJac, cov, s);
        ekfwj.mu(2) = fmod(ekfwj.mu(2), M_PI * 2);
        if(ekfwj.mu(2) < double(0))
            ekfwj.mu(2) += M_PI * 2;
        for(const Meas& m : s.landmarks) {
            unsigned idx = 3 + (2 * (m.id - 1));
            if(seen_landmark[m.id - 1]) {
                ekfwj.updateWithJacobian(measLand,measLandJac,m.cov,m.pos,idx);
            } else {
                ekfwj.mu.segment<2>(idx) = m.pos;
                ekfwj.sigma.block<2,2>(idx, idx) = m.cov;
                ekfwj.predictWithJacobian(initLand,initLandJac,cov,idx);
                seen_landmark[m.id - 1] = true;
            }
        }
        if(log)
            slam_ekf_pos << ekfwj.mu[0] << ";" << ekfwj.mu[1] << ";" << ekfwj.mu[2] << std::endl;
    }

    auto now = std::chrono::high_resolution_clock::now();
    auto mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - epoch).count();
    std::cout << "ekf with jacobians: " << mseconds << " ms" << std::endl;

    std::fill(seen_landmark,seen_landmark+MaxLandmarks,false);

    epoch = std::chrono::high_resolution_clock::now();


    for(unsigned i = 0; i < MaxSteps; i++) {
        const Step& s = steps[i];

        //std::cout << i+1 << "/" << MaxSteps << std::endl;

        Cov cov = Cov::Zero(StateSize, StateSize);
        cov.topLeftCorner<3,3>() = s.cov;
        ekf.predict(stepDyn, cov, s);
        ekf.mu(2) = fmod(ekf.mu(2), M_PI * 2);
        if(ekf.mu(2) < double(0))
            ekf.mu(2) += M_PI * 2;

        for(const Meas& m : s.landmarks) {
            unsigned idx = 3 + (2 * (m.id - 1));
            if(seen_landmark[m.id - 1]) {
                ekf.update(measLand,m.cov,m.pos,idx);
            } else {
                ekf.mu.segment<2>(idx) = m.pos;
                ekf.sigma.block<2,2>(idx, idx) = m.cov;
                ekf.predict(initLand,cov,idx);
                seen_landmark[m.id - 1] = true;
            }
        }
        if(log)
            slam_adekf_pos << ekf.mu[0] << ";" << ekf.mu[1] << ";" << ekf.mu[2] << std::endl;
    }

    now = std::chrono::high_resolution_clock::now();
    mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - epoch).count();
    std::cout << "ekf: " << mseconds << " ms" << std::endl;

    std::fill(seen_landmark,seen_landmark+MaxLandmarks,false);

    epoch = std::chrono::high_resolution_clock::now();
    for(unsigned i = 0; i < MaxSteps; i++) {
        const Step& s = steps[i];

        Cov cov = Cov::Zero(StateSize, StateSize);
        cov.topLeftCorner<3,3>() = s.cov;
        ukf.predict(std::bind(stepDynMTK, std::placeholders::_1, s), [&]{return cov;});
        ukf.mu_(2) = fmod(ukf.mu_(2), M_PI * 2);
        if(ukf.mu_(2) < double(0))
            ukf.mu_(2) += M_PI * 2;

        for(const Meas& m : s.landmarks) {
            if(seen_landmark[m.id - 1]) {
                ukf.update(m.pos,std::bind(measLandMTK, std::placeholders::_1, m.id),m.cov);
            } else {
                unsigned idx = 3 + (2 * (m.id - 1));
                ukf.mu_.segment<2>(idx) = m.pos;
                ukf.sigma_.block<2,2>(idx, idx) = m.cov;
                ukf.predict(std::bind(initLandMTK, std::placeholders::_1, m),[&]{return cov;});
                seen_landmark[m.id - 1] = true;
            }
        }
        if(log)
            slam_ukf_pos << ukf.mu_[0] << ";" << ukf.mu_[1] << ";" << ukf.mu_[2] << std::endl;
    }
    now = std::chrono::high_resolution_clock::now();
    mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - epoch).count();
    std::cout << "ukf: " << mseconds << " ms" << std::endl;

    std::cout.precision(17);
    std::cout << ekfwj.mu.head<3>().format(IOFormat(FullPrecision)) << std::endl << std::endl;
    std::cout << ekf.mu.head<3>().format(IOFormat(FullPrecision)) << std::endl << std::endl;
    std::cout << ukf.mu_.head<3>().format(IOFormat(FullPrecision)) << std::endl;
    std::cout << "----------END TEST SLAM----------" << std::endl;

    if(log && enableLandmarks) {
        for(unsigned idx = 3; idx < StateSize; idx = idx + 2) {
            slam_ekf_map << ekfwj.mu[idx] << ";" << ekfwj.mu[idx+1] << std::endl;
            slam_adekf_map << ekf.mu[idx] << ";" << ekf.mu[idx+1] << std::endl;
            slam_ukf_map << ukf.mu_[idx] << ";" << ukf.mu_[idx+1] << std::endl;
        }
    }

    if(log) {
        if(enableLandmarks) {
            slam_ekf_pos.close();
            slam_adekf_pos.close();
            slam_ukf_pos.close();

            slam_ekf_map.close();
            slam_adekf_map.close();
            slam_ukf_map.close();
        } else {
            slam_ekf_pos.close();
            slam_adekf_pos.close();
            slam_ukf_pos.close();
        }
    }

}

int main() {
    testHanddataset(false, 1);
    testSLAM(true,false);
}
