//
// Created by tomlucas on 24.03.20.
//

#ifndef CYLINDEREXAMPLE_ADEKF_VIZ_H
#define CYLINDEREXAMPLE_ADEKF_VIZ_H
#include "PoseRenderer.h"
#include "HeatMap.h"
#include "LinePlot.h"
#include <thread>
#include <memory>
#include <QApplication>

#include <boost/asio/io_service.hpp>
namespace adekf::viz{
    inline std::shared_ptr<QApplication> qwidget;
    void initGuis(int & argc, char * argv[]){
        qwidget=std::make_shared<QApplication>(argc,argv);
        PoseRenderer::initGui();
        LinePlot::initPlots();
    }

    void finishGuis(){
        PoseRenderer::disposeWindow();
        HeatMap::disposePlots();
        LinePlot::disposePlots();
    }

    void runGuis(){
        std::cout << "Starting gui" << std::endl;
        LinePlot::startUpdateThread();
        while(!PoseRenderer::isDone()) {
            PoseRenderer::updateWindow();
            HeatMap::updatePlots();
            LinePlot::ioService.poll();
            qwidget->processEvents();
// std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        finishGuis();
    }

    void plotVector(const Eigen::VectorXd & vector, const char *title, size_t buffer_size, const char * legend){
        LinePlot::plotVector(vector,title,buffer_size,legend);
    }
    template<class EstimatorType>
    static void displayCovariance(EstimatorType *  estimator, const char *  title, double  min=-2., double  max=10.){
        HeatMap::displayCovariance(estimator,title,min,max);
    }
    template<class EstimatorType>
    static void displayPosition(EstimatorType * estimator, const char * color){
        PoseRenderer::displayPosition(estimator,color);
    }
    template<class EstimatorType>
    static void displayPose(EstimatorType * estimator, const char * color){
        PoseRenderer::displayPose(estimator,color);
    }


}

#endif //CYLINDEREXAMPLE_ADEKF_VIZ_H
