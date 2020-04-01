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



}

#endif //CYLINDEREXAMPLE_ADEKF_VIZ_H
