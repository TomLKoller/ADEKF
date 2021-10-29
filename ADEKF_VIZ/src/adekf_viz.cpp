//
// Created by tomlucas on 07.04.20.
//

#include "adekf_viz.h"
#include "LinePlot.h"
namespace adekf::viz{
    void initGuis(int &argc, char *argv[]) {
        qwidget =new QApplication(argc, argv);
        PoseRenderer::initGui();
    }
    void finishGuis() {
        PoseRenderer::disposeWindow();
        HeatMap::disposePlots();
        LinePlot::disposePlots();
    }

    void updateAll(){
        PoseRenderer::updateWindow();
        HeatMap::updatePlots();
        LinePlot::ioService.poll();
        LinePlot::ioService.reset();
        LinePlot::updatePlots();
        qwidget->processEvents();
    }
    void runGuis() {
        while (!PoseRenderer::isDone()) {
           updateAll();
        }
        finishGuis();
    }
    void plotVector(const Eigen::VectorXd &vector, const char *title, size_t buffer_size, const char *legend,size_t stride) {
        LinePlot::plotVector(vector, title, buffer_size, legend,stride);
    }

    void plotMatrix(const Eigen::MatrixXd &whole_matrix, const char *title,  const char *legend){
        LinePlot::plotMatrix(whole_matrix,title,legend);
    }

}