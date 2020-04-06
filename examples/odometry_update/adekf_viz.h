//
// Created by tomlucas on 24.03.20.
//

#ifndef ADEKF_ADEKF_VIZ_H
#define ADEKF_ADEKF_VIZ_H

#include "PoseRenderer.h"
#include "HeatMap.h"
#include "LinePlot.h"
#include <thread>
#include <memory>
#include <QApplication>

#include <boost/asio/io_service.hpp>

namespace adekf::viz {
    inline std::shared_ptr<QApplication> qwidget;

    /**
     * Initialize all GUIs (2D and 3D plots)
     * @param argc the argc argument of the main function
     * @param argv the argv argument of the main function
     */
    void initGuis(int &argc, char *argv[]) {
        qwidget = std::make_shared<QApplication>(argc, argv);
        PoseRenderer::initGui();
        LinePlot::initPlots();
    }
    /**
     * Disposes all GUIs and frees memory (2D and 3D)
     */
    void finishGuis() {
        PoseRenderer::disposeWindow();
        HeatMap::disposePlots();
        LinePlot::disposePlots();
    }
    /**
     * Runs all GUIs (2D and 3D)
     *
     * Has to be called at the main thread
     * Also calls finishGuis() when the 3D window is closed.
     * @see finishGuis
     */
    void runGuis() {
        LOG_STREAM << "Starting GUIs" LOG_END
        LinePlot::startUpdateThread();
        while (!PoseRenderer::isDone()) {
            PoseRenderer::updateWindow();
            HeatMap::updatePlots();
            LinePlot::ioService.poll();
            qwidget->processEvents();
// std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        finishGuis();
    }

    /**
     * Plots the given vector in a line graph.
     *
     * Call this function whenever you want to add a new vector to the plot
     * The title has to be unique.
     *
     * Just passes argument to @see adekf::viz::LinePlot::plotVector
     * @param vector The vector to plot
     * @param title The title of the plot. Has to be a unique identifier
     * @param buffer_size The amount of vectors of the past to be displayed
     * @param legend The legend entries of the vector axis. requires a char for each dimension of the vector
     * @param stride Only each strideth vector is shown
     */
    void plotVector(const Eigen::VectorXd &vector, const char *title, size_t buffer_size, const char *legend) {
        LinePlot::plotVector(vector, title, buffer_size, legend);
    }

    /**
     * Display the covariance of the given estimator
     *
     * Just passes arguments to @see adekf::viz::HeatMap::displayCovariance
     * @tparam Estimator class of the passed estimator (auto deduced)
     * @param estimator pointer to the estimator
     * @param title The title of the plot
     * @param min The minimal expected covariance value
     * @param max The maximal expected covariance value
     */
    template<class EstimatorType>
    static void displayCovariance(EstimatorType *estimator, const char *title, double min = -2., double max = 10.) {
        HeatMap::displayCovariance(estimator, title, min, max);
    }

    template<class EstimatorType>
    static void displayPosition(EstimatorType *estimator, const char *color) {
        PoseRenderer::displayPosition(estimator, color);
    }

    template<class EstimatorType>
    static void displayPose(EstimatorType *estimator, const char *color) {
        PoseRenderer::displayPose(estimator, color);
    }


}

#endif //ADEKF_ADEKF_VIZ_H
