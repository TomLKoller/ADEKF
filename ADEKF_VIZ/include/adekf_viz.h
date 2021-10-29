//
// Created by tomlucas on 24.03.20.
//

#ifndef ADEKF_ADEKF_VIZ_H
#define ADEKF_ADEKF_VIZ_H

#include "PoseRenderer.h"
#include "HeatMap.h"
#include <memory>
#include <QApplication>

#include <boost/asio/io_service.hpp>

namespace adekf::viz {
    inline QApplication * qwidget;

    /**
     * Initialize all GUIs (2D and 3D plots)
     *
     * @Issue apparently a call to this function changes your c-locale (this may intefere with calls  atof, possibly others)
     *  This behaviour may result in atof to use "," instead of the "." delimiter for floating point numbers
     *
     * @param argc the argc argument of the main function
     * @param argv the argv argument of the main function
     */
    void initGuis(int &argc, char *argv[]) ;
    /**
     * Disposes all GUIs and frees memory (2D and 3D)
     */
    void finishGuis() ;
    /**
     * Runs all GUIs (2D and 3D)
     *
     * Has to be called at the main thread
     * Also calls finishGuis() when the 3D window is closed.
     * @see finishGuis
     */
    void runGuis() ;


    /**
     * Runs a single update step of all GUIs
     * Exposed for own update loop
     * @see runGuis()
     */
     void updateAll();

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
    void plotVector(const Eigen::VectorXd &vector, const char *title, size_t buffer_size, const char *legend, size_t stride=1) ;


    /**
        * Plots the given content in a line graph at once.
        *
        * Call this function whenever you want to plot a whole matrix at once.
        * The title has to be unique
        *
        * Just passes argument to @see adekf::viz::LinePlot::plotMatrix
        * @param whole_matrix The content to plot
        * @param title The title of the plot. Has to be a unique identifier
        * @param legend The legend entries of the vector axis. requires a char for each row dimension of the content
        */
    void plotMatrix(const Eigen::MatrixXd &whole_matrix, const char *title,  const char *legend);

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


    /**
        * Display the estimated position of the given estimator as a box.
        *
        * Passthrough to @see PoseRenderer::displayPosition
        * Requires that the estimator has a field mu (estimated state) with a field position (estimated position)
        * Stores the pointer (no ownership) to the estimator and reads the position with a @see GenericPositionReader
        * Call @see initGuis() before
        * @tparam EstimatorType The type of the estimator
        * @param estimator The pointer to the estimator
        * @param color The color of the displayed box
        */
    template<class EstimatorType>
    static void displayPosition(EstimatorType *estimator, const char *color) {
        PoseRenderer::displayPosition(estimator, color);
    }


    /**
        * Display the estimated pose of the given estimator as a box.
        *
        * Passthrough to @see PoseRenderer::displayPose
        * Requires that the estimator has a field mu (estimated state) with a field position (estimated position) and  a field orientation (estimated orientation)
        * Stores the pointer (no ownership) to the estimator and reads the position with a @see GenericPoseReader
        * Call @see initGuis() before
        * @tparam EstimatorType The type of the estimator
        * @param estimator The pointer to the estimator
        * @param color The color of the displayed box
        */
    template<class EstimatorType>
    static void displayPose(EstimatorType *estimator, const char *color) {
        PoseRenderer::displayPose(estimator, color);
    }


    /**
     * Color palette to choose colors from by index up to 6
     */
    static inline const char * colorpalette[]={"red","blue","yellow","green","violet","orange","brown"};

}

#endif //ADEKF_ADEKF_VIZ_H
