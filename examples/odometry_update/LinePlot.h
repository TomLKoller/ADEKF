//
// Created by tomlucas on 24.03.20.
//

#ifndef ADEKF_LINEPLOT_H
#define ADEKF_LINEPLOT_H

#include <Eigen/Core>
#include <memory>
#include <map>
#include "jkqtplotter/jkqtplotter.h"
#include "jkqtplotter/graphs/jkqtpscatter.h"
#include <functional>
#include <ADEKF/ADEKFUtils.h>

#include <boost/asio/io_service.hpp>
#include <memory>
namespace adekf::viz {

    /**
     * Class to Buffer Vectors in  a RingBuffer
     * With this, the user can store the last n vectors to display it without increasing the memory usage in each time step
     */
    class VectorRingBuffer {
        //! contains a double array for each of the vector dimensions
        std::vector<std::shared_ptr<double[]>> buffer;
        //! maximum size of the buffer, current element, counter to implement the stride
        size_t buffer_size, current, stride_counter;
        //the stride (only each strideth vector is shown)
        const size_t stride;
    public:

/**
 * Constructor
 * @param buffer_size the maximum size of the buffer (how many past vectors will be shown)
 * @param first_vector  the first vector to show to determine the size of the vector
 * @param stride only each strideth vector will be shown 1 to show all
 */
        VectorRingBuffer(size_t buffer_size, const Eigen::VectorXd &first_vector, size_t stride);
/**
 * Appends an vector to the buffers.
 * Splits the vector into several double buffers (required for plotting)
 * @param vector the vector to add
 */
        void append(const Eigen::VectorXd &vector);

/**
 * Get the buffer to the corresponding element of the vector
 * @param index  vector element index
 * @return the double * which contains the past values of the vector element
 */
        double *getBuffer(size_t index);
    };
    /**
     * Class to plot line graphs.
     * The class is not intended to be instantiated outside of it
     * all relevant data is stored inside the static variables.
     *
     * First call adekf::viz::initGuis() (recommended) or
     * -create a QApplication  qwidget
     * -Call initPlots() before you draw the first vector
     *
     *
     * Call plotVector() to push the vector to the graph. Graphs are identified by title so do not use the same title twice!
     *
     * to update the plots  call adekf::viz:runGuis() (recommended) or
     * -startUpdateThread() which will update the plots every 40ms
     * -LinePlot::ioService.poll(); The plots allow to be created and filled from other threads than the main thread. This call is created to properly create plots
     * -qwidget->processEvents()  which will catch the update events of the update thread and finally update the plots
     *
     * Call adekf::viz::finishGuis() (recommended) at the end or
     * -disposePlots()
     */
    class LinePlot : public JKQTPlotter {
    Q_OBJECT
        //Maps titles to data
        static inline std::map<const char *, VectorRingBuffer> plot_data;
        //List of all plots
        static inline std::list<std::shared_ptr<JKQTPlotter> > plots;
        //Object to enable QPlot timer signals
        static inline std::shared_ptr<LinePlot> class_object;
        /**
         * Sets up the QT part of the plot (not the data buffering)
         * @param size  the vector size
         * @param title the title of the plot
         * @param buffer_size maximum buffer size of the plot
         * @param legend legend entries of the lines
         */
        static void createPlot(size_t size, const char *title, size_t buffer_size, const char *legend);

        /**
         * Private constructor to disable object instantiation from outside
         */
        LinePlot() {
        }

    public:
        // boost io Service object to call functions in the gui thread instead of other threads
        static inline boost::asio::io_service ioService;

        /**
         * Initializes the Lineplot static object
         */
        static void initPlots();

        /**
         * Plots the given vector in a line graph.
         *
         * Call this function whenever you want to add a new vector to the plot
         * The title has to be unique
         * @param vector The vector to plot
         * @param title The title of the plot. Has to be a unique identifier
         * @param buffer_size The amount of vectors of the past to be displayed
         * @param legend The legend entries of the vector axis. requires a char for each dimension of the vector
         * @param stride Only each strideth vector is shown
         */
        static void plotVector(const Eigen::VectorXd &vector, const char *title, size_t buffer_size, const char *legend,
                               size_t stride = 1);
        /**
         * Clears all plots and frees memory
         */
        static void disposePlots();
        /**
         * Starts the update timer
         */
        static void startUpdateThread();

    public slots:
        /**
         * Signal routine to update all plots
         */
        void updatePlots();
    };
}

#endif //ADEKF_LINEPLOT_H
