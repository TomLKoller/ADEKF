//
// Created by tomlucas on 01.04.20.
//
#include "LinePlot.h"
#include <ADEKF/ADEKFUtils.h>

namespace adekf::viz {

    VectorRingBuffer::VectorRingBuffer(size_t buffer_size, const Eigen::VectorXd &first_vector, size_t stride)
            : buffer_size(buffer_size), current(0),  stride_counter(-1), stride(stride) {
        for (size_t i = 0; i < first_vector.rows(); i++) {
            //Create buffers and fill them with 0s
            buffer.push_back(std::shared_ptr<double[]>(new double[buffer_size]));
            std::fill_n(buffer.back().get(), buffer_size, 0);
        }
        append(first_vector);
    }

    VectorRingBuffer::VectorRingBuffer(const Eigen::MatrixXd & whole_matrix)
            : buffer_size(whole_matrix.cols()), current(0),  stride_counter(-1), stride(1) {
        for (size_t i = 0; i < whole_matrix.rows(); i++) {
            //Create buffers and fill them with 0s
            buffer.push_back(std::shared_ptr<double[]>(new double[buffer_size]));
        }
        setBuffer(whole_matrix);
    }

    void VectorRingBuffer::setBuffer(const Eigen::MatrixXd & whole_matrix){
        for (size_t i = 0; i < whole_matrix.rows(); i++) {
            (Eigen::Map<Eigen::VectorXd>(buffer[i].get(),whole_matrix.cols()))=whole_matrix.row(i);
        }
    }


    void VectorRingBuffer::append(const Eigen::VectorXd &vector) {
        //only add each strideth element
        if (stride > 1) {
            stride_counter = (stride_counter + 1) % stride;
            if (stride_counter !=0) {
                return;
            }
        }
        //add elements to buffers
        for (size_t i = 0; i < vector.rows(); i++) {
            buffer[i][current] = vector(i);
        }
        current++;
        if (current == buffer_size)
            current = current % buffer_size;
    }

    double *VectorRingBuffer::getBuffer(size_t index) {
        return buffer[index].get();
    }



    void LinePlot::updatePlots() {
        for (auto &plot: plots) {
            plot->zoomToFit(false, true);
            plot->redrawPlot();
        }
    }

    void LinePlot::createPlot(size_t size, const char *title, size_t buffer_size, const char *legend) {
        LOG_STREAM << "Creating plot " << title LOG_END
        auto plot = std::make_shared<JKQTPlotter>();
        plots.push_front(plot);
        JKQTPDatastore *ds = plot->getDatastore();
        plot->setWindowTitle(QString(title));
        //add time axis
        size_t columnT = ds->addLinearColumn(buffer_size, 0, buffer_size, "t");
        // add ones graph for each element of the vector
        for (size_t i = 0; i < size; i++) {
            size_t columnV = ds->addColumn(plot_data.find(title)->second.getBuffer(i), buffer_size, QString(legend[i]));
            JKQTPXYLineGraph *graph1 = new JKQTPXYLineGraph(plot.get());
            graph1->setXColumn(columnT);
            graph1->setYColumn(columnV);
            graph1->setTitle(QString(legend[i]));
            graph1->setSymbolType(JKQTPNoSymbol);
            plot->addGraph(graph1);
        }
        //disable antialiasing for faster graphics
        plot->getPlotter()->setUseAntiAliasingForGraphs(false);
        plot->getPlotter()->setUseAntiAliasingForSystem(false);
        plot->getPlotter()->setUseAntiAliasingForText(false);

        //set intervall of the time axis
        plot->setX(0, buffer_size);
        // show plotter and make it a decent size
        plot->show();
        plot->resize(600, 400);
    }

    void LinePlot::disposePlots() {
        plot_data.clear();
        plots.clear();
    }


    void
    LinePlot::plotVector(const Eigen::VectorXd &vector, const char *title, size_t buffer_size, const char *legend,size_t stride) {
        //Check if title already exists
        auto it = plot_data.find(title);
        //create new plot if not, append vector otherwise
        if (it == plot_data.end()) {
            plot_data.emplace(title, VectorRingBuffer(buffer_size, vector,stride));
            size_t size = vector.rows(); // to not copy vector
            ioService.post([=]() { createPlot(size, title, buffer_size, legend); });
        } else {
            it->second.append(vector);
        }
    }


    void LinePlot::plotMatrix(const Eigen::MatrixXd &whole_matrix, const char *title,  const char *legend){
        //Check if title already exists
        auto it = plot_data.find(title);
        //create new plot if not, append vector otherwise
        if (it == plot_data.end()) {
            plot_data.emplace(title, VectorRingBuffer(whole_matrix));
            size_t size = whole_matrix.rows(); // to not copy vector
            ioService.post([=]() { createPlot(size, title, whole_matrix.cols(), legend); });
        } else {
            it->second.setBuffer(whole_matrix);
        }
    }
}

