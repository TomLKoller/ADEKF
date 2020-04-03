//
// Created by tomlucas on 01.04.20.
//
#include "LinePlot.h"
#include <ADEKF/ADEKFUtils.h>

namespace adekf::viz {


    VectorRingBuffer::VectorRingBuffer(size_t buffer_size, const Eigen::VectorXd &first_vector, size_t stride)
            : buffer_size(buffer_size), current(0), current_size(0), stride_counter(-1), stride(stride) {
        for (size_t i = 0; i < first_vector.rows(); i++) {
            buffer.push_back((double *) calloc(buffer_size, sizeof(double)));
            std::fill_n(buffer.back(), buffer_size, 0);
        }
        append(first_vector);
    }

    void VectorRingBuffer::append(const Eigen::VectorXd &vector) {
        if (stride > 1) {
            stride_counter = (stride_counter + 1) % stride;
            if (stride_counter !=0) {
                return;
            }
        }

        for (size_t i = 0; i < vector.rows(); i++) {
            buffer[i][current] = vector(i);
        }
        if (current_size < buffer_size)
            current_size++;
        current++;
        if (current > current_size)
            current = current % current_size;
    }

    double VectorRingBuffer::operator()(size_t element, size_t index) {
        return buffer[index][(element - current + buffer_size) % buffer_size];
    }

    double *VectorRingBuffer::getBuffer(size_t index) {
        return buffer[index];
    }

    size_t VectorRingBuffer::getBufferSize() {
        return buffer_size;
    }

    size_t VectorRingBuffer::getCurrentSize() {
        return current_size;
    }

    BufferReader::BufferReader(VectorRingBuffer &buffer, size_t index) : buffer(buffer), index(index) {
    }

    inline double BufferReader::operator()(size_t element) {
        return buffer(element, index);
    }


    void LinePlot::initPlots() {
        class_object = std::shared_ptr<LinePlot>(new LinePlot());
    }

    void LinePlot::updatePlots() {
        for (auto &plot: plots) {
            plot->zoomToFit(false, true);
            plot->redrawPlot();
        }
        startUpdateThread();
    }

    void LinePlot::createPlot(size_t size, const char *title, size_t buffer_size, const char *legend) {
        plots.push_back(std::make_shared<JKQTPlotter>());
        JKQTPlotter *plot = plots.back().get();
        JKQTPDatastore *ds = plot->getDatastore();
        plot->setWindowTitle(QString(title));
        size_t columnT = ds->addLinearColumn(buffer_size, 0, buffer_size, "t");
        auto test = BufferReader(plot_data.find(title)->second, 0);
        for (size_t i = 0; i < size; i++) {
            size_t columnV = ds->addColumn(plot_data.find(title)->second.getBuffer(i), buffer_size, QString(legend[i]));
            JKQTPXYLineGraph *graph1 = new JKQTPXYLineGraph(plot);
            graph1->setXColumn(columnT);
            graph1->setYColumn(columnV);
            graph1->setTitle(QString(legend[i]));
            graph1->setSymbolType(JKQTPNoSymbol);
            // 5. add the graph to the plot, so it is actually displayed
            plot->addGraph(graph1);
        }

        plot->getPlotter()->setUseAntiAliasingForGraphs(false);
        plot->getPlotter()->setUseAntiAliasingForSystem(false);
        plot->getPlotter()->setUseAntiAliasingForText(false);
        // 6. autoscale the plot so the graph is contained
        //plot->zoomToFit();
        plot->setX(0, buffer_size);
        plot->setY(-20, 20);
        // show plotter and make it a decent size
        plot->show();
        plot->resize(600, 400);
    }

    void LinePlot::disposePlots() {
        plot_data.clear();
        plots.clear();
    }

    void LinePlot::startUpdateThread() {
        QTimer::singleShot(40, class_object.get(), SLOT(updatePlots()));
    }

    void
    LinePlot::plotVector(const Eigen::VectorXd &vector, const char *title, size_t buffer_size, const char *legend,size_t stride) {
        auto it = plot_data.find(title);
        if (it == plot_data.end()) {
            plot_data.emplace(title, VectorRingBuffer(buffer_size, vector,stride));
            size_t size = vector.rows(); // to not copy vector
            ioService.post([=]() { createPlot(size, title, buffer_size, legend); });
        } else {
            it->second.append(vector);
        }
    }
}

