//
// Created by tomlucas on 01.04.20.
//
#include "LinePlot.h"
#include <ADEKF/ADEKFUtils.h>

namespace adekf::viz {

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
    LinePlot::plotVector(const Eigen::VectorXd &vector, const char *title, size_t buffer_size, const char *legend) {
        auto it = plot_data.find(title);
        if (it == plot_data.end()) {
            plot_data.emplace(title, VectorRingBuffer(buffer_size, vector));
            size_t size = vector.rows(); // to not copy vector
            ioService.post([=]() { createPlot(size, title, buffer_size, legend); });
        } else {
            it->second.append(vector);
        }
    }
}

