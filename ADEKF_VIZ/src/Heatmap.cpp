//
// Created by tomlucas on 03.04.20.
//


#include "HeatMap.h"

namespace adekf::viz{
    void HeatMap::displayCovariance(double * data, int DOF, const char *  title, double  min, double  max) {
        // Create plot
        auto plot=std::make_shared<JKQTPlotter>();
        plots.push_front(plot);
        // Disable anti aliasing for spped
        plot->getPlotter()->setUseAntiAliasingForGraphs(false); // nicer (but slower) plotting
        plot->getPlotter()->setUseAntiAliasingForSystem(false); // nicer (but slower) plotting
        plot->getPlotter()->setUseAntiAliasingForText(false); // nicer (but slower) text rendering
        JKQTPDatastore* ds=plot->getDatastore();
        //add matrix data
        size_t imageC=ds->addImageColumn(data,DOF,DOF,QString(title));
        JKQTPColumnMathImage* graph=new JKQTPColumnMathImage(plot.get());
        graph->setTitle(QString(title));
        // image column with the data
        graph->setImageColumn(imageC);
        // set size of the data
        graph->setNx(DOF);
        graph->setNy(DOF);
        // where does the image start in the plot, given in plot-axis-coordinates (top-left corner)
        graph->setX(0);
        graph->setY(0);
        // width and height of the image in plot-axis-coordinates
        graph->setWidth(DOF);
        graph->setHeight(DOF);
        // color-map is "MATLAB"
        graph->setColorPalette(JKQTPMathImageMATLAB);
        // get coordinate axis of color-bar and set its label
        graph->getColorBarRightAxis()->setAxisLabel("Covariance");
        // Set min max value of the colorscale
        graph->setAutoImageRange(false);
        graph->setImageMin(min);
        graph->setImageMax(max);

        //  add the graphs to the plot, so it is actually displayed
        plot->addGraph(graph);

        //Invert y-axis to have classic matrix view with the origin in the top left
        plot->getYAxis()->setInverted(true);
        //  set axis labels
        plot->getXAxis()->setAxisLabel("");
        plot->getYAxis()->setAxisLabel("");


        // fix axis and plot aspect ratio to 1
        plot->getPlotter()->setMaintainAspectRatio(true);
        plot->getPlotter()->setMaintainAxisAspectRatio(true);

        //  autoscale the plot so the graph is contained
        plot->zoomToFit();

        // show plotter and make it a decent size
        plot->show();
        plot->resize(600,600);
        plot->setWindowTitle(title);
    }



    void HeatMap::updatePlots(){
        for(auto plot: plots){
            plot->redrawPlot();
        }
    }

    void HeatMap::disposePlots(){
        plots.clear();
    }
}
