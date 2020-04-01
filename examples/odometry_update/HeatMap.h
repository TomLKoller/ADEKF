//
// Created by tomlucas on 24.03.20.
//

#ifndef CYLINDEREXAMPLE_HEATMAP_H
#define CYLINDEREXAMPLE_HEATMAP_H

#include <Eigen/Core>
#include <memory>
#include <jkqtplotter/jkqtplotter.h>
#include <jkqtplotter/graphs/jkqtpimage.h>
namespace adekf::viz{
    class CovarianceReader{
    public:
        virtual ~CovarianceReader(){}
        virtual void updatePlot()=0;
    };

    template<class EstimatorType>
    class GenericCovarianceReader : public CovarianceReader{
          EstimatorType *estimator;
          std::shared_ptr<JKQTPlotter> plot;
    public:
        void updatePlot(){
            plot->redrawPlot();
          }
        GenericCovarianceReader(EstimatorType * estimator, const char  * title, double min, double max) :estimator(estimator),plot(std::make_shared<JKQTPlotter>()){
            int DOF=estimator->sigma.rows();
            plot->getPlotter()->setUseAntiAliasingForGraphs(false); // nicer (but slower) plotting
            plot->getPlotter()->setUseAntiAliasingForSystem(false); // nicer (but slower) plotting
            plot->getPlotter()->setUseAntiAliasingForText(false); // nicer (but slower) text rendering
            JKQTPDatastore* ds=plot->getDatastore();
            size_t imageC=ds->addImageColumn(estimator->sigma.data(),DOF,DOF,QString(title));
            JKQTPColumnMathImage* graph=new JKQTPColumnMathImage(plot.get());
            graph->setTitle("");
            // image column with the data
            graph->setImageColumn(imageC);
            // set size of the data (the datastore does not contain this info, as it only manages 1D columns of data and this is used to assume a row-major ordering
            graph->setNx(DOF);
            graph->setNy(DOF);
            // where does the image start in the plot, given in plot-axis-coordinates (bottom-left corner)
            graph->setX(0);
            graph->setY(0);
            // width and height of the image in plot-axis-coordinates
            graph->setWidth(DOF);
            graph->setHeight(DOF);
            // color-map is "MATLAB"
            graph->setPalette(JKQTPMathImageMATLAB);
            // get coordinate axis of color-bar and set its label
            graph->getColorBarRightAxis()->setAxisLabel("Covariance");
            // determine min/max of data automatically and use it to set the range of the color-scale
            graph->setAutoImageRange(false);
            graph->setImageMin(min);
            graph->setImageMax(max);

            // 5. add the graphs to the plot, so it is actually displayed
            plot->addGraph(graph);

            plot->getYAxis()->setInverted(true);
            // 6. set axis labels
            plot->getXAxis()->setAxisLabel("");
            plot->getYAxis()->setAxisLabel("");


            // 7. fix axis and plot aspect ratio to 1
            plot->getPlotter()->setMaintainAspectRatio(true);
            plot->getPlotter()->setMaintainAxisAspectRatio(true);

            // 8 autoscale the plot so the graph is contained
            plot->zoomToFit();

            // show plotter and make it a decent size
            plot->show();
            plot->resize(600,600);
            plot->setWindowTitle(title);

        }
    };

    class HeatMap {
        inline static std::list<CovarianceReader *> covarianceReaders;
    public:

        template<class EstimatorType>
       static void displayCovariance(EstimatorType *  estimator, const char *  title, double  min=-2., double  max=10.) {
            covarianceReaders.push_back(new GenericCovarianceReader<EstimatorType>(estimator,title,min,max));
        }
        static void updatePlots(){
            for(auto reader: covarianceReaders){
               reader->updatePlot();
            }
        }

        static void disposePlots(){
            for(auto reader: covarianceReaders){
                delete reader;
            }
            covarianceReaders.clear();
        }
    };
}

#endif //CYLINDEREXAMPLE_HEATMAP_H
