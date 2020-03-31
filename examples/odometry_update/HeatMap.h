//
// Created by tomlucas on 24.03.20.
//

#ifndef CYLINDEREXAMPLE_HEATMAP_H
#define CYLINDEREXAMPLE_HEATMAP_H

#include <Eigen/Core>
#include <memory>
namespace adekf::viz{
    class CovarianceReader{
    public:
        virtual ~CovarianceReader(){}
        virtual void updatePlot()=0;
    };

    template<class EstimatorType>
    class GenericCovarianceReader : public CovarianceReader{
          EstimatorType *estimator;
    public:
        void updatePlot(){
          }
        GenericCovarianceReader(EstimatorType * estimator, const char  * title) :estimator(estimator){


// set up the QCPColorMap:

        }
    };

    class HeatMap {
        inline static std::list<CovarianceReader *> covarianceReaders;
    public:

        template<class EstimatorType>
       static void displayCovariance(EstimatorType * estimator, const char * title) {
            covarianceReaders.push_back(new GenericCovarianceReader<EstimatorType>(estimator,title));
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
