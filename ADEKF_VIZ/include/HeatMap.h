//
// Created by tomlucas on 24.03.20.
//

#ifndef ADEKF_HEATMAP_H
#define ADEKF_HEATMAP_H

#include <Eigen/Core>
#include <memory>
#include <jkqtplotter/jkqtplotter.h>
#include <jkqtplotter/graphs/jkqtpimage.h>
namespace adekf::viz{
    /**
     * Class to plot HeatMaps of covariance matrices
     *
     * First call adekf::viz::initGuis() (recommend) or
     * create a QApplication
     *
     * call displayCovariance once on an estimator or a double * to the matrix data
     * The class does not take ownership of passed pointers
     *
     * call adekf::viz::runGuis() to update the plots (recommended) or
     * call updatePlots() from the gui Thread
     *
     * In the end call adekf::viz::finishGuis() (recommended) or
     * disposePlots()
     */
    class HeatMap {
        //All Heatmap Plots
        inline static std::list<std::shared_ptr<JKQTPlotter> >plots;
    public:
        /**
         * Display the given DOFxDOF matrix in a heatmap
         * @param data pointer to the data.
         * @param DOF the size of the matrix
         * @param title The title of the plot
         * @param min The minimal expected covariance value
         * @param max The maximal expected covariance value
         */
        static void displayCovariance(double * data, int DOF, const char *  title, double  min=-2., double  max=10.) ;

        /**
         * Display the covariance of the given estimator
         * @tparam Estimator class of the passed estimator (auto deduced)
         * @param estimator pointer to the estimator
         * @param title The title of the plot
         * @param min The minimal expected covariance value
         * @param max The maximal expected covariance value
         */
        template<class Estimator>
        static void displayCovariance(Estimator * estimator,const char *  title, double  min=-2., double  max=10.){
            displayCovariance(estimator->sigma.data(),estimator->sigma.rows(),title,min,max);
        }
        /**
         * Update all Heatmaps
         *
         * Needs to be called repitetively in the GUI Thread
         */
        static void updatePlots();

        /**
         * Deletes all plots
         */
        static void disposePlots();
    };
}

#endif //ADEKF_HEATMAP_H
