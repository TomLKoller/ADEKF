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
namespace adekf::viz {

class VectorRingBuffer{
        std::vector<double *> buffer;
        size_t buffer_size,current,current_size,stride_counter;
        const size_t stride ;
    public:
        VectorRingBuffer(size_t buffer_size,const Eigen::VectorXd & first_vector,size_t stride);
        void append(const Eigen::VectorXd &vector);
        double operator()(size_t element,size_t index);

        double * getBuffer(size_t index);

        size_t getBufferSize();
        size_t getCurrentSize();
    };
class BufferReader{
    VectorRingBuffer & buffer;
    size_t index;
public:
    BufferReader(VectorRingBuffer & buffer, size_t index);
    inline double operator()(size_t element);

};


    class LinePlot :public JKQTPlotter {
        Q_OBJECT
    static inline std::map<const char*, VectorRingBuffer > plot_data;
    static inline std::vector<std::shared_ptr<JKQTPlotter> > plots;
    static inline std::shared_ptr<LinePlot> class_object;
    static void createPlot(size_t size,const char *title, size_t buffer_size, const char * legend);
    LinePlot(){
    }
    public:
        static inline boost::asio::io_service ioService;
        static void initPlots();
        static void plotVector(const Eigen::VectorXd & vector, const char *title, size_t buffer_size, const char * legend,size_t stride=1) ;
        static void disposePlots();
        static void startUpdateThread();

    public slots:
        void updatePlots();
    };
}

#endif //ADEKF_LINEPLOT_H
