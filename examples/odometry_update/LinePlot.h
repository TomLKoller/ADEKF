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
        size_t buffer_size,current,current_size;
    public:
        VectorRingBuffer(size_t buffer_size,const Eigen::VectorXd & first_vector):buffer_size(buffer_size),current(0),current_size(0){
            for(size_t i=0; i < first_vector.rows(); i++) {
                buffer.push_back((double *)calloc(buffer_size, sizeof(double)));
                std::fill_n(buffer.back(),buffer_size,0);
            }
            append(first_vector);
        }
        void append(const Eigen::VectorXd &vector){
            for(size_t i=0; i < vector.rows(); i++) {
                buffer[i][current] = vector(i);
            }
            if(current_size < buffer_size)
                current_size++;
            current++;
            if(current > current_size)
                current=current % current_size;
        }
        double operator()(size_t element,size_t index){
            return buffer[index][(element-current+buffer_size)%buffer_size];
        }

        double * getBuffer(size_t index){
            return buffer[index];
        }

        size_t getBufferSize(){
            return buffer_size;
        }
        size_t getCurrentSize(){
            return current_size;
        }
    };
class BufferReader{
    VectorRingBuffer & buffer;
    size_t index;
public:
    BufferReader(VectorRingBuffer & buffer, size_t index) :buffer(buffer),index(index){
    }
    double operator()(size_t element){
            return buffer(element,index);
    }

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
        static void plotVector(const Eigen::VectorXd & vector, const char *title, size_t buffer_size, const char * legend) ;
        static void disposePlots();
        static void startUpdateThread();

    public slots:
        void updatePlots();
    };
}

#endif //ADEKF_LINEPLOT_H
