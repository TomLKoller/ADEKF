#pragma once

/**
 * @brief Create the necessary operator() for ceres and the getLocalParametrization Function
 * Enables easy portation of this manifold to ceres
 * 
 */
#ifdef MANIFOLD_WITH_CERES
#include <ceres/local_parameterization.h>   
#include <ceres/autodiff_local_parameterization.h>  
#define LOCAL_PARAMETRISATION(TYPE_NAME, GLOB_SIZE)   \
	static constexpr unsigned GLOBAL_SIZE=GLOB_SIZE; \
    template<typename T>                            \
    bool operator()(const T* x,const T* delta,T* x_plus_delta) const {                  \
        (TYPE_NAME <T>(x)+Eigen::Map<const Eigen::Matrix<T,DOF,1>>(delta)).toPointer(x_plus_delta);          \
        return true;                    \
    }                                   \
    static ceres::LocalParameterization * getLocalParameterization(){                   \
        return new ceres::AutoDiffLocalParameterization<TYPE_NAME,GLOBAL_SIZE,DOF>{};   \
    } 
 #else
   #define LOCAL_PARAMETRISATION(TYPE_NAME,GLOB_SIZE) static constexpr unsigned GLOBAL_SIZE=GLOB_SIZE;
#endif 
   