cmake_minimum_required (VERSION 3.0)
project(ADEKF)

add_library(ADEKF INTERFACE)

option(ADEKF_WITH_BOOST "Boost enables Helper Macros to define Manifolds" ON)
option(ADEKF_BUILD_TESTS "Build Tests" ON)


target_include_directories(ADEKF INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/common>
  $<INSTALL_INTERFACE:include>)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_include_directories(ADEKF INTERFACE ${EIGEN3_INCLUDE_DIR})

if (ADEKF_WITH_BOOST)
    target_include_directories(ADEKF INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/boost_utils>
            $<INSTALL_INTERFACE:include>)
    find_package(Boost 1.65 REQUIRED)
    target_link_libraries(ADEKF INTERFACE Boost::boost)
endif ()

export(TARGETS ADEKF FILE ADEKFConfig.cmake)

if (ADEKF_BUILD_TESTS)
    add_executable(testADEKF test/testADEKF.cpp ${MATLAB_SOURCES})
    target_link_libraries(testADEKF ADEKF)
    target_compile_options(testADEKF PRIVATE -Wall -Wextra -pedantic -march=native)
    set_property(TARGET testADEKF PROPERTY CXX_STANDARD 17)
endif ()

