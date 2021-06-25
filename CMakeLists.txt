cmake_minimum_required (VERSION 3.1)
project("ADEKF" VERSION 1.0.0
                DESCRIPTION "Automatic Differentiated Extended Kalman Filter (ADEKF) - This is a generic EKF Implementation that uses automatic differentiation to get rid of the need to define Jacobians."
                LANGUAGES CXX)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

option(${PROJECT_NAME}_BUILD_EXAMPLES "Build examples" OFF)
option(${PROJECT_NAME}_BUILD_TESTS "Build Tests" OFF)
#option(${PROJECT_NAME}_BUILD_RESEARCH "Build Research" OFF)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ADEKF>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost 1.37 REQUIRED)

target_include_directories(${PROJECT_NAME} INTERFACE ${EIGEN3_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} INTERFACE Boost::boost)

install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ADEKF DESTINATION include)


if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ADEKF_VIZ)
    message("Building Vizualisation library")
    add_subdirectory(ADEKF_VIZ/)
endif()

if (${PROJECT_NAME}_BUILD_EXAMPLES)
    message("Building Examples")
    add_subdirectory(examples/slam_and_orientation examples/slam_and_orientation/build)
    add_subdirectory(examples/velocity_without_orientation_tracking examples/velocity_without_orientation_tracking/build)
    #add_subdirectory(examples/testingField examples/testingField/build)
endif ()


#if (${PROJECT_NAME}_BUILD_RESEARCH)
#    message("Building Research")
    #add_subdirectory(research/odometry_update research/odometry_update/build)
#    add_subdirectory(research/IMM research/IMM/build)
    #add_subdirectory(research/DutchTrackCycling research/DutchTrackCycling/build)
#endif ()


if(${PROJECT_NAME}_BUILD_TESTS)
    message("Building Tests")
    add_subdirectory(tests)
endif()

