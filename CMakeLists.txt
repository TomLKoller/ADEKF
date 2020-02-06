cmake_minimum_required (VERSION 3.10)
project("ADEKF" VERSION 1.0.0
                DESCRIPTION "Automatic Differentiated Extended Kalman Filter (ADEKF) - This is a generic EKF Implementation that uses automatic differentiation to get rid of the need to define Jacobians."
                LANGUAGES CXX)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

option(${PROJECT_NAME}_BUILD_TESTS "Build Tests" ON)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ADEKF>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost 1.65 REQUIRED)

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

if (${PROJECT_NAME}_BUILD_TESTS)
    #add_subdirectory(examples/slam_and_orientation examples/slam_and_orientation/build)
    #add_subdirectory(examples/velocity_without_orientation_tracking examples/velocity_without_orientation_tracking/build)
    add_subdirectory(examples/testingField examples/testingField/build)
endif ()

