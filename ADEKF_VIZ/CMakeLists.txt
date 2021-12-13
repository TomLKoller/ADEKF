cmake_minimum_required (VERSION 3.10)
project("ADEKF_VIZ" VERSION 1.0.0
                DESCRIPTION "Automatic Differentiated Extended Kalman Filter (ADEKF) Visualization - This is a library to visualize the estimates of the ADEKF."
                LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTORCC ON)
#set(CMAKE_AUTOUIC ON)


include(GNUInstallDirs)

FILE(GLOB SOURCES src/*.cpp)
FILE(GLOB HEADERS include/*.h)
message(${SOURCES})
find_package(Eigen3 REQUIRED)
add_library(${PROJECT_NAME} SHARED ${SOURCES})
target_compile_definitions(${PROJECT_NAME} PRIVATE ADEKF_VIZ_BUILD_LIB)

set_target_properties(${PROJECT_NAME} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1)

find_package(Boost 1.65.1 COMPONENTS system  REQUIRED)

# Configure project for usage of Qt5
find_package(Qt5 COMPONENTS Core Gui Widgets PrintSupport Svg Xml  REQUIRED)

# include JKQTPlotter
find_package(JKQTCommonSharedLib REQUIRED)
find_package(JKQTMathTextSharedLib REQUIRED)
find_package(JKQTPlotterSharedLib REQUIRED)

find_package(VTK COMPONENTS
        vtkCommonColor
        vtkCommonCore
        vtkCommonDataModel
        vtkFiltersSources
        vtkInteractionStyle
        vtkRenderingContextOpenGL2
        vtkRenderingCore
        vtkRenderingFreeType
        vtkRenderingGL2PSOpenGL2
        vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
    message("Skipping: ${VTK_NOT_FOUND_MESSAGE}")
    return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

find_package(Eigen3)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/ ${EIGEN3_INCLUDE_DIR})

# include all components
target_link_libraries(${PROJECT_NAME} PUBLIC ${VTK_LIBRARIES} ${Boost_LIBRARIES}  ${CMAKE_THREAD_LIBS_INIT}  Qt5::Core Qt5::Widgets Qt5::Gui Qt5::PrintSupport Qt5::Svg Qt5::Xml JKQTPlotterSharedLib )

target_compile_features(${PROJECT_NAME}  PUBLIC  cxx_std_17 )
vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)



install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
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


install(FILES ${HEADERS} DESTINATION include/ADEKF/viz)
