project(ADEKF_EXAMPLES)
cmake_minimum_required(VERSION 3.10)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(ADEKF CONFIG REQUIRED)
find_package(Boost 1.65.1 COMPONENTS system  REQUIRED)

# Configure project for usage of Qt5
find_package(Qt5 COMPONENTS Core Gui Widgets PrintSupport Svg Xml  REQUIRED)

# include JKQTPlotter
find_package(JKQTCommonSharedLib REQUIRED)
find_package(JKQTMathTextSharedLib REQUIRED)
find_package(JKQTPlotterSharedLib REQUIRED)

find_package(Threads REQUIRED)
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
  message("Skipping CylinderExample: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

# include all components
add_executable(CylinderExample MACOSX_BUNDLE odometry.cpp LinePlot.cpp Heatmap.cpp)
target_include_directories(CylinderExample PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(CylinderExample PUBLIC ${VTK_LIBRARIES} ${Boost_LIBRARIES}  ADEKF::ADEKF ${CMAKE_THREAD_LIBS_INIT}  Qt5::Core Qt5::Widgets Qt5::Gui Qt5::PrintSupport Qt5::Svg Qt5::Xml JKQTPlotterSharedLib )
# vtk_module_autoinit is needed
#target_link_libraries(qplot  Qt5::Core
#        Qt5::Gui
#        Qt5::Widgets Qt5::PrintSupport)

target_compile_features(CylinderExample PUBLIC  cxx_std_17 )
vtk_module_autoinit(
        TARGETS CylinderExample
        MODULES ${VTK_LIBRARIES}
)
