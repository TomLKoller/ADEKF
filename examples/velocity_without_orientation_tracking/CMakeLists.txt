project(ADEKF_EXAMPLES)
cmake_minimum_required(VERSION 3.0)

find_package(ADEKF CONFIG REQUIRED)
find_package(Boost 1.37 REQUIRED)

add_executable(velocity_without_orientation_tracking velocity_without_orientation_tracking.cpp)
target_link_libraries(velocity_without_orientation_tracking PRIVATE ADEKF::ADEKF)

