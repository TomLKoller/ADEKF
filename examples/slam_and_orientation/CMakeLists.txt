add_executable(slam_and_orientation testADEKF.cpp)
target_link_libraries(slam_and_orientation PUBLIC ADEKF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#Enable Optimizations
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

file(GLOB_RECURSE DATASETS Datasets/*)
file(COPY ${DATASETS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
