add_executable(slam_and_orientation testADEKF.cpp)
target_link_libraries(slam_and_orientation PUBLIC ADEKF)
file(GLOB_RECURSE DATASETS Datasets/*)
file(COPY ${DATASETS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
