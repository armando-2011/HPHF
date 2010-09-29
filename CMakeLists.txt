cmake_minimum_required(VERSION 2.8)
project (HPHF)
add_executable(HPHF main.cpp)


set (CMAKE_MODULE_PATH ${HPHF_SOURCE_DIR}/cmake/modules)
find_package(Eigen2 REQUIRED)

INCLUDE_DIRECTORIES(${EIGEN2_INCLUDE_DIR})
