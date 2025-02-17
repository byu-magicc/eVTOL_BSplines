#This is the top level CMakeLists.txt file which puts together everything else in the submodules

cmake_minimum_required(VERSION 3.16)
project(eVTOL_BSplines
        VERSION 0.0.1
        LANGUAGES CXX)


SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED 1)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)


FIND_PACKAGE(GTest REQUIRED)
FIND_PACKAGE(Eigen3 3.3 REQUIRED)

link_directories(${Eigen_INCLUDE_DIRS})
include_directories(${GTEST_INCLUDE_DIRS})

add_subdirectory(source/submodules/path_generator/PathObjectivesAndConstraints)
add_subdirectory(source/submodules/trajectory_generator/trajectory_generation/constraint_functions/TrajectoryConstraintsCCode)


enable_testing()