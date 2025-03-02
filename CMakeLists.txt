cmake_minimum_required(VERSION 3.10)
project(OptimLight LANGUAGES CXX)

# Optimization Library

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options for selecting the linear algebra library
option(USE_ARMADILLO "Use Armadillo for linear algebra" ON)
option(USE_EIGEN "Use Eigen for linear algebra" OFF)

# Include directories
# include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Add library
add_library(OptimLight ${SOURCES})

# Link libraries based on the selected option
if(USE_ARMADILLO)
    find_package(Armadillo REQUIRED)
    target_link_libraries(OptimLight PUBLIC armadillo)
    target_compile_definitions(OptimLight PUBLIC USE_ARMADILLO)
elseif(USE_EIGEN)
    find_package(Eigen3 REQUIRED)
    target_link_libraries(OptimLight PUBLIC Eigen3::Eigen)
    target_compile_definitions(OptimLight PUBLIC USE_EIGEN)
endif()

# Add example executable
# add_executable(example ${PROJECT_SOURCE_DIR}/examples/main.cpp)
# target_link_libraries(example PRIVATE OptimLight)