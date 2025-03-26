cmake_minimum_required(VERSION 3.10)
project(OptimLight  
  VERSION 0.1.0
  DESCRIPTION "OptimLight  is a manifold optimization library."
  HOMEPAGE_URL "https://github.com/xxxxx"
  LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# version
set (OPTIMLIGHT_MAJOR_VERSION "0")
set (OPTIMLIGHT_MINOR_VERSION "0")
set (OPTIMLIGHT_BUGFIX_VERSION "0")
set (OPTIMLIGHT_VERSION_STRING ${OPTIMLIGHT_MAJOR_VERSION}.${OPTIMLIGHT_MINOR_VERSION}.${OPTIMLIGHT_BUGFIX_VERSION})
message (STATUS "OptimLight  version ${OPTIMLIGHT_VERSION_STRING}")

list (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# Build options
option(USE_ARMADILLO "Use Armadillo for linear algebra" ON)
option(USE_EIGEN "Use Eigen for linear algebra" OFF)
option(BUILD_TESTING "Build OptimLight  unit tests" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)



# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# target_include_directories(OptimLight 
#     PUBLIC 
#         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#         $<INSTALL_INTERFACE:include>
#     PRIVATE
#         ${CMAKE_CURRENT_SOURCE_DIR}/src
# )

# Add source files
set (OPTIMLIGHT_SOURCES
    src/types.cpp
    # src/optimlight.cpp
)

# Create the main library
add_library(OptimLight ${OPTIMLIGHT_SOURCES})

# Linear algebra backend configuration
if(USE_ARMADILLO)
    find_package(Armadillo REQUIRED)
    target_link_libraries(OptimLight  PUBLIC armadillo)
    target_compile_definitions(OptimLight  PUBLIC USE_ARMADILLO)
elseif(USE_EIGEN)
    find_package(Eigen3 REQUIRED)
    target_link_libraries(OptimLight  PUBLIC Eigen3::Eigen)
    target_compile_definitions(OptimLight  PUBLIC USE_EIGEN)
endif()



# Testing configuration
if(BUILD_TESTING)
    include(CTest)
    enable_testing()
    find_package(GTest HINTS /usr/local/lib/ ${GTEST_DIR})
    if(NOT ${GTest_FOUND})
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG "release-1.12.1"
            GIT_SHALLOW TRUE
            GIT_PROGRESS TRUE)
        FetchContent_MakeAvailable(googletest)
    endif()
    add_subdirectory(tests)
endif()

# Installation configuration
include(GNUInstallDirs)
install(TARGETS OptimLight 
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

# Export targets
install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)