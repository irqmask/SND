project(SND)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(SDL2 REQUIRED)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/ CACHE STRING
    "Do not change this setting."
    FORCE)

add_subdirectory(sndtest)
