project(sndtest)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11) 

find_package(SDL2 REQUIRED)
if(WIN32)
    file(GLOB SDL_DLLS ${SDL2_LIBRARY_DIR}/*.dll)
endif(WIN32)

# add_compile_options(-Wall)
include_directories(${SDL2_INCLUDE_DIRS})

set(sndtest_SRC
    main.cpp
)

add_executable(sndtest ${sndtest_SRC})

target_link_libraries(sndtest ${SDL2_LIBRARIES})

if(WIN32)
    add_custom_command(TARGET sndtest POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${SDL_DLLS} ${EXECUTABLE_OUTPUT_PATH})
endif(WIN32)
