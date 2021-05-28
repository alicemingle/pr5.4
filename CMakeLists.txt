cmake_minimum_required(VERSION 3.1)

project(zadacha5)

set(PROG_SRC main_5.cpp)
set(PROG_HDR CImg.h)

add_executable(${PROJECT_NAME} ${PROG_SRC} ${PROG_HDR})

configure_file(conf.txt conf.txt COPYONLY)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lX11 -lpthread")
