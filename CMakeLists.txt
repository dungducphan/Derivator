cmake_minimum_required(VERSION 3.23)
project(gradient_descent)

set(CMAKE_CXX_STANDARD 17)
include_directories(gradient_descent_lib)

add_library(Derivator gradient_descent_lib/Derivator.cpp gradient_descent_lib/Derivator.h)
add_executable(gradient_descent main.cpp)
target_link_libraries(gradient_descent Derivator)
