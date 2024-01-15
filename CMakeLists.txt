cmake_minimum_required(VERSION 3.10)

project(Double_Pendulum)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(DOUBLE_PENDULUM doublependulum.h doublependulum_Main.cpp)

