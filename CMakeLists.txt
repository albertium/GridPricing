cmake_minimum_required(VERSION 3.9)
project(HullWhite)

set(CMAKE_CXX_STANDARD 11)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(HullWhite main.cpp Curve.cpp Curve.h GridPricer.cpp GridPricer.h DiffusionPDE.cpp DiffusionPDE.h PricerHelper.cpp PricerHelper.h Options.cpp Options.h Types.h)
target_link_libraries (HullWhite Eigen3::Eigen)