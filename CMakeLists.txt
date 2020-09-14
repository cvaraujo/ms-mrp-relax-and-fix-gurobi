cmake_minimum_required(VERSION 3.5)
project(MaxService)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/gurobi902/linux64/include /opt/boost_1_73_0)

link_directories(/opt/gurobi902/linux64/lib /opt/boost_1_73_0)

set(SOURCE_FILES main.cpp headers/Graph.h src/Graph.cpp headers/Include.h src/Model.cpp headers/Model.h src/Arc.cpp headers/Arc.h)

add_executable(MaxService ${SOURCE_FILES})

target_link_libraries(MaxService gurobi_c++ gurobi90 pthread m)
