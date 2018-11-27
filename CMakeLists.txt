cmake_minimum_required(VERSION 3.0.0)
project(PercMPI VERSION 0.0.0)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(SOURCE
    src/main.cpp src/percolationloader.cpp
    src/datablock.cpp
    src/unionfind.cpp
    src/unionfindblock.cpp
    )
add_executable(${PROJECT_NAME} ${SOURCE})
target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES})
