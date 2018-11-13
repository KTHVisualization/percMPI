cmake_minimum_required(VERSION 3.0.0)
project(PercMPI VERSION 0.0.0)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(SOURCE main.cpp watershed.cpp)
add_executable(${PROJECT_NAME} ${SOURCE})
target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES})
