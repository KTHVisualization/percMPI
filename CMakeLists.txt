cmake_minimum_required(VERSION 3.0.0)
project(HelloWorld VERSION 0.0.0)

set(SOURCE main.cpp watershed.cpp)
add_executable(${PROJECT_NAME} ${SOURCE})

#set(CPACK_PROJECT_NAME ${PROJECT_NAME})
#set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
#include(CPack)
