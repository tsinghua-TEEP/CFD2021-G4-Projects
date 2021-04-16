cmake_minimum_required(VERSION 3.1)
project(CFD2021Projects)
set (CMAKE_CXX_STANDARD 11)

option(CFD2021PROJECTS_DEBUG "Output detailed results for debugging." OFF)

if(CFD2021PROJECTS_DEBUG)
   add_definitions(-D_DEBUG_)
endif()

include_directories(include)

aux_source_directory(. SRC)

source_group(SOURCE\ FILES FILES ${SRC})

add_executable(CFD2021Projects ${SRC}
