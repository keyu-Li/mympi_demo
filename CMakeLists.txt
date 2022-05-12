cmake_minimum_required(VERSION 2.8)
project(mympi_demo C)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories("src")
#add_subdirectory("src")

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER  mpicc)
set(CMAKE_C_STANDARD 99)

add_executable(mympi_demo main.c src/helloFromProcess.c src/helloPro.c)

# openMP 配置
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

