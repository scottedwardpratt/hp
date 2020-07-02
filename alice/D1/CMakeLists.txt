cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(b3d_fromhydro CXX)
project(hydro2uds CXX)

add_custom_target(extern_b3d_alice COMMAND make WORKING_DIRECTORY ../software/b3d2)
add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ../../../commonutils/software)

#include(${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/CMakeLists.txt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/include ${CMAKE_CURRENT_SOURCE_DIR}/../software/hydrobalance/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/include")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/lib ${CMAKE_CURRENT_SOURCE_DIR}/../software/hydrobalance/lib ${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/lib)
message("Libraries: b3d_alice ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )


find_library(${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/lib b3d_alice)
find_library(${CMAKE_CURRENT_SOURCE_DIR}/../software/hydrobalance/lib hydro2uds_alice)
find_library(${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/lib commonutils)


set(SOURCE_FILES ./b3d_fromhydro.cc)
SET_SOURCE_FILES_PROPERTIES(
      b3d_fromhydro.cc PROPERTIES OBJECT_DEPENDS
    ../software/b3d2/lib/libb3d_alice.a
 )
add_executable(b3d_fromhydro ${SOURCE_FILES})
add_dependencies(b3d_fromhydro extern_b3d_alice extern_commonutils)
target_link_libraries(b3d_fromhydro b3d_alice commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

set(SOURCE_FILES ./hydro2uds.cc)
SET_SOURCE_FILES_PROPERTIES(
      hydro2uds.cc PROPERTIES OBJECT_DEPENDS
    ../software/b3d2/lib/libb3d_alice.a
 )
 SET_SOURCE_FILES_PROPERTIES(
       hydro2uds.cc PROPERTIES OBJECT_DEPENDS
     ../software/hydrobalance/lib/libhydro2uds_alice.a
  )
add_executable(hydro2uds ${SOURCE_FILES})
add_dependencies(hydro2uds extern_b3d_alice extern_commonutils)
target_link_libraries(hydro2uds hydro2uds_alice b3d_alice commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
