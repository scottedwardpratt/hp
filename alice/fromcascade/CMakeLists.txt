
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(b3d_fromcascade CXX)
project(balance_fromcascade CXX)

add_custom_target(extern_b3d_alice COMMAND make WORKING_DIRECTORY ../software/b3d2)
add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ../../../commonutils/software)

#include(${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/CMakeLists.txt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/include")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/lib ${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/lib)
message("Libraries: b3d_alice ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

find_library(${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/lib b3d_alice)
find_library(${CMAKE_CURRENT_SOURCE_DIR}/../../../commonutils/software/lib commonutils)

set(SOURCE_FILES ./b3d_fromcascade.cc)
SET_SOURCE_FILES_PROPERTIES(
      b3d_fromcascade.cc PROPERTIES OBJECT_DEPENDS
    ../software/b3d2/lib/libb3d_alice.a
 )
add_executable(b3d_fromcascade ${SOURCE_FILES})
add_dependencies(b3d_fromcascade extern_b3d_alice extern_commonutils)
target_link_libraries(b3d_fromcascade b3d_alice commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

set(SOURCE_FILES ./balance_fromcascade.cc)
SET_SOURCE_FILES_PROPERTIES(
      balance_fromcascade.cc PROPERTIES OBJECT_DEPENDS
    ../software/b3d2/lib/libb3d_alice.a
 )
add_executable(balance_fromcascade ${SOURCE_FILES})
add_dependencies(balance_fromcascade extern_b3d_alice extern_commonutils)
target_link_libraries(balance_fromcascade b3d_alice commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

