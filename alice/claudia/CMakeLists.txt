cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(claudia CXX)
project(bfnorm CXX)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_b3d_alice COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/hp/alice/software/b3d2)
add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/include ${CMAKE_CURRENT_SOURCE_DIR}/../software/hydrobalance/include ${GITHOME_MSU}/commonutils/software/include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/b3d2/include ${GITHOME_MSU}/commonutils/software/include")

link_directories(${GITHOME_MSU}/hp/alice/software/b3d2/lib ${GITHOME_MSU}/hp/alice/software/hydrobalance/lib ${GITHOME_MSU}/commonutils/software/lib)
message("Libraries: b3d_alice ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

find_library(${GITHOME_MSU}/hp/alice/software/b3d2/lib b3d_alice)
find_library(${GITHOME_MSU}/hp/alice/software/hydrobalance/lib hydro2uds_alice)
find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)

set(SOURCE_FILES ./claudia.cc)
SET_SOURCE_FILES_PROPERTIES(
      claudia.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a
 )
add_executable(claudia ${SOURCE_FILES})
add_dependencies(claudia extern_b3d_alice extern_commonutils)
target_link_libraries(claudia b3d_alice  commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

set(SOURCE_FILES ./bfnorm.cc)
SET_SOURCE_FILES_PROPERTIES(
      bfnorm.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a
 )
add_executable(bfnorm ${SOURCE_FILES})
add_dependencies(bfnorm extern_b3d_alice extern_commonutils)
target_link_libraries(bfnorm b3d_alice  commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})