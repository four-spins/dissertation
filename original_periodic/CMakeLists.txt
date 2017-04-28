cmake_minimum_required (VERSION 2.8)
PROJECT(GoniPlaq)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})

add_definitions(-std=c++11 -DDLIB_ISO_CPP_ONLY )
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -DVERBOSE_WEIGHTS")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -march=native")

set (DLIB_ROOT "${PROJECT_SOURCE_DIR}/../external_libraries/dlib")

include_directories ( 
  include,
  "${DLIB_ROOT}"
  )

add_executable( goni3d_rec_muca.x src/main.cpp )

