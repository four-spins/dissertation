cmake_minimum_required (VERSION 2.8)
PROJECT(GoniPlaq)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})

add_definitions(-std=c++11 -DDLIB_ISO_CPP_ONLY )
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -DVERBOSE_WEIGHTS")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -march=native")

set (DLIB_ROOT "${PROJECT_SOURCE_DIR}/../external_libraries/dlib")
set (GZSTREAM_ROOT "${PROJECT_SOURCE_DIR}/../external_libraries/gzstream")

include_directories ( 
  include,
  "${DLIB_ROOT}"
  "${GZSTREAM_ROOT}"
  )

add_subdirectory(gzstream)
add_executable( reweight_timeseries.x src/reweight_timeseries.cpp )
target_link_libraries( reweight_timeseries.x gzstream )

