cmake_minimum_required (VERSION 2.8)
PROJECT(Anisotropic_Ashkin_Teller)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})

add_definitions(-std=c++11 -DDSFMT_MEXP=19937 -DMUCA -DDLIB_ISO_CPP_ONLY)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -DVERBOSE_WEIGHTS")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -march=native")

set (DSFMTLIB_ROOT "${PROJECT_SOURCE_DIR}/../external_libraries/dSFMT-src-2.2.3")
set (DLIB_ROOT "${PROJECT_SOURCE_DIR}/../external_libraries/dlib-19.2")

include_directories (
  "${PROJECT_SOURCE_DIR}/include"
  "${DSFMTLIB_ROOT}"
  "${DLIB_ROOT}"
  )

add_library(dsfmt "${DSFMTLIB_ROOT}/dSFMT.c")
add_library(dlib "${DLIB_ROOT}/dlib/all/source.cpp")

add_executable( AniAshkTell.x  src/main.cpp )

target_link_libraries (AniAshkTell.x dsfmt dlib)
