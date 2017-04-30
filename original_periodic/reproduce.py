#!/usr/bin/env python
from __future__ import print_function
import os
import subprocess
import shutil
import random
import sys

#Ls          = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]   # linear dimension for a 'L**d' cubes
Ls    = range(4,28)   # linear dimension for a 'L**d' cubes
seeds = [ 877999, 793557, 887979, 105099, 545336, 684037, 607688, 595393, 795316, 965577, 484040, 271837, 604557, 978687, 38116 , 745300, 569718, 447785, 620240, 760296, 321213, 562871, 61145 , 850676 ]

work_dir = os.getcwd() + '/'

cmake_file_content =\
    """cmake_minimum_required (VERSION 2.8)
PROJECT(GoniPlaq)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})

add_definitions(-std=c++11 -DDLIB_ISO_CPP_ONLY )
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -DVERBOSE_WEIGHTS")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O2 -march=native")

set (DLIB_ROOT "${PROJECT_SOURCE_DIR}/../../external_libraries/dlib")

include_directories ( 
  include,
  "${DLIB_ROOT}"
  )

add_executable( goni3d_rec_muca.x  "${PROJECT_BINARY_DIR}/src/main.cpp" )
"""


for i, L in enumerate(Ls):
    config_hpp = \
    """// created via script (arbeitsamt.py)
# ifndef CONFIG_HPP 
# define CONFIG_HPP 

// Simulation parameters: 
const int L = {0:d};                // linear lattice size in all directions
const int V = L*L*L;                // volume of the cubic lattice (or number of spins)
const int measurements = 131072;    // number of measurements to write out
const int thermalization = L*L;     // number of thermalization steps before calculating weights
const int measure_every = V;        // measure every "measure_every" sweeps
const int layer_updates = 0;        // number of layer updates after each sweep
const int checkpoint_every = 1024;  // make checkpoints every 'this number' of MEASUREMENTS

// restrictions on the energy interval
const double minimal_energy = -3.0/2.0*V;
const double maximal_energy = 0.0;
const double binsize = 2.0;

// desired parameters of recursive muca algorithm
//
// maximum number of multicanical iterations 
// (ignored, if tunnel_events are encountered first)
const int max_muca_iterations = V;  
// maximum number of tunnel events (system travelling from minimal_energy 
// to maximal_energy and back again)
const int max_tunnel_events = 20;
// sweep in each iteration
const int sweeps_per_iteration = V;
// seed for the random number generator
const int seed = {1:d};

# endif // CONFIG_HPP
""".format(L, seeds[i])
    
    # preparing source
    path = work_dir + 'L%d/'%L
    os.makedirs(path)
    shutil.copytree( work_dir + 'src', path + "src")

    with open(path + 'CMakeLists.txt', 'w') as f:
      f.write(cmake_file_content)

    with open(path + "/src/config.hpp", 'w') as f:
        f.write(config_hpp)

    # building
    buildpath = path# + 'build'
    os.chdir(buildpath)
    os.system('cmake . -DCMAKE_BUILD_TYPE=release')
    subprocess.call(['make'], shell=True)

