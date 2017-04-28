/******************************************************************************
 *
 * @file: config.hpp
 *
 * @date: 03/01/2013 04:32:52 PM (CET)
 *
 * @author: Marco Müller <muelma@gmail.com>
 *
 ******************************************************************************/
# ifndef CONFIG_HPP 
# define CONFIG_HPP 

// -----
// Simulation parameters: 
const int L = 8;                   // linear lattice size in all directions
const int V = L*L*L;                // volume of the cubic lattice (or number of spins)
const int measurements = 100000;    // number of measurements to write out
const int thermalization = 0;     // number of thermalization steps before calculating weights
const int measure_every = 1;        // measure every "measure_every" sweeps

// restrictions on the energy interval
//const double minimal_energy = -1.5*V;
const double minimal_energy = -3.0/2.0*((L+1)*(L+1)*(L+1) - (L+1)*(L+1));
const double maximal_energy = 0.5*minimal_energy; //0.0;
const double binsize = 2.0;

// desired parameters of recursive muca algorithm
//
// maximum number of multicanical iterations 
// (ignored, if tunnel_events are encountered first)
const int max_muca_iterations = 100;  
// maximum number of tunnel events (system travelling from minimal_energy 
// to maximal_energy and back again)
const int max_tunnel_events = 20;
// sweep in each iteration
const int sweeps_per_iteration = V;
// seed for the random number generator
const int seed = 341234;

# endif // CONFIG_HPP
