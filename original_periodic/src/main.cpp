/******************************************************************************
 *
 * @file: main.cpp
 *
 * @date: 02/27/2013 11:35:09 AM (CET)
 *
 * @author: Marco Müller <muelma@gmail.com>
 *
 ******************************************************************************/

// comments starting with //leg are legacy lines used for publications in 2014
# include <iostream>
# include <iomanip>
# include <vector>
# include <cassert>
# include <cmath>
# include <fstream>
# include <ctime>
# include <sstream>
# include <csignal>

// random number generator
# include "rngclass.mar.hpp"
// histograms and such
# include "discrete_function.hpp"
// type for logarithmic calculations
# include "logval.hpp"
// stores all the necessary parameters
# include "config.hpp" 

// signalling
volatile sig_atomic_t defer_signal; // for sensible code, defer_signal++ <code> defer_signal--
volatile sig_atomic_t signal_pending; // remembers a signal, if defer_signal != 0
void signal_callback_handler(int signum);

// histogram storage
void hists2file();

// helpers
logval ln_total_updates(0.0);     // total number of updates
logval ln_total_updates_for_weights(0.0); // total number of updates for weights
logval ln_accepted_updates(0.0);  // number of accepted updates

int current_m(0);  // current measurement
int current_m_skip(0); // current count for skipping measurements
int checkpoint_m(-1); // the identifier for the measuremnt when the checkpoint
                     // created

unsigned long checkpoints(0);
int current_checkpoint_counter(0);

int nr_updates_layer(0);
int nr_accepted_updates_layer(0);

double total_energy;              // current total energy of the system
rngclass rng(seed);

std::vector<int> state(V, +1);  // stores the current configuration
// lower and upper bound for histogram-like datatypes
const double low_bound  = minimal_energy - binsize/2.0;
const double high_bound = maximal_energy + binsize/2.0;
// stores the histogram of all measurements
discrete_function<logval> lnhisto            (low_bound, high_bound, binsize);
// stores the histogram of all sweeps
discrete_function<logval> lnhisto_all_sweeps (low_bound, high_bound, binsize);
// stores weights w(E) in p(E) /propto exp(w[E])
discrete_function<logval> weights            (low_bound, high_bound, binsize);
// helper storage for the full statistics for recursive estimation of weights
discrete_function<logval> statistics         (low_bound, high_bound, binsize);
// helper storage for the transition rates between consecutive weight - bins
discrete_function<logval> transition         (low_bound, high_bound, binsize);

// these two flags are used to register tunnel events
// they are altered by void update() and used in calculate_weights(),
// therefore globally used
bool touched_min_energy = false;
bool touched_max_energy = false;

# ifdef STORE_CARTESIAN
std::vector<int> spin_x(V, -1);  // stores the x-coordinate of the spin
std::vector<int> spin_y(V, -1);  // stores the y-coordinate of the spin
std::vector<int> spin_z(V, -1);  // stores the z-coordinate of the spin
# endif

// to collect some statistics on the runtime of parts of the simulation
time_t time_start, time_before_weights;
time_t time_before_measurement, time_after_measurement;
// this will be the overall logfile
std::ofstream logfile; 
// measurement file
std::ofstream meas_file;
// measurment of the full energy time series
std::ofstream energy_meas_file; 

void checkpoint();
void from_checkpoint(unsigned nr);

// creates the header for measurement files
std::string get_header(){
  std::stringstream header;
  header << "# recursive Multicanonical simulation\n";
  header << "# Gonihedric Ising Model, plaquette only\n";
  time_t rawtime; time( &rawtime );
  header << "# " << ctime( &rawtime );
  header << "#\n";
  header << "# random number generator = " << rng.get_type() << "\n";
  header << "# random number seed = " << rng.get_seed() << "\n";
  header << "# boundary conditions = periodic\n";
  header << "# beta = 0\n";
  header << "# grid dimensions = " << L << ", " << L << ", " << L << "\n";
  header << "# volume = " << L*L*L << "\n";
  header << "# E_min = " << minimal_energy << "\n";
  header << "# E_max = " << maximal_energy << "\n";
  header << "# number of measurements = " << measurements << "\n";
  header << "# measured every [this number] sweeps = " << measure_every << "\n";
  header << "# thermalization sweeps = " << thermalization << "\n";
  return header.str();
}

// creates statistics on simulation properties for output in measurement files
// this function is only to be called after measurements, when global statistics
// variables like ln_total_updates, ln_accepted_updates, times are set
std::string get_statistics(){
  std::stringstream statistics;
  statistics << "# total number of updates = " << ln_total_updates.toDouble() << "\n";
  statistics << "# number of accepted updates = " << ln_accepted_updates.toDouble() << "\n";
  logval frac = ln_accepted_updates/ln_total_updates;
  statistics << "# acceptance ratio = " << frac.toDouble() << "\n";
  statistics << "# Timings:\n";
  statistics << "# thermalization (seconds) = " << time_before_weights - time_start << "\n";
  statistics << "# weights calculation (seconds) = " << time_before_measurement - time_before_weights << "\n";
  statistics << "# measurements (seconds) = " << time_after_measurement - time_before_measurement << "\n";
  statistics << "# layer flip acceptance ratio = " << static_cast<double>(nr_accepted_updates_layer)/nr_updates_layer << "\n";
  statistics << "# layer flips = " << nr_accepted_updates_layer << "\n";
  statistics << "# layer flip attempts = " << nr_updates_layer << "\n";
  return statistics.str();
}

// we use an index for the storage of system configurations, so 
// 1) it is easier to iterate over spins
// 2) we do not need to draw random numbers for each dimension
//    when choosing spins randomly in the system
// 
// injective mapping of three cartesian coordinates to an 1D-index
inline int cartesian_to_index(int x, int y, int z){
  assert( x >= 0 && x < L );
  assert( y >= 0 && y < L );
  assert( z >= 0 && z < L );
  return x + L*y + L*L*z;
}
// calculate the x-coordinate from an index
inline int index_to_x(int index){
# ifdef STORE_CARTESIAN
  return spin_x[index];
# else
  return index % L;
# endif
}
// calculate the y-coordinate from an index
inline int index_to_y(int index){
# ifdef STORE_CARTESIAN
  return spin_y[index];
# else
  return (index/L) % L;
# endif
}
// calculate the z-coordinate from an index
inline int index_to_z(int index){
# ifdef STORE_CARTESIAN
  return spin_z[index];
# else
  return index/(L*L);
# endif
}

// delivers the value of the spin at position (x,y,z)
// this function is supposed to handle boundary conditions
inline int spin(int x, int y, int z){
  // periodic boundary conditions
  if( x < 0 ) x += L;
  if( y < 0 ) y += L;
  if( z < 0 ) z += L;
  if( x >= L ) x -= L;
  if( y >= L ) y -= L;
  if( z >= L ) z -= L;
  assert( x >= 0 && x < L );
  assert( y >= 0 && y < L );
  assert( z >= 0 && z < L );
  return state[cartesian_to_index(x,y,z)];
}
// get the local energy of a spin (that is the sum over all contributions
// in the Hamiltonian which contain this particular spin)
inline double local_energy_cartesian( int x, int y, int z ){
  double sum = 0;
  // plaquettes parallel to xy-plane
  sum += spin(x+1,y,z) * spin(x+1,y+1,z) * spin(x,y+1,z);
  sum += spin(x-1,y,z) * spin(x-1,y+1,z) * spin(x,y+1,z); 
  sum += spin(x-1,y,z) * spin(x-1,y-1,z) * spin(x,y-1,z); 
  sum += spin(x+1,y,z) * spin(x+1,y-1,z) * spin(x,y-1,z); 

  // plaquettes parallel to xz-plane
  sum += spin(x+1,y,z) * spin(x+1,y,z+1) * spin(x,y,z+1);
  sum += spin(x-1,y,z) * spin(x-1,y,z+1) * spin(x,y,z+1); 
  sum += spin(x-1,y,z) * spin(x-1,y,z-1) * spin(x,y,z-1); 
  sum += spin(x+1,y,z) * spin(x+1,y,z-1) * spin(x,y,z-1); 
  
  // plaquettes parallel to yz-plane
  sum += spin(x,y+1,z) * spin(x,y+1,z+1) * spin(x,y,z+1); 
  sum += spin(x,y-1,z) * spin(x,y-1,z-1) * spin(x,y,z-1);
  sum += spin(x,y-1,z) * spin(x,y-1,z+1) * spin(x,y,z+1);
  sum += spin(x,y+1,z) * spin(x,y+1,z-1) * spin(x,y,z-1);

  sum *= spin(x,y,z);
  return -0.5*sum;
}
inline double local_energy(int index){
  int x = index_to_x(index);
  int y = index_to_y(index);
  int z = index_to_z(index);
  return local_energy_cartesian(x,y,z);
}

inline double energy(){
  double sum = 0;
  for( int x = 0; x < L; x++ )
    for( int y = 0; y < L; y++ )
      for( int z = 0; z < L; z++ ) {
        // periodic boundary conditions: look only in "plus"-directions
        double sumi = 0;
        sumi += spin(x+1,y,z) * spin(x+1,y+1,z) * spin(x,y+1,z); // first plaquette parallel to xy-plane
        sumi += spin(x+1,y,z) * spin(x+1,y,z+1) * spin(x,y,z+1); // first plaquette parallel to xz-plane
        sumi += spin(x,y+1,z) * spin(x,y+1,z+1) * spin(x,y,z+1); // first plaquette parallel to yz-plane
        sumi *= spin(x,y,z);
        sum += sumi;
      }
  return -0.5*sum; 
}


bool test_indices_calculations() {
  bool succeed = true;
  for( int x = 0; x < L; x++ )
    for( int y = 0; y < L; y++ )
      for( int z = 0; z < L; z++ ) {
        int index = cartesian_to_index( x, y, z);
        int xb = index_to_x(index);
        int yb = index_to_y(index);
        int zb = index_to_z(index);
        if ( xb != x or yb != y or zb != z)
          succeed = false;
      }
  return succeed;
}
inline void flip_layer(int layer, int dir) {
  if ( dir == 0 ){
    for (int x=0; x<L;x++){
      for (int y=0; y<L;y++){
        state.at(cartesian_to_index(layer, x,y)) *= -1;
      }
    }
  }
  if ( dir == 1 ){
    for (int x=0; x<L;x++){
      for (int y=0; y<L;y++){
        state.at(cartesian_to_index(x,layer,y)) *= -1;
      }
    }
  }
  if ( dir == 3 ){
    for (int x=0; x<L;x++){
      for (int y=0; y<L;y++){
        state.at(cartesian_to_index(x,y,layer)) *= -1;
      }
    }
  }
}
//inline void update_layer(int layer, int dir){
//  assert(layer >= 0 && layer < L);
//  assert(dir >= 0 && dir < 3);
//  
//  nr_updates_layer++;
//  flip_layer(layer, dir);
//  nr_accepted_updates_layer++;
//
///*
//  const double old_energy = energy();
//  flip_layer(layer, dir);
//  const double new_energy = energy();
//  const double dE = new_energy - old_energy; 
//  if ( old_energy != new_energy ){
//    std::cout << " err " << std::endl;
//  }
//
//  // test whether we are in the energy interval of interest
//  if ( new_energy - minimal_energy > -1e-6 && new_energy - maximal_energy < 1e-6 ) {
//    // check for acceptance of the new state
//    const double r = rng.rand();
//    const double frac = exp(weights[new_energy].value() - weights[total_energy].value());
//    if ( r <= frac ) {
//      // accept
//      total_energy += dE;
//      nr_accepted_updates_layer++;
//    } 
//    else {
//      // revert to previous state
//      flip_layer(layer, dir);
//    }
//  }
//
//  assert(energy() == total_energy);
//*/
//}
inline void update(int site){
  const double dE = -2.0*local_energy(site);
  const double new_energy = total_energy + dE;
  ln_total_updates++;

  // test whether we are in the energy interval of interest
  if ( new_energy - minimal_energy > -1e-6 && new_energy - maximal_energy < 1e-6 ) {
    // check for acceptance of the new state
    const double r = rng.rand();
    const double frac = exp(weights[new_energy].value() - weights[total_energy].value());
//    if (std::abs(new_energy - total_energy) > 1e-15)
//      std::cout << "new: " << new_energy << ", old: " << total_energy << ", frac: " << std::scientific << frac << std::endl; 
    if ( r <= frac ) {
      state[site] *= -1;
      total_energy += dE;
      ln_accepted_updates++;
    } 
  }

  if ( new_energy <= minimal_energy ){
    touched_min_energy = true;
  } else if ( new_energy >= maximal_energy ){
    touched_max_energy = true;
  }
  
  assert(energy() == total_energy);
}
inline void sweep(){
  for ( int i = 0; i < V; i++ ) {
    // choose a spin, uniformly distributed
    double r = rng.rand(); // in [0,1)
    int site = floor( r * V );
    assert(site >= 0 and site < V);
    update(site);  
  }
//  for ( int i=0; i < layer_updates; i++){
//    int layer = floor( rng.rand() * L );
//    int dir = floor( rng.rand() * 3 );
//    update_layer(layer, dir);
//  }
}

inline void calculate_weights(){
  // recursive determination of multicanonical weights
  //
  logval exponent; // kappa in the original paper of Janke
  logval fraction; // helper variable for better readability
  int tunnel_events = 0; // current number of tunnel events
  int muca_iterations = 0;
  // last registered minimal energy
  double last_energy_min = 1.0e100;
  double last_energy_max = -1.0e100;
  while ( tunnel_events < max_tunnel_events 
      && muca_iterations < max_muca_iterations ) {
    muca_iterations++;
# ifdef VERBOSE_WEIGHTS
    std::cout << "# iterations; tunnel events: " 
      << muca_iterations << "/" << max_muca_iterations << " ; "
      << tunnel_events << "/" << max_tunnel_events << std::endl;
# endif
    // clear histogram
    std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
    // reset current energy minima and maxima
    double cur_energy_min = 1.0e100;
    double cur_energy_max = -1.0e100;

//leg    for ( int m = 0; m < thermalization; m++ ) {
//leg      sweep();
//leg    }
    for ( int iteration = 0; iteration < sweeps_per_iteration; iteration++ ){
      // make a sweep
      for ( int step = 0; step < V; step++ ){
        double r = rng.rand(); // in [0,1)
        int site = floor( r * V );
        assert(site >= 0 and site < V);
        update(site);
        lnhisto[total_energy] += logval(1.0);

        cur_energy_min = std::min(total_energy, cur_energy_min);
        cur_energy_max = std::max(total_energy, cur_energy_max);

        if ( touched_min_energy && touched_max_energy ){
          tunnel_events++;
          // reset state
          touched_max_energy = false;
          touched_min_energy = false;
          // std::cout << "# registered tunnel event, tunnel_events = " << tunnel_events << "\n";
          if ( tunnel_events >= max_tunnel_events ){
            break;
          }
        }
//leg   } // sweep
        if ( tunnel_events >= max_tunnel_events ){
          break;
        }
      } // sweep //!leg
      for ( int i=0; i < layer_updates; i++){
        int layer = floor( rng.rand() * L );
        int dir = floor( rng.rand() * 3 );
        // update_layer(layer, dir);
        lnhisto[total_energy] += logval(1.0);
      }
    } // sweeps per iteration

    // break out of while before another update of weights is performed!
    if ( tunnel_events >= max_tunnel_events ){
      break;
    }

    // update of weights:
    last_energy_min = std::min(cur_energy_min, last_energy_min);
    last_energy_max = std::max(cur_energy_max, last_energy_max);
    if ( last_energy_min > minimal_energy ){
      size_t e_min_i = lnhisto.get_index( last_energy_min ); // index of the minimal energy found
      size_t e_max_i = lnhisto.get_index( last_energy_max ); // index of maximal energy found
      assert( e_max_i > e_min_i );
      double diff = weights[e_max_i].value() - weights[e_min_i].value();
      double a = diff/static_cast<double>(e_max_i - e_min_i);
      double n = weights[e_max_i].value() - a * static_cast<double>(e_max_i);
      for ( size_t i = 0; i < e_min_i; i++ ){
        weights[i] = tologvalExp( a* static_cast<double>(i) + n );
      }
    } // linear interpolation
    // accumulated statistics:
    for ( size_t i = 0; i < lnhisto.size() - 1; i++ ){
      // special case with no histogram entries so far
      if ( lnhisto[i] == logval(0.0) || lnhisto[i+1] == logval(0.0) ){
        exponent = logval(1.0);
        fraction = logval(1.0);
      } else {
        // no division by zero possible here
        fraction = lnhisto[i]*lnhisto[i+1]/(lnhisto[i]+lnhisto[i+1]);
        statistics[i] += fraction;        
        // here, division by zero can occur
        exponent = fraction/statistics[i];
        fraction = lnhisto[i] / lnhisto[i+1];
      }
      // check for NaNs
      assert( fraction == fraction );
      logval trans = weights[i+1]/weights[i];
      transition[i] = trans*pow(fraction.toDouble(), exponent.toDouble());
    } // accumulated statistics
    for ( size_t i = 0; i < lnhisto.size() - 1; i++ ){
      weights[i+1] = transition[i]*weights[i];
    }
    // weights are updated -> count current steps
    ln_total_updates_for_weights = ln_total_updates;
  } // while

  // TODO:
  // statistics -- how many iterations needed to store the weights (notice that while can be breaked out early)
  logfile << "# -- number of muca iterations: " << muca_iterations << " of " << max_muca_iterations << "\n";
  logfile << "# -- encountered tunnel events: " << tunnel_events << " of " <<  max_tunnel_events << "\n";

}

struct fuke_nuke {
  int fuke_nuke_abs;
  int fuke_nuke_sq;
};
inline fuke_nuke fuke_nuke_x(){
  /// calculates the fuke-nuke order parameter with fixed x-coordinates
  /// ie. \f[\sum_\mathrm{all\ yz-planes}\left\langle\left|\sum_
  /// \mathrm{single\ plane}\sigma_i\sigma_{i+\hat{e}_x}\right|\right\rangle \f]
  int full_sum_abs = 0;
  int full_sum_sq  = 0;
  // sum over yz-planes
  for ( int x = 0; x < L; x++ ){
    // sum over a single plane
    int sum_of_products = 0;
    for ( int y = 0; y < L; y++ ) {
      for ( int z = 0; z < L; z++ ) {
        sum_of_products += spin(x,y,z) * spin(x+1,y,z);
      }
    }
    full_sum_abs += std::abs( sum_of_products );
    full_sum_sq  += sum_of_products*sum_of_products;
  }
  fuke_nuke ret;
  ret.fuke_nuke_abs = full_sum_abs;
  ret.fuke_nuke_sq  = full_sum_sq;
  return ret;
}
struct fuke_nuke_struct {
  int abs_x;
  int sq_x;
  int abs_y;
  int sq_y;
  int abs_z;
  int sq_z;
  int M;
};
inline fuke_nuke_struct fuke_nuke_all(){
  /// calculates all fuke-nuke order parameters at once
  /// and the magnetisation
  int full_sum_abs_x = 0;
  int full_sum_sq_x  = 0;
  int full_sum_abs_y = 0;
  int full_sum_sq_y  = 0;
  int full_sum_abs_z = 0;
  int full_sum_sq_z  = 0;
  int M = 0;
  for ( int a = 0; a < L; a++ ){
    int sum_of_products_x = 0;
    int sum_of_products_y = 0;
    int sum_of_products_z = 0;
    for ( int b = 0; b < L; b++ ) {
      for ( int c = 0; c < L; c++ ) {
        sum_of_products_x += spin(a,b,c) * spin(a+1,b,c);
        sum_of_products_y += spin(b,a,c) * spin(b,a+1,c);
        sum_of_products_z += spin(b,c,a) * spin(b,c,a+1);
        M += spin(a,b,c);
      }
    }
    full_sum_abs_x += std::abs( sum_of_products_x );
    full_sum_sq_x  += sum_of_products_x*sum_of_products_x;
    full_sum_abs_y += std::abs( sum_of_products_y );
    full_sum_sq_y  += sum_of_products_y*sum_of_products_y;
    full_sum_abs_z += std::abs( sum_of_products_z );
    full_sum_sq_z  += sum_of_products_z*sum_of_products_z;
  }
  fuke_nuke_struct ret;
  ret.abs_x = full_sum_abs_x;
  ret.sq_x  = full_sum_sq_x;
  ret.abs_y = full_sum_abs_y;
  ret.sq_y  = full_sum_sq_y;
  ret.abs_z = full_sum_abs_z;
  ret.sq_z  = full_sum_sq_z;
  ret.M = M;
  return ret;
}
void initialize_state(){
  for( int x = 0; x < L; x++ )
    for( int y = 0; y < L; y++ )
      for( int z = 0; z < L; z++ ) {
        int index = cartesian_to_index( x, y, z);
        if ( x < L/2 ) {
          state[index] = +1;
        }
        else {
          state[index] = -1;
        }
        /*
        if ( rng.rand() < 0.5 )
          state[index] = 1;
        else
          state[index] = -1;
        */
  }
  // check validity of the state!
}
int main(int argc, char** argv) {

  if ( argc > 2 ){
    std::cout << "Error: parsing commandline Usage: [binary] {optional: checkpoint_nr}\n";
    exit(1);
  }

  std::signal(SIGINT, signal_callback_handler);
  std::signal(SIGTERM, signal_callback_handler);

  if (argc == 2){
    logfile.open("simulation.dat", std::ios::app);
    from_checkpoint(std::stoul(argv[1]));
  }
  else {
    logfile.open("simulation.dat");
    logfile << get_header();
    // test whether calculation of index and cartesian coordinates are 
    // at least consistently wrong ;-)
    assert(test_indices_calculations());
    // initialization
    std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
    // clear containers for muca iterations 
    std::fill(weights.begin(), weights.end(), logval(1.0));
    std::fill(transition.begin(), transition.end(), logval(0.0));
    std::fill(statistics.begin(), statistics.end(), logval(0.0));

# ifdef STORE_CARTESIAN
    for ( int index = 0; index < V; index++ ){
      spin_x[index] = index % L;
      spin_y[index] = (index/L) % L;
      spin_z[index] = (index/(L*L));
    }
# endif

    initialize_state();
    total_energy = energy();

    assert( total_energy >= minimal_energy and total_energy <= maximal_energy );
    if( total_energy < minimal_energy or total_energy > maximal_energy ){
      std::cerr << "Error: initialized state is not in the desired energy interval!\n";
      return 1;
    }


# ifdef STORE_CARTESIAN
    logfile << "# stored cartesian coordinates\n";
# endif

    logfile << "# Total energy after state initialization: ";
    logfile << total_energy << std::endl;
    
    time_start = std::time(NULL);
    // equilibrate the system
    for ( int m = 0; m < thermalization; m++ ) {
      sweep();
    }
    
    // calculate weights
    time_before_weights = std::time(NULL);

    logfile << "# weights..." << std::endl;
    calculate_weights();
  
    // write weights to file
    std::ofstream weights_file; weights_file.open("weights.dat");
    weights_file << get_header();
    weights_file << get_statistics();
    weights_file << weights;
    weights_file.close();

    logfile << "# total number of updates stored in weights = " << ln_total_updates_for_weights.toDouble() << "\n";
    logfile << "# total number of updates conducted to construct weights = " << ln_total_updates.toDouble() << "\n";
    // measurements
    logfile << "# measurements..." << std::endl;

    // reset counters messed up by calculation of the weights
    ln_total_updates = logval(0.0);
    ln_accepted_updates = logval(0.0);
    std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
    std::fill(lnhisto_all_sweeps.begin(), lnhisto_all_sweeps.end(), logval(0.0));
    checkpoint();
  }

  
  // prepare measurements file
  meas_file.open("measurements.dat");
  meas_file << get_header();
  meas_file << "# energy, magnetisation, 6 * fuke_nuke order parameters\n";
  meas_file << "#";
  meas_file << std::setw(21) << "E";
  meas_file << std::setw(22) << "M";
  meas_file << std::setw(22) << "M^\\mathrm{abs}_x";
  meas_file << std::setw(22) << "M^\\mathrm{abs}_y";
  meas_file << std::setw(22) << "M^\\mathrm{abs}_z";
  meas_file << std::setw(22) << "M^\\mathrm{sq}_x";
  meas_file << std::setw(22) << "M^\\mathrm{sq}_y";
  meas_file << std::setw(22) << "M^\\mathrm{sq}_z";
  meas_file << "\n";
  
  // to store the energy for every sweep
  energy_meas_file.open("measurements.full.dat");
  energy_meas_file << "# full timeseries, measured every sweep\n";
  energy_meas_file << get_header();
  energy_meas_file << "# energy\n";
  energy_meas_file << "# ";
  energy_meas_file << "E\n";

  // for statistics take current time
  time_before_measurement = std::time(NULL);
  while( current_m < measurements ){
    if ( current_checkpoint_counter == checkpoint_every ){
      current_checkpoint_counter = 0;
      checkpoint();
    }
    while (current_m_skip < measure_every){
      defer_signal++;
      sweep();
      energy_meas_file << total_energy << "\n";
      lnhisto_all_sweeps[total_energy]++;
      current_m_skip++;
      defer_signal--;
      if ( signal_pending != 0 ){
        raise( signal_pending );
      }
    }
    lnhisto[total_energy]++;
    fuke_nuke_struct fuke = fuke_nuke_all();
    meas_file << std::setw(22) << total_energy;
    meas_file << std::setw(22) << fuke.M;
    meas_file << std::setw(22) << fuke.abs_x;
    meas_file << std::setw(22) << fuke.abs_y;
    meas_file << std::setw(22) << fuke.abs_z;
    meas_file << std::setw(22) << fuke.sq_x;
    meas_file << std::setw(22) << fuke.sq_y;
    meas_file << std::setw(22) << fuke.sq_z;
    meas_file << "\n";
//    for ( int i=0; i < layer_updates; i++){
//      int layer = floor( rng.rand() * L );
//      int dir = floor( rng.rand() * 3 );
//      update_layer(layer, dir);
//      fuke_nuke_struct fuke = fuke_nuke_all();
//    }
    current_checkpoint_counter++;
    current_m++;
    current_m_skip = 0;
  }
  time_after_measurement = std::time(NULL);
  meas_file.close();
  energy_meas_file.close();
  
  // write histograms
  hists2file();

  logfile << get_statistics();
  logfile.close();

  return 0;
} 

void checkpoint(){
  if ( current_m != checkpoint_m ){
    checkpoints++;
    logfile << "# Creating checkpoint nr = " << std::setw(10) << checkpoints
      << " at measurement = "  << std::setw(10)
      << current_m 
      << " ..." ;
    try {
      std::ostringstream fname;
      fname << "checkpoint_";
      fname << checkpoints;
      fname << ".dump";
      const std::string filename(fname.str());

      std::ofstream dump(filename.c_str(), std::ios::binary);
      dlib::serialize(checkpoints, dump);
      dlib::serialize(current_checkpoint_counter, dump);
      dlib::serialize(current_m, dump);
      dlib::serialize(current_m_skip, dump);
      dlib::serialize(nr_updates_layer, dump);
      dlib::serialize(nr_accepted_updates_layer, dump);
      dlib::serialize(total_energy, dump);
      dlib::serialize(state, dump);
    //  dlib::serialize(low_bound, dump);
    //  dlib::serialize(high_bound, dump);
      serialize(ln_total_updates, dump);
      serialize(ln_total_updates_for_weights, dump);
      serialize(ln_accepted_updates, dump);
      serialize(rng, dump);
      serialize(lnhisto, dump);
      serialize(lnhisto_all_sweeps, dump);
      serialize(weights, dump);
      serialize(statistics, dump);
      serialize(transition, dump);
      dlib::serialize(touched_min_energy, dump);
      dlib::serialize(touched_max_energy, dump);
# ifdef STORE_CARTESIAN
      dlib::serialize(spin_x, dump);
      dlib::serialize(spin_y, dump);
      dlib::serialize(spin_z, dump);
# endif
      dump.close();
    }
    catch (std::exception& e){
      logfile << " FAILED!" << std::endl;
      logfile.close();
      std::cerr << "Error writing checkpoint number " << checkpoints << ". ";
      std::cerr << e.what() << std::endl;
    }
    logfile << " success!" << std::endl;
  }
}
void from_checkpoint(unsigned nr){
  std::ostringstream fname;
  fname << "checkpoint_";
  fname << nr;
  fname << ".dump";
  const std::string filename(fname.str());

  logfile << "# reading dump called " << filename << " ...";

  try {
    std::ifstream dump(filename.c_str(), std::ios::binary);
    dlib::deserialize(checkpoints, dump);
    dlib::deserialize(current_checkpoint_counter, dump);
    dlib::deserialize(current_m, dump);
    checkpoint_m = current_m;
    dlib::deserialize(current_m_skip, dump);
    dlib::deserialize(nr_updates_layer, dump);
    dlib::deserialize(nr_accepted_updates_layer, dump);
    dlib::deserialize(total_energy, dump);
    dlib::deserialize(state, dump);
  //  dlib::deserialize(low_bound, dump);
  //  dlib::deserialize(high_bound, dump);
    deserialize(ln_total_updates, dump);
    deserialize(ln_total_updates_for_weights, dump);
    deserialize(ln_accepted_updates, dump);
    deserialize(rng, dump);
    deserialize(lnhisto, dump);
    deserialize(lnhisto_all_sweeps, dump);
    deserialize(weights, dump);
    deserialize(statistics, dump);
    deserialize(transition, dump);
    dlib::deserialize(touched_min_energy, dump);
    dlib::deserialize(touched_max_energy, dump);
# ifdef STORE_CARTESIAN
    dlib::deserialize(spin_x, dump);
    dlib::deserialize(spin_y, dump);
    dlib::deserialize(spin_z, dump);
# endif
    dump.close();
  }
  catch (std::exception& e){
    logfile << " FAILED!" << std::endl;
    std::cerr << "Error reading checkpoint number " << checkpoints << ". ";
    std::cerr << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }
  logfile << "# ... success!" << std::endl;
}


void
signal_callback_handler(int signum)
{
  if ( defer_signal ){
    signal_pending = signum;
  }
  else {
    logfile << "# Caught signal " << signum << "\n";
    checkpoint();
    if (logfile.is_open()){ 
      logfile << std::flush;
      logfile.close();
    }
    // measurement file
    if (meas_file.is_open()){
      meas_file << std::flush;
      meas_file.close();
    }
    // measurment of the full energy time series
    if (energy_meas_file.is_open()){
      energy_meas_file << std::flush;
      energy_meas_file.close();
    }
    hists2file();

    // Terminate program
    exit(signum);
  }
}

void hists2file(){
  std::ofstream hist_file; hist_file.open("energy_histogram.simu.dat");
  hist_file << get_header();
  hist_file << get_statistics();
  hist_file << lnhisto;
  hist_file.close();

  hist_file.open("energy_histogram.simu.full.dat");
  hist_file << "# full timeseries, measured every sweep\n";
  hist_file << get_header();
  hist_file << get_statistics();
  hist_file << lnhisto_all_sweeps;
  hist_file.close();

  // canonical histogram
  discrete_function<logval> lnhisto_cano(minimal_energy - 1.0, maximal_energy + 1.0, binsize);
  discrete_function<logval> lnhisto_cano_all_sweeps(minimal_energy - 1.0, maximal_energy + 1.0, binsize);
  logval norm (0.0);
  logval norm_all_sweeps (0.0);
  for( size_t i = 0; i < lnhisto_cano.size(); i++ ){
    lnhisto_cano[i] = lnhisto[i]/weights[i]; 
    lnhisto_cano_all_sweeps[i] = lnhisto_all_sweeps[i]/weights[i]; 
    norm += lnhisto_cano[i];
    norm_all_sweeps += lnhisto_cano_all_sweeps[i];
  }
  for( size_t i = 0; i < lnhisto_cano.size(); i++ ){
    lnhisto_cano[i] = lnhisto_cano[i]/norm; 
    lnhisto_cano_all_sweeps[i] = lnhisto_cano_all_sweeps[i]/norm_all_sweeps;
  }

  hist_file.open("energy_histogram.cano.dat");
  hist_file << get_header();
  hist_file << get_statistics();
  hist_file << lnhisto_cano;
  hist_file.close();

  hist_file.open("energy_histogram.cano.full.dat");
  hist_file << "# full timeseries, measured every sweep\n";
  hist_file << get_header();
  hist_file << get_statistics();
  hist_file << lnhisto_cano_all_sweeps;
  hist_file.close();

}
