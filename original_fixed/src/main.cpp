/******************************************************************************
 *
 * @file: main.cpp
 *
 * @date: 02/27/2013 11:35:09 AM (CET)
 *
 * @author: Marco Müller <muelma@gmail.com>
 *
 ******************************************************************************/
# include <iostream>
# include <iomanip>
# include <vector>
# include <cassert>
# include <cmath>
# include <fstream>
# include <ctime>
# include <sstream>
// random number generator
# include "rngclass.mar.hpp"
// histograms and such
# include "discrete_function.hpp"
// type for logarithmic calculations
# include "logval.hpp"

// stores all the necessary parameters
# include "config.hpp" 
// helpers
logval ln_total_updates(0.0);     // total number of updates
logval ln_total_updates_for_weights(0.0);

logval ln_accepted_updates(0.0);  // number of accepted updates
double total_energy;              // current total energy of the system

rngclass rng(seed);

std::vector<int> state(V, +1);  // stores the current configuration
// stores the histogram and weights
discrete_function<logval> lnhisto(minimal_energy - 1.0, maximal_energy + 1.0, binsize);
discrete_function<logval> weights(minimal_energy - 1.0, maximal_energy + 1.0, binsize);
discrete_function<logval> statistics(minimal_energy - 1.0, maximal_energy + 1.0, binsize);
discrete_function<logval> transition(minimal_energy - 1.0, maximal_energy + 1.0, binsize);

double tEmin = -1.0; // step when the minimal energy was first touched
double tmax;

# ifdef STORE_CARTESIAN
std::vector<int> spin_x(V, -1);  // stores the x-coordinate of the spin
std::vector<int> spin_y(V, -1);  // stores the x-coordinate of the spin
std::vector<int> spin_z(V, -1);  // stores the x-coordinate of the spin
# endif
time_t time_start, time_before_weights;
time_t time_before_measurement, time_after_measurement;
std::ofstream logfile; 

std::string get_header(){
  std::stringstream header;
  header << "# recursive Multicanonical simulation\n";
  header << "# Gonihedric Ising Model, plaquette only\n";
  time_t rawtime; time( &rawtime );
  header << "# " << ctime( &rawtime );
  header << "#\n";
  header << "# random number generator = " << rng.get_type() << "\n";
  header << "# random number seed = " << rng.get_seed() << "\n";
  header << "# boundary conditions = fixed+\n";
  header << "# beta = 0\n";
  header << "# grid dimensions = " << L << ", " << L << ", " << L << "\n";
  header << "# volume = " << L*L*L << "\n";
  header << "# E_min = " << minimal_energy << "\n";
  header << "# E_max = " << maximal_energy << "\n";
  header << "# number of measurements = " << measurements << "\n";
  header << "# measured every [this number] sweeps = " << measure_every << "\n";
  header << "# thermalization sweeps = " << thermalization << "\n";
//  header << "# wang landau tunnel events per iteration = " << tunnel_events << "\n";
  return header.str();
}
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
  return statistics.str();
}
inline int cartesian_to_index(int x, int y, int z){
  assert( x >= 0 && x < L );
  assert( y >= 0 && y < L );
  assert( z >= 0 && z < L );
  return x + L*y + L*L*z;
}
inline int index_to_x(int index){
# ifdef STORE_CARTESIAN
  if (spin_x[index] < 0)
    spin_x[index] = index % L;
  return spin_x[index];
# else
  return index % L;
# endif
}
inline int index_to_y(int index){
# ifdef STORE_CARTESIAN
  if (spin_y[index] < 0)
    spin_y[index] = (index/L) % L;
  return spin_y[index];
# else
  return (index/L) % L;
# endif
}
inline int index_to_z(int index){
# ifdef STORE_CARTESIAN
  if (spin_z[index] < 0)
    spin_z[index] = (index/(L*L));
  return spin_z[index];
# else
  return index/(L*L);
# endif
}
inline int spin(int x, int y, int z){
  // do not change order of the next two if's 
  // vanishing spins:
  if ( x < -1 or y < -1 or z < -1 or x > L or y > L or z > L ){
    return 0;
  }
  // fixed spins
  if ( x == -1 or y == -1 or z == -1 or x == L or y == L or z == L ){
    return +1;
    // uncomment to force plane 
    //if ( x < L/2 ) {
    //  return +1;
    //}
    //else {
    //  return -1;
    //}
  }
  // if you got here, assert that we are in the box:
  assert( x >= 0 && x < L );
  assert( y >= 0 && y < L );
  assert( z >= 0 && z < L );
  return state[cartesian_to_index(x,y,z)];
}
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
  // fixed boundary conditions mess up the clean calculation of energies:
  // total energy is calculated by walking through the lattice and checking
  // whether the plaquette would contain only fixed spins
  double sum = 0;
  for( int x = 0; x < L; x++ )
    for( int y = 0; y < L; y++ )
      for( int z = 0; z < L; z++ ) {
        double sumi = 0;
        // all plaquettes in direction "+"
        sumi += spin(x+1,y,z) * spin(x+1,y+1,z) * spin(x,y+1,z); // plaquette parallel to xy-plane
        sumi += spin(x+1,y,z) * spin(x+1,y,z+1) * spin(x,y,z+1); // plaquette parallel to xz-plane
        sumi += spin(x,y+1,z) * spin(x,y+1,z+1) * spin(x,y,z+1); // plaquette parallel to yz-plane
        //sumi *= spin(x,y,z);
       
        if ( x == 0 ) {
          // add plaquettes non-parallel to yz-plane in "-x"--direction and "+"--direction otherwise
          // parallel to xy-plane
          sumi += spin(x-1,y,z) * spin(x-1,y+1,z) * spin(x,y+1,z); 
          // parallel to xz-plane
          sumi += spin(x-1,y,z) * spin(x-1,y,z+1) * spin(x,y,z+1); 
        }

        if ( y == 0 ) {
          // add plaquettes non-parallel to xz-plane in "-y"--direction and "+"--direction otherwise
          // parallel to xy-plane
          sumi += spin(x+1,y,z) * spin(x+1,y-1,z) * spin(x,y-1,z); 
          // parallel to yz-plane
          sumi += spin(x,y-1,z) * spin(x,y-1,z+1) * spin(x,y,z+1);
        }

        if ( z == 0 ) {
          // add plaquettes non-parallel to xy-plane in "-z"--direction and "+"--direction otherwise
          // parallel to xz-plane
          sumi += spin(x+1,y,z) * spin(x+1,y,z-1) * spin(x,y,z-1); 
          // parallel to yz-plane
          //sumi += spin(x,y-1,z) * spin(x,y-1,z-1) * spin(x,y,z-1);
          sumi += spin(x,y+1,z) * spin(x,y+1,z-1) * spin(x,y,z-1);
        }

        // take care of the axes:
        if ( x == 0 and y == 0 ) {
          sumi += spin(x-1,y,z) * spin(x-1,y-1,z) * spin(x,y-1,z); 
        }

        if ( x == 0 and z == 0 ) {
          sumi += spin(x-1,y,z) * spin(x-1,y,z-1) * spin(x,y,z-1); 
        }

        if ( y == 0 and z == 0 ) {
          sumi += spin(x,y-1,z) * spin(x,y-1,z-1) * spin(x,y,z-1);
        }

        sumi *= spin(x,y,z);

        // std::cout << " (x,y,z) = ( " << x << " , " << y << " , " << z << " ): " << sumi << "\n";
        sum += sumi;
      }
  // now we forgot to count the plaquettes in "-" -- direction on the boundary!
  const double E = -0.5*sum;
  return E;
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
}

inline void calculate_weights(){
  // recursive determination of multicanonical weights
  //
  logval exponent; // kappa in the original paper of Janke
  logval fraction; // helper variable for better readability
  int tunnel_events = 0; // current number of tunnel events
  int muca_iterations = 0;
  bool touched_min_energy = false;
  bool touched_max_energy = false;
  // last registered minimal energy
  double last_energy_min = 1.0e100;
  double last_energy_max = -1.0e100;
  while ( tunnel_events < max_tunnel_events 
      && muca_iterations < max_muca_iterations ) {
    muca_iterations++;
    // clear histogram
    std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
    // reset current energy minima and maxima
    double cur_energy_min = 1.0e100;
    double cur_energy_max = -1.0e100;

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

        if ( total_energy <= minimal_energy ){
          // std::cout << "# touched minimal energy\n";
          touched_min_energy = true;
        } else if ( total_energy >= maximal_energy ){
          // std::cout << "# touched maximal energy\n";
          touched_max_energy = true;
        }
        if ( touched_min_energy && touched_max_energy ){
          tunnel_events++;
          // reset state
          touched_max_energy = false;
          touched_min_energy = false;
          std::cout << "# registered tunnel event, tunnel_events = " << tunnel_events << "\n";
          if ( tunnel_events >= max_tunnel_events ){
            break;
          }
        }
        if ( tunnel_events >= max_tunnel_events ){
          break;
        }
      } // sweep
    } // sweeps per iteration

    // break out of while before another update of weights is performed!
    if ( tunnel_events >= max_tunnel_events ){
      break;
    }

    // update of weights:
    last_energy_min = std::min(cur_energy_min, last_energy_min);
    last_energy_max = std::max(cur_energy_max, last_energy_max);
    if ( last_energy_min > minimal_energy ){
//      std::cout << "# linear interpolation\n";
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

struct fuki_nuke {
  int fuki_nuke_abs;
  int fuki_nuke_sq;
};
struct fuki_nuke_struct {
  int abs_x;
  int sq_x;
  int abs_y;
  int sq_y;
  int abs_z;
  int sq_z;
  int M;
};
inline fuki_nuke_struct fuki_nuke_all(){
  /// calculates all fuki-nuke order parameters at once
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
  fuki_nuke_struct ret;
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
  logfile.open("simulation.dat");
  // test whether calculation of index and cartesian coordinates are 
  // at least consistently wrong ;-)
  assert(test_indices_calculations());

  // initialization
  std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));

  // clear containers for muca iterations 
  std::fill(weights.begin(), weights.end(), logval(1.0));
  std::fill(transition.begin(), transition.end(), logval(0.0));
  std::fill(statistics.begin(), statistics.end(), logval(0.0));

  initialize_state();
  total_energy = energy();

  logfile << get_header();
  logfile << "# Total energy after state initialization: ";
  logfile << total_energy << std::endl;
  

  assert( total_energy >= minimal_energy and total_energy <= maximal_energy );
  if( total_energy < minimal_energy or total_energy > maximal_energy ){
    std::cerr << "Error: initialized state is not in the desired energy interval!\n";
    exit(1);
  }

  time_start = std::time(NULL);
  // equilibrate the system
  for ( int m = 0; m < thermalization; m++ ) {
    sweep();
  }
  
  // calculate weights
  time_before_weights = std::time(NULL);
  logfile << "# weights..." << std::endl;
  calculate_weights();
  logfile << "# total number of updates stored in weights = " << ln_total_updates_for_weights.toDouble() << "\n";
  logfile << "# total number of updates conducted to construct weights = " << ln_total_updates.toDouble() << "\n";

  // measurements
  logfile << "# measurements..." << std::endl;
  // reset counters messed up by calculation of the weights
  ln_total_updates = logval(0.0);
  ln_accepted_updates = logval(0.0);
  std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
  // prepare measurements file
  std::ofstream meas_file; meas_file.open("measurements.dat");
  meas_file << get_header();
  meas_file << "# energy, magnetisation, 6 * fuki_nuke order parameters\n";
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

  // for statistics take current time
  time_before_measurement = std::time(NULL);
  for ( int m = 0; m < measurements; m++ ) {
    for ( int skip = 0; skip < measure_every; skip ++){
      sweep();
    }
    lnhisto[total_energy]++;
    fuki_nuke_struct fuki = fuki_nuke_all();
    meas_file << std::setw(22) << total_energy;
    meas_file << std::setw(22) << fuki.M;
    meas_file << std::setw(22) << fuki.abs_x;
    meas_file << std::setw(22) << fuki.abs_y;
    meas_file << std::setw(22) << fuki.abs_z;
    meas_file << std::setw(22) << fuki.sq_x;
    meas_file << std::setw(22) << fuki.sq_y;
    meas_file << std::setw(22) << fuki.sq_z;
    meas_file << "\n";
  }
  time_after_measurement = std::time(NULL);
  meas_file.close();

  std::ofstream hist_file; hist_file.open("energy_histogram.simu.dat");
  hist_file << get_header();
  hist_file << get_statistics();
  hist_file << lnhisto;
  hist_file.close();

  std::ofstream weights_file; weights_file.open("weights.dat");
  weights_file << get_header();
  weights_file << get_statistics();
  weights_file << weights;
  weights_file.close();

  // canonical histogram
  discrete_function<logval> lnhisto_cano(minimal_energy - 1.0, maximal_energy + 1.0, binsize);
  logval norm (0.0);
  for( size_t i = 0; i < lnhisto_cano.size(); i++ ){
    lnhisto_cano[i] = lnhisto[i]/weights[i]; 
    norm += lnhisto_cano[i];
  }
  for( size_t i = 0; i < lnhisto_cano.size(); i++ ){
    lnhisto_cano[i] = lnhisto_cano[i]/norm; 
  }

  std::ofstream hist_can_file; hist_can_file.open("energy_histogram.cano.dat");
  hist_can_file << get_header();
  hist_can_file << get_statistics();
  hist_can_file << lnhisto_cano;
  hist_can_file.close();

  logfile << get_statistics();
  logfile.close();

  return 0;
} 

