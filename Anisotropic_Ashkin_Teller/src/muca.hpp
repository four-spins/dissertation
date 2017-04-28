/**
 * @brief implementation of the Metropolis algorithm
 *
 * @author Marco Mueller
 * @date 18.10.2011
 *
 */

# ifndef __MUCA_HPP__
# define __MUCA_HPP__

// include constant stuff
# include "const.hpp"
// simulation parameters
# include "simulation.hpp"
// neighbour tables
# include "lattice.hpp"
// energy calculations
# include "hamiltonian.hpp"

// random number generator
# include "rngclass.dsfmt.hpp" // marsaglia rng

// arithmetics on logarithmic values
# include "logval.hpp"

class mc_algorithm
{
    public:
        mc_algorithm(simulation const &parameters_, lattice const &geometry_);
        double get_acceptance_ratio();

        inline void sweep(std::ofstream&);
        inline void update(uint32);

        const char * get_rng_type(){  return rng.get_type(); }
        void energy_histogram_to_file(std::ofstream &);
        void weights_to_file(std::ofstream &);
        void muca_energy_histogram_to_file(std::ofstream &);
        bool converged_flatness();
        bool converged_range();
        int32 get_number_of_iterations();

    private:
        simulation parameters;
        lattice geometry;
        hamiltonian energy;

        // random number generator
        rngclass rng;

        // delivered by parameters - instance of simulation class, 
        // but for speed also stored here
        // interaction weight

        // inverse temperature
        double beta;
        // accepted steps
        logval accepted_steps;
        logval total_steps;
        
        // number of spins in the system
        uint32 number_of_spins;

        // micro states are stored in a 1d vector
        int32_1d sigma;
        int32_1d tau;
        void initialize_states();
        
        double total_energy;

       
        // MuCa-Simulation stuff follows: 
        inline void calculate_muca_weights();
        inline void muca_sweep();
        inline void muca_update(uint32 site, bool is_sigma);
        double minimal_energy_abs;
        uint32 number_of_bins;
        inline int32 delta_energy_index(double dE);
        inline int32 energy_index(double E);
        uint32 dim;
        int32 maximal_energy_difference;
        int32 max_tunnel_events;

        // minimal energy of interest (tests for flatness ...)
        double E_min;
        // maximal energy of interest
        double E_max;

        std::vector<logval> lnhisto;
        std::vector<logval> statistics;
        std::vector<logval> weights;
        std::vector<logval> transition;
        std::vector<logval> energy_histogram;

        // for printing out timeseries of MuCa-Simulation
        // parameters for tests on convergency ...
        std::ofstream weights_file;
        std::ofstream transition_file;
        std::ofstream histo_file;
        
        // stores whether muca converged according to
        // flatness and full energy range
        bool muca_flatness;
        bool muca_full_range;
        int32 number_of_iterations;
        bool weights_over_time;

};

mc_algorithm::mc_algorithm(simulation const &parameters_, lattice const &geometry_)
    : parameters(parameters_), geometry(geometry_), energy(parameters_, geometry_),
      rng(parameters.get_seed()), beta(parameters.get_beta()), 
      accepted_steps(0.0), total_steps(0.0), number_of_spins(parameters.get_number_of_spins()),
      max_tunnel_events(parameters.get_max_tunnel_events()),
      E_min(parameters.get_E_min()), E_max(parameters.get_E_max()), 
      muca_flatness(false), muca_full_range(false),
      number_of_iterations(0), weights_over_time(parameters.get_weights_over_time())
{
    const int32 max_init = 100;
    int32 init = 0;
    do
    {
        if (init >= max_init)
        {
            std::cerr << "Error: Bad energy interval: could not create initial " 
                      << "state in the desired range using " << init
                      << " attempts.\n";
            exit(1);
        }
        initialize_states();
        total_energy = energy.calculate_total_energy(sigma, tau);
        init ++;
    }
    while (total_energy < E_min || total_energy > E_max);

    dim = parameters.get_dimension();
    // Every spin is member of 2d(d-1) plaquettes, where d is the dimension
    // of the lattice. Since every plaquette consists of 4 spins,
    // the total number of different plaquettes is one quarter of
    // the total number of spins times 2d(d-1), thus the total number of 
    // plaquettes x reads:
    //      x = d*(d - 1)/2*number_of_spins;
    // Every plaquette can contribute either an energy of +1 or -1, therefore
    // the energy range is [-x, x].
    // Now a single spin flip changes the sign of all 2d(d-1) plaquettes
    // touching the spin. So the step size of the energy is dE = 2*2d(d-1),
    // thus leading to the total number of bins: 2x/dE:
    minimal_energy_abs = 1.0/4.0*dim*(dim-1)*number_of_spins;
    
    // if interface conditions are set, the total number of spins contributing to
    // the energy is +2 for each dimension
    if (parameters.get_boundary_conditions() == INTERFACE)
    {
        uint32_1d dimensions;
        parameters.get_dimensions(dimensions);
        uint32 volume = 1;
        for(size_t i = 0; i < dimensions.size(); i++)
        {
            volume *= (dimensions[i] + 2);
        }
        minimal_energy_abs = 1.0/4.0*dim*(dim-1)*volume;
    }

    number_of_bins = ceil(minimal_energy_abs) + 1; 
        // parameters.get_number_of_bins(); //100; //ceil(minimal_energy_abs) + 1;
    
    maximal_energy_difference = 2*2*static_cast<int>(dim*(dim - 1));
     
    // histogram for every recursion step
    // using 1 as initialization introduces a statistical error,
    // but this should be neglectable for the calculation of weights
    lnhisto.resize(number_of_bins, logval(0.0));
    weights.resize(number_of_bins, logval(1.0));
    
    energy_histogram.resize(number_of_bins, logval(0.0));
    // the transition vector has one element less than the histogram
    transition.resize(number_of_bins - 1, logval(1.0));
    // so has the vector collecting statistics
    statistics.resize(number_of_bins - 1, logval(0.0));

    calculate_muca_weights();

}

int32 mc_algorithm::get_number_of_iterations()
{ return number_of_iterations; }

bool mc_algorithm::converged_flatness()
{ return muca_flatness; }

bool mc_algorithm::converged_range()
{ return muca_full_range; }

inline int32 mc_algorithm::delta_energy_index(double dE)
{
    return energy_index(total_energy+dE);
}

inline int32 mc_algorithm::energy_index(double E)
{
    //int32 index = static_cast<int32>(E + minimal_energy_abs)/2;
    //double binwidth = (2.0*minimal_energy_abs)/(number_of_bins-1);
    //std::cout << binwidth << "\t";
    int32 index = static_cast<int32>(E + minimal_energy_abs)/2;
    
    // remove prohibited states
//    if (index > 1) index--;
/*
    assert(E >= -1 * minimal_energy_abs);
    assert(index >= 0);
    if (index >= static_cast<int32>(number_of_bins))
    {
        std::cout << "E = " << E << ", i=" << index;// << ", b = " << binwidth;
        std::cout << ", E + Emin: " << E + minimal_energy_abs << "\n";
    }
    assert(index < static_cast<int32>(number_of_bins));
*/
    return index;
}

void mc_algorithm::calculate_muca_weights()
{
    std::cout << "# calculation of multicanonical weights started.\n";
    
    logval exponent; // kappa in the original source
    logval fraction;

    weights_file.open("weights.muca.time.dat");
    histo_file.open("histo.muca.time.dat");
    transition_file.open("transition.muca.time.dat");
    if (! weights_over_time)
    {
        weights_file << "# was disabled in configuration file.\n";
        histo_file << "# was disabled in configuration file.\n";
        transition_file << "# was disabled in configuration file.\n";
    }


    std::vector<logval> weights_;
    std::vector<double> energies_;
    if(parameters.weights_from_file())
    {
        // retrieve weights from file and check for the right size
        parameters.get_weights(weights_, energies_);
        std::fill(weights.begin(), weights.end(), logval(1.0));
        for (size_t i = 0; i < energies_.size(); i++)
        {
            weights[energy_index(energies_[i])] = weights_[i];
        }
        for(size_t i=0; i < transition.size(); i++)
            transition[i] =  weights[i+1]/weights[i];
        weights_file    << join(weights)  << std::endl;
    }

// TODO: Fine tuning of the following:
    // allow a bit more than a random walk through energy space
    // every sweep touches all spins, a random walk needs (at least) sweeps**2
    // old version:
    // uint32 sweeps_per_muca_step = number_of_spins*number_of_spins;

    // 19.08.2011: Allow for a random walk in the histogram instead of 
    // state space (actually this expression should be the very same order
    // as the version above):
    // uint32 sweeps_per_muca_step = 10*number_of_bins*number_of_bins;

    // maximum number of recursions to be conducted
    int32 muca_steps = parameters.get_muca_steps();
    // are there entries in the energy range [- E_min_muca, E_max_muca]
    bool energy_range = false;
    bool flatten_histogram = false;
    int tunnel_events = 0;
    // whether the histogram is "flat" (ie. H[i] < exp(3)*<H> for all i)
    bool flat_histogram = false;
    bool touched_minE = false;
    bool touched_maxE = false;
    // tests for a larger energy range:
    double last_energy_min = 1.0e100;
    double last_energy_max = -1.0e100;

    int32 E_min_index = energy_index(E_min);
    int32 E_max_index = energy_index(E_max);
    
    // number of sweeps in this iteration
    uint32 iteration_sweeps = 0;
    while ((muca_steps > 0)  /*&& 
           (!energy_range || !flat_histogram)*/
           && (tunnel_events < max_tunnel_events))
    {
        double cur_energy_min = 1.0e100;
        double cur_energy_max = -1.0e100;

        std::cout << "# iteration " 
                  << 1 + parameters.get_muca_steps() - muca_steps
                  << " started." << std::endl;

        // clear histogram 
        std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
        // now conduct a number of sweeps:
        if ( parameters.weights_from_file() )
        {
            iteration_sweeps = static_cast<int>(weights.size()*weights.size());
        }
       
        iteration_sweeps = 1000;
        if (flatten_histogram)
        {
            iteration_sweeps = std::pow(
                energy_index(last_energy_max) - energy_index(last_energy_min), 2);
        }
        //iteration_sweeps = number_of_bins*number_of_bins;
        muca_steps --;

        if ( !parameters.weights_from_file() 
            || !(muca_steps == parameters.get_muca_steps()))
        {
            for(uint32 step = 0; step < iteration_sweeps; step++)
            {
                for (uint32 i = 0; i < 2*number_of_spins; i++)
                {
                    // choose a site randomly (to avoid additional correlations)
                    double r = rng.rand(); // in (0,1)
                    uint32 site = static_cast<uint32>(floor(r * 2*number_of_spins)); 

                    // Update either sigma or tau
                    if (site >= number_of_spins) // flip a sigma spin
                    {
                        muca_update(site - number_of_spins, true);
                    }
                    else
                    {
                        muca_update(site, false);
                    }

                    double curE = total_energy; // current energy to be added
                    cur_energy_min = std::min(curE, cur_energy_min);
                    cur_energy_max = std::max(curE, cur_energy_max);
                    lnhisto[delta_energy_index(0)] ++;
                }
            }
        }
        double number_of_measurements = iteration_sweeps*number_of_spins;

        // assume the requirements to abort are met
        energy_range = true;
        flat_histogram = true;

        // assume no new energy was found
        last_energy_min = std::min(cur_energy_min, last_energy_min);
        last_energy_max = std::max(cur_energy_max, last_energy_max);

        if ( last_energy_min > E_min || last_energy_max < E_max)
        {
            energy_range = false;
        }
        
        std::cout << energy_index(cur_energy_min) << " : " << energy_index(E_min) << "\n";
        std::cout << energy_index(cur_energy_max) << " : " << energy_index(E_max) << "\n";
        if ( energy_index(cur_energy_min) <= energy_index(E_min) )
            touched_minE = true;
        
        if ( energy_index(cur_energy_max) >= energy_index(E_max) )
            touched_maxE = true;

        if (touched_minE && touched_maxE)
        {
            tunnel_events ++;
            touched_minE = false;
            touched_maxE = false;
            //if (tunnel_events > max_tunnel_events-2)
                //flatten_histogram = true;
        }

        std::cout << "# tunnel events: " << tunnel_events << "\n";
        // accumulate statistics
        for (uint32 i = 0; i < number_of_bins - 1; i++)
        {
            if(lnhisto[i] == logval(0.0) || lnhisto[i+1] == logval(0.0))
            {
                exponent = logval(1.0);
                fraction = logval(1.0);
            }
            else
            {
                // logval histosum = lnhisto[i] + lnhisto[i+1];
                // double tester = lnhisto[i].value() + lnhisto[i+1].value() 
                //    - histosum.value();
                fraction = lnhisto[i]*lnhisto[i+1]/(lnhisto[i]+lnhisto[i+1]);
                // assert(abs(fraction.value() - tester) < 0.0000001);
                statistics[i] += fraction;
                exponent = fraction/statistics[i];
                fraction = lnhisto[i] / lnhisto[i+1];
                double flat = 
                    std::abs(
                        std::log(
                            static_cast<double>(
                                energy_index(cur_energy_max) - energy_index(cur_energy_min))
                        ) + lnhisto[i].value() - 
                        std::log(number_of_measurements)
                    );


                if (static_cast<int32>(i) >= E_min_index && static_cast<int32>(i) <= E_max_index && flat > 3)
                {
                    //std::cout << "flatness: " << flat << "\n";
                    flat_histogram = false;
                }
            }

            assert(fraction == fraction); // check for NaNs
            transition[i] = transition[i]*pow(fraction.toDouble(), exponent.toDouble());
            weights[i+1] = transition[i]*weights[i];
        }
        
        number_of_iterations = parameters.get_muca_steps() - muca_steps; 

        // std::cout << "# Emin " << cur_energy_min << ", Emax " << cur_energy_max << "\n";
        // std::cout << "# Emin* " << E_min << ", Emax* " << E_max << "\n";
        // std::cout << "# ind(Emin) " << energy_index(cur_energy_min) << ", ind(Emax) " << energy_index(cur_energy_max) << "\n";
        // std::cout << "# ind(Emin*) " << energy_index(E_min) << ", ind(Emax*) " << energy_index(E_max) << "\n";
        // std::cout << "# "; 
        if (flat_histogram) 
        { 
           //  std::cout << "flat "; 
            muca_flatness = true; 
        } 
        if (energy_range) 
        { 
           //  std::cout << "filled "; 
            muca_full_range = true;     
        } 
        // std::cout << std::endl; 

        if (weights_over_time)
        {
            weights_file    << join(weights)        << std::endl;
            transition_file << join(transition)     << std::endl;
            histo_file      << join(lnhisto)        << std::endl;
        }
    }
    weights_file.close();
    transition_file.close();
    histo_file.close();
    accepted_steps = logval(0.0);
    total_steps = logval(0.0);
}

double mc_algorithm::get_acceptance_ratio()
{
    logval ratio = accepted_steps/total_steps;
    return ratio.toDouble();
}

inline void mc_algorithm::sweep(std::ofstream &out)
{
    muca_sweep();
    out << total_energy << "\n";
}

inline void mc_algorithm::muca_sweep()
{
    for (uint32 i = 0; i < 2*number_of_spins; i++)
    {
        // choose a site randomly (to avoid additional correlations)
        double r = rng.rand(); // in (0,1)
        uint32 site = static_cast<uint32>(floor(r * 2*number_of_spins));

        // Update either sigma or tau
        if (site >= number_of_spins) // flip a sigma spin
        {
            muca_update(site - number_of_spins, true);
        }
        else
        {
            muca_update(site, false);
        }
    }
    energy_histogram[energy_index(total_energy)]++;
}


inline void mc_algorithm::muca_update(uint32 site, bool is_sigma)
{
    double dEloc = energy.get_local_energy_difference(site, sigma, tau, is_sigma);

#   ifdef DEBUG
        // calculate energy difference slowly for comparison of correct implementation
        // of local spin flip:
        // flip
        if (is_sigma) sigma[site] *= -1;
        else tau[site] *= -1;
        double new_energy = energy.calculate_total_energy(sigma, tau);
        // reflip (original state)
        if (is_sigma) sigma[site] *= -1;
        else tau[site] *= -1;
        assert(std::abs((new_energy - total_energy) - dEloc) < 1e-6);
#   endif // DEBUG

    // check whether new state falls into energy interval of interest
    // accept with multicanonical probability if so, else reject new state
    if (total_energy + dEloc - E_min > -1e-6 && total_energy + dEloc - E_max < 1e-6)
    {
        double r = rng.rand();
        logval frac = weights[delta_energy_index(dEloc)] / weights[delta_energy_index(0)];
        double p_muca = exp(-beta*dEloc)*frac.toDouble();
        // acceptance probability = min(p_muca, 1.0)
        // but r is in range (0,1] therefore the calculation can be omitted
        if (r < p_muca)
        {
            // accept
            if (is_sigma) sigma[site] *= -1;
            else tau[site] *= -1;

            total_energy += dEloc;
            // spin was already flipped --> (1)
            accepted_steps++;
        }
    }
    total_steps ++;
}

void mc_algorithm::energy_histogram_to_file(std::ofstream &out)
{
    out << "#" << std::setw(19) << "E" << std::setw(30) << "log(H(E))\n";
    
    int32 last_ind = static_cast<int32>(energy_index(E_min)) - 1;
    for (int32 E = static_cast<int32>(E_min); E <= static_cast<int32>(E_max); E+=2)
    {  
        int32 ind = energy_index(E); 
        // skipping prohibited energies
        if (ind != last_ind)
        {
            out.unsetf(std::ios::floatfield);
            out << std::setw(20) << E ;
            out.setf(std::ios::scientific, std::ios::floatfield);
            out.precision(16);
            out << std::setw(30) << energy_histogram[ind]/weights[ind];
            out << "\n";
        }
        last_ind = ind;
    }
}

void mc_algorithm::weights_to_file(std::ofstream &out)
{
    double binwidth = (2.0*minimal_energy_abs + 1)/(number_of_bins);
    out << "# binned weights:\n# number of bins: " << number_of_bins
        << "\n# binwidth: " << binwidth << "\n";
    out << "#" << std::setw(19) << "E" << std::setw(30) << "log(w(E))\n";

    int32 last_ind = static_cast<int32>(energy_index(E_min)) - 1;
    for (int32 E = static_cast<int32>(E_min); E <= static_cast<int32>(E_max); E+=2)
    {  
        int32 ind = energy_index(E); 
        // skipping prohibited energies
        if (ind != last_ind)
        {
            out.unsetf(std::ios::floatfield);
            out << std::setw(20) << E ;
            out.setf(std::ios::scientific, std::ios::floatfield);
            out.precision(16);
            out << std::setw(30) << weights[ind];
            out << "\n";
        }
        last_ind = ind;
    }
}

void mc_algorithm::muca_energy_histogram_to_file(std::ofstream &out)
{
    out << "#" << std::setw(19) << "E" << std::setw(30) << "log(H_muca(E))\n";
    int32 last_ind = static_cast<int32>(energy_index(E_min)) - 1;
    for (int32 E = static_cast<int32>(E_min); E <= static_cast<int32>(E_max); E+=2)
    {  
        int32 ind = energy_index(E); 
        // skipping prohibited energies
        if (ind != last_ind)
        {
            out.unsetf(std::ios::floatfield);
            out << std::setw(20) << E ;
            out.setf(std::ios::scientific, std::ios::floatfield);
            out.precision(16);
            out << std::setw(30) << energy_histogram[ind];
            out << "\n";
        }
        last_ind = ind;
    }
}


void mc_algorithm::initialize_states()
{
    if (!parameters.get_initial_state(sigma))
    {
        for (uint32 n = 0; n < number_of_spins; n++)
        {
/*
            int32 spin = (rng.rand() < 0.5) ? -1 : 1;
            sigma.push_back(spin);

            spin = (rng.rand() < 0.5) ? -1 : 1;
            tau.push_back(spin);
*/
            sigma.push_back(1);
            tau.push_back(1);
        }
        if(parameters.get_boundary_conditions() == INTERFACE)
        {
            // set virtual boundary spins
            // (see lattice.hpp for details on the implementation of 
            // interface boundary conditions)
            sigma.push_back(1);
            sigma.push_back(-1);
            tau.push_back(1);
            tau.push_back(-1);
            std::cerr << "No sense of having interface boundary conditions in dual representation.\n";
            exit(1);
        }
        if (parameters.get_boundary_conditions() == VANISH)
        {
            sigma.push_back(0);
            sigma.push_back(0);
            tau.push_back(0);
            tau.push_back(0);
        }
    } // if (!parameters.get_initial_state(..))
    else
    {
        std::cerr << "Cannot use initial state -- not implemented for dual representation.\n";
        exit(1);
    }
}

# endif
