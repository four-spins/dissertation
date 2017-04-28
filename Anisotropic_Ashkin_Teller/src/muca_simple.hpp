/**
 * @brief implementation of the Metropolis algorithm
 *
 * @author Marco Mueller
 * @date 11.05.2011
 *
 */

# ifndef __MUCA_SIMPLE_HPP__
# define __MUCA_SIMPLE_HPP__

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
        
        // number of spins in the system
        uint32 number_of_spins;

        // accepted steps
        logval accepted_steps;
        logval total_steps;

        // micro states are stored in a 1d vector
        int32_1d sigma;
        int32_1d tau;

        void initialize_states();
        
        double total_energy;

        typedef std::map<double, logval> histogram_type;
        histogram_type energy_histogram;
       
        // MuCa-Simulation stuff follows: 
        inline void calculate_muca_weights();
        inline void muca_sweep();
        inline void muca_update(uint32 site, bool is_sigma);
        inline int32 delta_energy_index(double dE);
        inline int32 energy_index(double E);

        uint32 number_of_bins;
        uint32 dim;
        // minimal energy that can be reached by the system
        double minimal_energy_abs;
        // biggest energy difference that can be conducted by an update
        int32 maximal_energy_difference;
        // minimal energy of interest (tests for flatness ...)
        double E_min;
        // maximal energy of interest
        double E_max;

        std::vector<logval> lnhisto;
        std::vector<logval> weights;

        // for printing out timeseries of MuCa-Simulation
        // parameters for tests on convergency ...
        std::ofstream weights_file;
        std::ofstream histo_file;
        
        // stores whether muca converged according to
        // flatness and full energy range
        bool muca_flatness;
        bool muca_full_range;
        int32 number_of_iterations;

};

mc_algorithm::mc_algorithm(
    simulation const &parameters_, lattice const &geometry_)
    : parameters(parameters_), geometry(geometry_), 
      energy(parameters_, geometry_), rng(parameters.get_seed()), 
      beta(parameters.get_beta()), 
      number_of_spins(parameters.get_number_of_spins()),
      accepted_steps(0.0), total_steps(0.0),
      dim(parameters.get_dimension()), E_min(parameters.get_E_min()),
      E_max(parameters.get_E_max()), muca_flatness(false),
      muca_full_range(false), number_of_iterations(0)
{

    const int32 max_init = 100;
    int32 init = 0;
    do
    {
        if (init >= max_init)
        {
            std::cerr << "Error: Bad energy interval: could not create initial "
                      << "state in the desired range using " << init 
                      << " attempts.";
            exit(1);
        }
        initialize_states();
        total_energy = energy.calculate_total_energy(sigma, tau);
        init ++;
    }
    while (total_energy < E_min || total_energy > E_max);

    
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
    minimal_energy_abs = 1.5*number_of_spins;
    
    number_of_bins = ceil(minimal_energy_abs) + 1;
     
    // histogram for every recursion step
    // using 1 as initialization introduces a statistical error,
    // but this should be neglectable for the calculation of weights
    lnhisto.resize(number_of_bins, logval(0.0));
    weights.resize(number_of_bins, logval(1.0));
    
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
    int32 index = static_cast<int32>(E + minimal_energy_abs)/2;
# ifdef DEBUG
    assert(E >= -1 * minimal_energy_abs);
    assert(index >= 0);
    assert(index < static_cast<int32>(number_of_bins));
# endif // DEBUG
    return index;
}

void mc_algorithm::calculate_muca_weights()
{
    std::cout << "# calculation of multicanonical weights started.\n";
    
    /* 
    // thermalize
        for(uint32 step = 0; step < 10000; step++)
        {
            for (uint32 i = 0; i < number_of_spins; i++)
            {
                // choose a site randomly (to avoid additional correlations)
                double r = rng.rand(); // in (0,1)
                uint32 site = static_cast<uint32>(floor(r * number_of_spins));
                // Monte Carlo Update
                muca_update(site);
            }
        }
    */
    weights_file.open("weights.muca.time.dat");
    histo_file.open("histo.muca.time.dat");

    std::vector<double> energies_;
    if(parameters.weights_from_file())
    {
        // retrieve weights from file and check for the right size
        parameters.get_weights(weights, energies_);
        assert(weights.size() == number_of_bins);           
        weights_file    << join(weights)  << std::endl;
    }

    // 19.08.2011: Allow for a random walk in the histogram instead of 
    // state space (actually this expression should be the very same order
    // as the version above):
    //  uint32 sweeps_per_muca_step = 10*number_of_bins*number_of_bins;

    // maximum number of recursions to be conducted
    int32 muca_steps = parameters.get_muca_steps();

    // are there entries in the energy range [- E_min, E_max]
    bool energy_range = false;

    // whether the histogram is "flat" (ie. H[i] < exp(3)*<H> for all i)
    bool flat_histogram = false;
    double flatfactor = 0.05;
    
    // number of sweeps in this iteration
    uint32 iteration_sweeps = 0;
    while ( (muca_steps > 0)  &&  !(flat_histogram && energy_range) )
    {

        double last_min_energy = 1e10;
        double last_max_energy = -1e10;
        std::cout << "# iteration " 
                  << 1 + parameters.get_muca_steps() - muca_steps
                  << " started. " << std::endl;

        // clear histogram 
        std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
        // now conduct a number of sweeps:
        // iteration_sweeps = pow(std::max(number_of_bins/muca_steps, (uint32)10), 2);
        iteration_sweeps = 
            std::pow(energy_index(E_max) - energy_index(E_min), 2);

        logval ln_number_of_measurements = logval(0.0);
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
                last_min_energy = std::min(total_energy, last_min_energy);
                last_max_energy = std::max(total_energy, last_max_energy);
            }
            lnhisto[delta_energy_index(0)] ++;
            ln_number_of_measurements ++;
        }
        // assume the requirements to abort are met
        energy_range = true;
        flat_histogram = true;

        if(last_min_energy > E_min || last_max_energy < E_max)
        {
            energy_range = false;
        }

        logval H_min = ln_number_of_measurements;
        logval H_max = logval(0.0);
//        logval H_mean = logval(0.0);
//        logval H_std = logval(0.0);

//        logval nr_covered_bins = logval(0.0);
        logval norm(0.0);
        for (int32 i = 0; i < static_cast<int32>(number_of_bins); i++)
        {
            // update weights
            if(lnhisto[i] > logval(0.5))
                weights[i] = weights[i]/lnhisto[i];
            norm += weights[i];

            if (i >= energy_index(E_min) && i <= energy_index(E_max))
            {

                if (lnhisto[i] > H_max)
                {
                    H_max = lnhisto[i];
                }
                
                if (lnhisto[i] > logval(0.5)) // Histogram not empty at i
                {
                    if (lnhisto[i] < H_min)
                    {
                        H_min = lnhisto[i];
                    }
//                    H_mean += lnhisto[i];
    //                H_std  += lnhisto[i]*lnhisto[i];
//                    nr_covered_bins += logval(1.0);
                }
            }
        }
//        H_mean /= nr_covered_bins;
        logval nr_covered_bins = logval(energy_index(E_max) - energy_index(E_min));
        logval H_mean = ln_number_of_measurements/nr_covered_bins;      
/*
        H_mean = ln_number_of_measurements/number_of_bins;
        H_std = (H_std - nr_covered_bins*H_mean*H_mean)/nr_covered_bins; // biased variance
        H_std = logval(std::sqrt(H_std.toDouble()));
        std::cout << "mean, var of H(E): " << H_mean.toDouble() << " , " << H_std.toDouble() << "\n";
        std::cout << "ten percent of mean: " << H_mean.toDouble()*0.1 << "\n";
        std::cout << "H_min: " << H_min.toDouble() << "\n";
        std::cout << "H_max: " << H_max.toDouble() << "\n";
*/
//        if ( H_max.value() - H_min.value() > 1.5 )
        if (H_min < H_mean*flatfactor 
            || H_max*flatfactor > H_mean)
        {
            flat_histogram = false;
            // std::cout << "Not flat.\n";
        }


        for (uint32 i=0; i < number_of_bins; i++)
        {
            weights[i] /= norm;
            lnhisto[i] /= ln_number_of_measurements;
        }
        
        number_of_iterations = parameters.get_muca_steps() - muca_steps; 

        if (flat_histogram)
        {
            muca_flatness = true;
        }
        if (energy_range)
        {
            muca_full_range = true;    
        }

# ifdef DEBUG
        std::cout << "# ";
        if (flat_histogram)
            std::cout << "flat ";
        if (energy_range)
            std::cout << "filled ";
        std::cout << std::endl;
# endif //DEBUG

        weights_file    << join(weights)        << std::endl;
        histo_file      << join(lnhisto)        << std::endl;

        muca_steps --;
    } //while
    weights_file.close();
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
    energy_histogram[total_energy]++;
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
    for (histogram_type::const_iterator it = energy_histogram.begin();
            it != energy_histogram.end(); ++it)
    {   
        out.unsetf(std::ios::floatfield);
        out << std::setw(20) << it->first;
        out.setf(std::ios::scientific, std::ios::floatfield);
        out.precision(16);
        out << std::setw(30) << it->second/weights[energy_index(it->first)];
        out << "\n";
    }
}

void mc_algorithm::weights_to_file(std::ofstream &out)
{
    out << "#" << std::setw(19) << "E" << std::setw(30) << "log(w(E))\n";
/*
    for (histogram_type::const_iterator it = energy_histogram.begin();
            it != energy_histogram.end(); ++it)
    {
        out.unsetf(std::ios::floatfield);
        out << std::setw(20) << it->first;
        out.setf(std::ios::scientific, std::ios::floatfield);
        out.precision(16);
        out << std::setw(30) << weights[energy_index(it->first)];
        out << "\n";
    }
*/

    for (size_t i=0; i < weights.size(); i++)
    {
        out.unsetf(std::ios::floatfield);
        out << std::setw(20) << static_cast<int32>(2*i) - minimal_energy_abs;
        out.setf(std::ios::scientific, std::ios::floatfield);
        out.precision(16);
        out << std::setw(30) << weights[i];
        out << "\n";
    }
}

void mc_algorithm::muca_energy_histogram_to_file(std::ofstream &out)
{
    out << "#" << std::setw(19) << "E" << std::setw(30) << "log(H_muca(E))\n";
    for (histogram_type::const_iterator it = energy_histogram.begin();
            it != energy_histogram.end(); ++it)
    {
        out.unsetf(std::ios::floatfield);
        out << std::setw(20) << it->first;
        out.setf(std::ios::scientific, std::ios::floatfield);
        out.precision(16);
        out << std::setw(30) << it->second;
        out << "\n";
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
