/**
 * @brief implementation of the Metropolis algorithm
 *
 * @author Marco Mueller
 * @date 05.10.2011
 *
 */

# ifndef __WANG_LANDAU_HPP__
# define __WANG_LANDAU_HPP__

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

        inline void sweep(std::ofstream&);

        void energy_histogram_to_file(std::ofstream &out);
        void weights_to_file(std::ofstream &out);

        void wl_energy_histogram_to_file(std::ofstream &out);
        uint32 get_number_of_iterations(){ return number_of_iterations; }
        double get_acceptance_ratio();

        const char * get_rng_type(){  return rng.get_type(); }

    private:
        simulation parameters;
        lattice geometry;
        hamiltonian energy;

        // random number generator
        rngclass rng;

        // inverse temperature
        double beta;
        
        // number of spins in the system
        uint32 number_of_spins;
        logval total_steps;
        logval accepted_steps;

        // micro states are stored in a 1d vector
        int32_1d sigma;
        int32_1d tau;

        void initialize_states();
        void calculate_weights();

        double factor; // factor for wang landau algorithm (f in the paper)
        
        double total_energy;

        typedef std::map<double, logval> histogram_type;
        histogram_type energy_histogram;
       
        // MuCa-Simulation stuff follows: 
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

        int32 number_of_iterations;

        std::ofstream weights_file;
        std::ofstream histo_file;
        // update that alters the weights
        inline void wl_update(uint32 site, bool is_sigma);
        // update without touching the weights
        inline void update(uint32 site, bool is_sigma);

        double accuracy;
        int32 iterations_test_flat;
};

mc_algorithm::mc_algorithm(
    simulation const &parameters_, lattice const &geometry_)
    : parameters(parameters_), geometry(geometry_), 
      energy(parameters_, geometry_), rng(parameters.get_seed()), 
      beta(parameters.get_beta()), number_of_spins(parameters.get_number_of_spins()),
      total_steps(0.0), accepted_steps(0.0),
      dim(parameters.get_dimension()), E_min(parameters.get_E_min()),
      E_max(parameters.get_E_max()), accuracy(parameters.get_wang_accuracy()),
      iterations_test_flat(parameters.get_wang_iterations_test_flatness())
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
    while (total_energy < E_min or total_energy > E_max);

    
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
            volume *= (dimensions.at(i) + 2);
        }
        minimal_energy_abs = 1.0/4.0*dim*(dim-1)*volume;
    }
    
    number_of_bins = ceil(minimal_energy_abs) + 1;
    maximal_energy_difference = 2*2*static_cast<int>(dim*(dim - 1));
     
    // histogram for every recursion step
    // using 1 as initialization introduces a statistical error,
    // but this should be neglectable for the calculation of weights
    lnhisto.resize(number_of_bins, logval(0.0));
    weights.resize(number_of_bins, logval(1.0));
    factor = exp(1);
    
    calculate_weights();
}

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

void mc_algorithm::calculate_weights()
{
    std::cout << "# wang landau density approximation started.\n";

    weights_file.open("weights.muca.time.dat");
    histo_file.open("histo.muca.time.dat");
    
    bool flat_histogram = false;
    uint32 iteration_counter = 0;
    while(factor > accuracy)
    {
        iteration_counter ++;
        std::cout << "Iteration " << iteration_counter << std::endl;

        int32 flat_test_counter = 0;

        std::fill(lnhisto.begin(), lnhisto.end(), logval(0.0));
        logval ln_updates = logval(0.0);

        while(!flat_histogram)
        {
            number_of_iterations ++;
            for(uint32 i = 0; i < 2*number_of_spins; i++)
            {
                // choose a site randomly (to avoid additional correlations)
                double r = rng.rand(); // in (0,1)
                uint32 site = static_cast<uint32>(floor(r * 2*number_of_spins));
                
                // update either sigma or tau
                if (site >= number_of_spins) // flip a sigma spin
                {
                    wl_update(site - number_of_spins, true);
                }
                else
                {
                    wl_update(site, false);
                }
                // Monte Carlo Update
                lnhisto[delta_energy_index(0)] ++;
                ln_updates ++;
            }

            // test every 10000 updates for flatness
            if(flat_test_counter >= iterations_test_flat)
            {
                flat_test_counter = 0;
                flat_histogram = true;
                logval H_min = ln_updates;
                logval H_max = logval(0.0);
                for (int32 i = energy_index(E_min); i < energy_index(E_max); i++)
                {
                    if (lnhisto[i] > logval(0.5) and lnhisto[i] < H_min)
                    {
                            H_min = lnhisto[i];
                    }
                    if (lnhisto[i] > H_max)
                    {
                        H_max = lnhisto[i];
                    }
                }
                if (H_min < H_max * 0.8)
                {
                    flat_histogram = false;
                }

                if (!flat_histogram)
                {
#               ifdef DEBUG
                    std::cout << " not flat when tested " << std::endl;
#               endif // DEBUG
                }
            } // flat_test_counter
            flat_test_counter ++;
        }
        factor = sqrt(factor);    
        flat_histogram = false;

        weights_file    << join(weights) << std::endl;
        histo_file      << join(lnhisto) << std::endl;
    }

    logval norm = logval(0.0);
    for (uint32 i = 0; i < number_of_bins; i++)
    {
        norm += weights.at(i);
    }

    for (uint32 i = 0; i < number_of_bins; i++)
    {
        weights.at(i)/=norm;
    }
    
    weights_file.close();
    histo_file.close();
    
    accepted_steps = logval(0.0);
    total_steps = logval(0.0);

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


void mc_algorithm::wl_energy_histogram_to_file(std::ofstream &out)
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


inline void mc_algorithm::sweep(std::ofstream &out)
{
    for(uint32 i = 0; i < 2*number_of_spins; i++)
    {
        // choose a site randomly (to avoid additional correlations)
        double r = rng.rand(); // in (0,1)
        uint32 site = static_cast<uint32>(floor(r * 2*number_of_spins));
        
        // update either sigma or tau
        if (site >= number_of_spins) // flip a sigma spin
        {
            update(site - number_of_spins, true);
        }
        else
        {
            update(site, false);
        }
    }
    energy_histogram[total_energy]++;
    out << total_energy << "\n";
}

inline void mc_algorithm::update(uint32 site, bool is_sigma)
{
    double dEloc = energy.get_local_energy_difference(site, sigma, tau, is_sigma);
    
    // check whether the new state is out of the energy interval of interest
    // reject if so or accept with wang landau probability if it is out of the
    // interval of interest 
    if (total_energy + dEloc - E_min > -1e-6 and total_energy + dEloc - E_max < 1e-6)
    {
        double r = rng.rand();
        logval frac = weights[delta_energy_index(dEloc)] / weights[delta_energy_index(0)];

        // acceptance probability = min(1.0, g(E_alt)/g(E_neu))
        // but r is in range (0,1] therefore the calculation can be omitted
        if (r <= frac.toDouble())
        {
            // accept
            if (is_sigma) sigma.at(site) *= -1;
            else tau.at(site) *= -1;
            total_energy += dEloc;
            accepted_steps++;
        }
    } 
    total_steps ++;
}

inline void mc_algorithm::wl_update(uint32 site, bool is_sigma)
{
    double dEloc = energy.get_local_energy_difference(site,sigma,tau,is_sigma);

    // check whether the new state is out of the energy interval of interest
    // reject if so or accept with wang landau probability if it is out of the
    // interval of interest 
    if (total_energy + dEloc - E_min > -1e-6 and total_energy + dEloc - E_max < 1e-6)
    {
        double r = rng.rand();
        logval frac = weights[delta_energy_index(dEloc)] / weights[delta_energy_index(0)];

        // acceptance probability = min(1.0, g(E_alt)/g(E_neu))
        // but r is in range (0,1] therefore the calculation can be omitted
        if (r <= frac.toDouble())
        {
            // accept
            if (is_sigma) sigma.at(site) *= -1;
            else tau.at(site) *= -1;
            total_energy += dEloc;
        }
    } 
    weights.at(delta_energy_index(0)) /= factor;
}

void mc_algorithm::initialize_states()
{
    if (!parameters.get_initial_state(sigma))
    {
        for (uint32 n = 0; n < number_of_spins; n++)
        {
            int32 spin = (rng.rand() < 0.5) ? -1 : 1;
            sigma.push_back(spin);

            spin = (rng.rand() < 0.5) ? -1 : 1;
            tau.push_back(spin);
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

double mc_algorithm::get_acceptance_ratio()
{
    logval ratio = accepted_steps/total_steps;
    return ratio.toDouble();
}


# endif
