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
        inline void sweeps(std::ofstream&, uint32);
        inline void update(uint32, bool);
        inline void update_no_measurements(uint32, bool);
        inline void thermalisation_sweeps(uint32);

        const char * get_rng_type(){  return rng.get_type(); }
        void set_beta(double);

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
        double minimal_energy_abs;
        uint32 dim;
};

mc_algorithm::mc_algorithm(simulation const &parameters_, lattice const &geometry_)
    : parameters(parameters_), geometry(geometry_), energy(parameters_, geometry_),
      rng(parameters.get_seed()), beta(0.0), accepted_steps(0.0), total_steps(0.0),
      number_of_spins(parameters.get_number_of_spins())
{
    initialize_states();
    total_energy = energy.calculate_total_energy(sigma, tau);
    dim = parameters.get_dimension();
    minimal_energy_abs = 1.0/4.0*dim*(dim-1)*number_of_spins;
}

double mc_algorithm::get_acceptance_ratio()
{
    logval ratio = accepted_steps/total_steps;
    return ratio.toDouble();
}

inline void mc_algorithm::sweep(std::ofstream &out)
{
    for (uint32 i = 0; i < 2*number_of_spins; i++)
    {
        // choose a site randomly (to avoid additional correlations)
        double r = rng.rand(); // in (0,1)
        uint32 site = static_cast<uint32>(floor(r * 2*number_of_spins));

        // Update either sigma or tau
        if (site >= number_of_spins) // flip a sigma spin
        {
            update(site - number_of_spins, true);
        }
        else
        {
            update(site, false);
        }
    }
    out << total_energy << "\n";
}

inline void mc_algorithm::sweeps(std::ofstream &out, uint32 nr_sweeps)
{
    for (uint32 k = 0; k < nr_sweeps; k++)
    {
        sweep(out);
    }
}

inline void mc_algorithm::thermalisation_sweeps(uint32 nr_sweeps)
{
    for (uint32 k = 0; k < nr_sweeps; k++)
    {
        for (uint32 i = 0; i < 2*number_of_spins; i++)
        {
            // choose a site randomly (to avoid additional correlations)
            double r = rng.rand(); // in (0,1)
            uint32 site = static_cast<uint32>(floor(r * 2*number_of_spins));

            // Update either sigma or tau
            if (site >= number_of_spins) // flip a sigma spin
            {
                update_no_measurements(site - number_of_spins, true);
            }
            else
            {
                update_no_measurements(site, false);
            }
        }
    }
}

inline void mc_algorithm::set_beta(double beta_)
{
    beta = beta_;
}

inline void mc_algorithm::update_no_measurements(uint32 site, bool is_sigma)
{
    double dEloc = energy.get_local_energy_difference(site, sigma, tau, is_sigma);
    double r = rng.rand();
    double p_flip = exp(-beta*dEloc);
    // acceptance probability = min(p_muca, 1.0)
    // but r is in range (0,1] therefore the calculation can be omitted
    if (r < p_flip)
    {
        // accept
        if (is_sigma) sigma[site] *= -1;
        else tau[site] *= -1;

        total_energy += dEloc;
        // spin was already flipped --> (1)
    }
}

inline void mc_algorithm::update(uint32 site, bool is_sigma)
{
    double dEloc = energy.get_local_energy_difference(site, sigma, tau, is_sigma);
    double r = rng.rand();
    double p_flip = exp(-beta*dEloc);
    // acceptance probability = min(p_muca, 1.0)
    // but r is in range (0,1] therefore the calculation can be omitted
    if (r < p_flip)
    {
        // accept
        if (is_sigma) sigma[site] *= -1;
        else tau[site] *= -1;

        total_energy += dEloc;
        // spin was already flipped --> (1)
        accepted_steps++;
    }
    total_steps ++;
}

void mc_algorithm::initialize_states()
{
    double p_init = parameters.get_p_init();
    for (uint32 n = 0; n < number_of_spins; n++)
    {
        int32 spin;

        spin = (rng.rand() < p_init) ? -1 : 1;
        sigma.push_back(spin);

        spin = (rng.rand() < p_init) ? -1 : 1;
        tau.push_back(spin);

    }
}

# endif
