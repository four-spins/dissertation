/**
 * @brief calculates the hamiltonian of the gonihedric ising model
 * 
 *  \f[
 *  H = -0.5\kappa\sum_{<i,j>_x}\sigma_i\sigma_j + 
 *      \ -0.5\kappa\sum_{<i,j>_y}\sigma_i\sigma_j + 
 *      \ -0.5\kappa\sum_{<i,j>_z}\sigma_i\sigma_j
 *  \f]
 * 
 * @author Marco Mueller
 * @date 11.05.2011
 *
 */

# ifndef __HAMILTONIAN_HPP__
# define __HAMILTONIAN_HPP__

// include constant stuff
# include "const.hpp"
// simulation parameters
# include "simulation.hpp"
// neighbour tables
# include "lattice.hpp"

class hamiltonian
{
    public:
        hamiltonian(simulation const &, lattice const &);
        /// calculates total energy of the system
        inline double calculate_total_energy(int32_1d &, int32_1d &);
        /// calculates local energy difference when changing the sign 
        /// of a specific spin
        inline double get_local_energy_difference(uint32, int32_1d &, int32_1d &, bool);

    private:
        simulation parameters;
        lattice geometry;

        double beta;
        uint32 number_of_spins;
        
};

hamiltonian::hamiltonian
    (simulation const &parameters_, lattice const &geometry_)
    : parameters(parameters_), geometry(geometry_), 
      beta(parameters.get_beta()), 
      number_of_spins(parameters.get_number_of_spins())
{
}


/// returns the hamiltonian of the system
inline double hamiltonian::calculate_total_energy
    (int32_1d& sigma, int32_1d& tau)
{
    double energy = 0;

    // store the interaction sums seperately
    double nearx = 0;
    double neary = 0;
    double nearz = 0;

    // for every spin in the lattice
    for (uint32 n = 0; n < number_of_spins; n++)
    {
        uint32 nnrx = geometry.get_near_neighbour(n, 0);
        uint32 nnlx = geometry.get_near_neighbour(n, 1);
        uint32 nnry = geometry.get_near_neighbour(n, 2);
        uint32 nnly = geometry.get_near_neighbour(n, 3);
        uint32 nnrz = geometry.get_near_neighbour(n, 4);
        uint32 nnlz = geometry.get_near_neighbour(n, 5);

        nearx += sigma[nnrx]*sigma[n]; 
        nearx += sigma[nnlx]*sigma[n]; 
        neary += tau[nnry]*tau[n]; 
        neary += tau[nnly]*tau[n]; 
        
        nearz += sigma[nnrz]*sigma[n] * tau[nnrz]*tau[n]; 
        nearz += sigma[nnlz]*sigma[n] * tau[nnlz]*tau[n]; 

        // fixed boundary conditions force the interaction on the boundary
        // to be calculated twice (as the Hamiltonian counts only interactions
        // and divides by a factor of two)
        if (nnrx >= number_of_spins) nearx += sigma[nnrx]*sigma[n];
        if (nnlx >= number_of_spins) nearx += sigma[nnlx]*sigma[n];
        if (nnry >= number_of_spins) neary += tau[nnry]*tau[n];
        if (nnly >= number_of_spins) neary += tau[nnly]*tau[n];
        if (nnrz >= number_of_spins) 
            nearz += sigma[nnrz]*sigma[n] * tau[nnrz]*tau[n]; 
        if (nnlz >= number_of_spins) 
            nearz += sigma[nnlz]*sigma[n] * tau[nnlz]*tau[n]; 
    }
    
    // correction of prefactors:
    // every nearest neighbour interaction was counted twice,
    nearx /= 2.0;
    neary /= 2.0;
    nearz /= 2.0;

    energy -= 0.5*(nearx + neary + nearz);

    return energy;
}

inline double hamiltonian::get_local_energy_difference(uint32 n, int32_1d &sigma, int32_1d &tau, bool is_sigma)
{
    // calculates the energy difference when flipping a spin on site n
    // if is_sigma == true the energy difference when flipping sigma is delivered
    // if is_sigma == false the energy difference when flippign tau is deliviered
    double dE = 0;

    // store the interaction sums seperately
    double nearx = 0;
    double neary = 0;
    double nearz = 0;

    uint32 nnrx = geometry.get_near_neighbour(n, 0);
    uint32 nnlx = geometry.get_near_neighbour(n, 1);
    uint32 nnry = geometry.get_near_neighbour(n, 2);
    uint32 nnly = geometry.get_near_neighbour(n, 3);
    uint32 nnrz = geometry.get_near_neighbour(n, 4);
    uint32 nnlz = geometry.get_near_neighbour(n, 5);

    nearx += sigma[nnrx]*sigma[n]; 
    nearx += sigma[nnlx]*sigma[n]; 
    neary += tau[nnry]*tau[n]; 
    neary += tau[nnly]*tau[n]; 
    nearz += sigma[nnrz]*sigma[n] * tau[nnrz]*tau[n]; 
    nearz += sigma[nnlz]*sigma[n] * tau[nnlz]*tau[n]; 
    
    if (is_sigma)
    {
        dE = nearx + nearz;
    }
    else
    {
        dE = neary + nearz;
    }

    return dE;
}

# endif
