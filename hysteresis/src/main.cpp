/**
* @mainpage Monte Carlo Simulation of the Gonihedric Ising Model,
* first dual representation
*
* @author Marco Mueller
* @date 21.09.2011
* @note 
*
* The simulation uses a config file called 'config.ini' to set
* the physical and geometrical properties, other configuration
* filenames can be passed to the binary adding -c [config_filename]
* as a command line parameter\n
*\n
* Monte Carlo simulations are a probabilistic attempt to get the
* density of states of a physical model. Therefore they use 
a pseudo random numbers to update a given physical state to 
* minimize the energy.\n
*/
# include <iostream>
# include <iomanip>
# include <fstream>

// holds config and commandlineparser for retrieving relevant information
// about the simulation parameters
# include "simulation.hpp"
// calculates geometry of the lattice with different boundary conditions
# include "lattice.hpp"
# include "canon.hpp"

void print_header(std::ostream&, simulation&, mc_algorithm& );

int main(int argc, char** argv)
{
    time_t time_start = std::time(NULL);
    logger_ofstream.open("run.log");

    // initialize simulation environment
    // by parsing commandline, config files and setting all parameters
    simulation parameters(argc, argv);
    lattice geometry(parameters);
    mc_algorithm algorithm(parameters, geometry);

    print_header(logger_ofstream, parameters, algorithm);

    time_t time_init = std::time(NULL);
    double beta_min = parameters.get_beta_min();
    double beta_max = parameters.get_beta_max();
    uint32 nr_betas = parameters.get_nr_betas();
    double stepsize;
    if (nr_betas > 1){
      stepsize = (beta_max - beta_min)/static_cast<double>(nr_betas-1);
    }
    else {
      stepsize = (beta_max + beta_min)/2.0;
    }

    if(std::abs(parameters.get_p_init() - 0.5) < 0.1)
    // start with a disorder configuration and cool down
    {
        algorithm.set_beta(beta_min);
        algorithm.thermalisation_sweeps(parameters.get_thermalisation_sweeps());
    }
    else
    // start with an ordered configuration and heat up
    {
        algorithm.set_beta(beta_max);
        algorithm.thermalisation_sweeps(parameters.get_thermalisation_sweeps());
    }

    double beta;
    for (uint32 step = 0; step < nr_betas; step ++)
    {
        if(std::abs(parameters.get_p_init() - 0.5) < 0.1)
        // start with a disorder configuration and cool down
        {
            beta = beta_min + step*stepsize;
        }
        else
        // start with an ordered configuration and heat up
        {
            beta = beta_max - step*stepsize;
        }

        std::string fname = "energy_beta_" + std::to_string(beta) + ".dat";
        std::ofstream meas_file(fname);
        meas_file <<  std::scientific << std::setprecision(16);
        meas_file << "# beta = " << beta << std::endl;
        algorithm.set_beta(beta);
        algorithm.thermalisation_sweeps(parameters.get_tspb());
        algorithm.sweeps(meas_file, parameters.get_mspb());        
    }
    time_t time_end = std::time(NULL);
    
    logger_ofstream << "# initialization time (seconds): " 
        << time_init - time_start << "\n";
    logger_ofstream << "# simulation time (seconds): "
        << time_end - time_init << "\n";
    logger_ofstream << "# total time (seconds): "
        << time_end - time_start << "\n";

    logger_ofstream.close();

}
void print_header(
    std::ostream& out, simulation& parameters, mc_algorithm& algorithm)
{
    time_t rawtime;
    time ( &rawtime );
    out
        << "# " << ctime (&rawtime)
        << "#----------\n"
        << "# Parameters:\n"
        << "# random number generator = " << algorithm.get_rng_type() << "\n"
        << "# random number seed = "  << parameters.get_seed() << "\n";
    uint32_1d grid_dim; parameters.get_dimensions(grid_dim);
    out << "# grid dimensions = " << join(grid_dim) << "\n";
    out << "#----------" << std::endl;
}
