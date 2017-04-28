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
*/

# include <unistd.h>
// created by the makefile, stores the git hash and time of compilation
# include "version.hpp"
// holds config and commandlineparser for retrieving relevant information
// about the simulation parameters
# include "simulation.hpp"
// calculates geometry of the lattice with different boundary conditions
# include "lattice.hpp"

# ifdef WANGLANDAU
#   include "wang_landau.hpp"
# endif

# ifdef MUCA
#   include "muca.hpp"
# endif

# ifdef MUCA_SIMPLE
#   include "muca_simple.hpp"
#   define MUCA
# endif


void print_header(std::ostream&, simulation&, mc_algorithm& );

int main(int argc, char** argv)
{
    time_t time_start = std::time(NULL);
    std::cout << "# Process id: " << getpid() << std::endl;
    
    // initialize simulation environment
    // by parsing commandline, config files and setting all parameters
    simulation parameters(argc, argv);
    lattice geometry(parameters);

    time_t time_init = std::time(NULL);
    mc_algorithm algorithm(parameters, geometry);

    // collect data and store timeseries:
    std::cout << "# measurements started" << std::endl;
    std::ofstream timeseries_ofstream;
    timeseries_ofstream.open("energy_time.dat");
    print_header(timeseries_ofstream, parameters, algorithm);
    for (int32 i = 0; i < parameters.get_sweeps(); i++)
    {
        algorithm.sweep(timeseries_ofstream);
    }
    timeseries_ofstream.close();
    time_t time_end = std::time(NULL);
 
    // some statistics
    std::ofstream logger_ofstream;
    logger_ofstream.open("simulation.log");
    print_header(logger_ofstream, parameters, algorithm);
    logger_ofstream << "# initialization time (seconds): " 
        << time_init - time_start << "\n";
    logger_ofstream << "# simulation time (seconds): "
        << time_end - time_init << "\n";
    logger_ofstream << "# total time (seconds): "
        << time_end - time_start << "\n";
    logger_ofstream << "# acceptance ratio: " 
        << algorithm.get_acceptance_ratio() << "\n";
    

    std::ofstream histogram_ofstream;
    histogram_ofstream.open("energy_histogram.dat");
        print_header(histogram_ofstream, parameters, algorithm);
        algorithm.energy_histogram_to_file(histogram_ofstream);
    histogram_ofstream.close();

# ifdef WANGLANDAU
    std::ofstream weights_ofstream;
    weights_ofstream.open("weights.muca.dat");
        print_header(weights_ofstream, parameters, algorithm);
        algorithm.weights_to_file(weights_ofstream);
    weights_ofstream.close();

    std::ofstream wl_histogram_ofstream;
    wl_histogram_ofstream.open("energy_histogram.muca.dat");
        print_header(wl_histogram_ofstream, parameters, algorithm);
        algorithm.wl_energy_histogram_to_file(wl_histogram_ofstream);
    wl_histogram_ofstream.close();
# endif

# ifdef MUCA
    std::ofstream weights_ofstream;
    weights_ofstream.open("weights.muca.dat");
        print_header(weights_ofstream, parameters, algorithm);
        algorithm.weights_to_file(weights_ofstream);
    weights_ofstream.close();

    std::ofstream muca_histogram_ofstream;
    muca_histogram_ofstream.open("energy_histogram.muca.dat");
        print_header(muca_histogram_ofstream, parameters, algorithm);
        algorithm.muca_energy_histogram_to_file(muca_histogram_ofstream);
    muca_histogram_ofstream.close();

# endif
}

void print_header(
    std::ostream& out, simulation& parameters, mc_algorithm& algorithm)
{

# ifdef WANGLANDAU
    out << "# Wang Landau Sampling:\n";
    out << "# number of iterations: " 
        << algorithm.get_number_of_iterations() << "\n";
    out << "# final factor: " << parameters.get_wang_accuracy() << "\n";
    out << "# sweeps between test of flatness: "
        << parameters.get_wang_iterations_test_flatness() << "\n";
# endif

# ifdef MUCA
    out << "# Multicanonical Simulation:\n";
    out << "# --- convergency of weights:\n";
    if (parameters.weights_from_file())
        out << "# weights read from file\n";
    if (parameters.get_muca_steps() == 0)
        out << "# no multicanonical steps, cannot determine quality of weights.\n";
    else
    {
        if (algorithm.converged_range())
            out << "# hardcoded energy range filled, ";
        else
            out << "# hardcoded energy range not filled, ";
        if (algorithm.converged_flatness())
            out << "flat";
        else
            out << "not flat";
    }
    out << "# \n";

    out << "# number of iterations: " 
        << algorithm.get_number_of_iterations() << "\n";
# endif

    out << "# --- git:\n";
    out << "# git hash = " << GIT_HASH << "\n";
    if (!GIT_COMMITED)
        out << "# *** UNCOMMITED CHANGES *** \n#\n";

    time_t rawtime;
    time ( &rawtime );
    out
        << "# " << ctime (&rawtime)
        << "#----------\n"
        << "# Parameters:\n"
        << "# random number generator = " << algorithm.get_rng_type() << "\n"
        << "# random number seed = "  << parameters.get_seed() << "\n"
        << "# beta = " << parameters.get_beta() << "\n"
        << "# boundary conditions = ";
    if (parameters.get_boundary_conditions() == PERIODIC)
        out << " periodic\n";
    else if(parameters.get_boundary_conditions() == INTERFACE)
        out << " interface\n";
    else if(parameters.get_boundary_conditions() == VANISH)
        out << " vanish\n";
    else if(parameters.get_boundary_conditions() == HELICAL)
        out << " helical\n";

        out << "# thermalization sweeps = " << 0 << "\n";
        out << "# sweeps = " << parameters.get_sweeps() << "\n";
        //<< "# metropolis sweeps per parallel tempering attempt = "
        //    << met_sweeps << "\n"
    uint32_1d grid_dim; parameters.get_dimensions(grid_dim);
    out << "# grid dimensions = " << join(grid_dim) << "\n";
    out << "#----------\n";
}

