# ifndef __SIMULATION_HPP__
# define __SIMULATION_HPP__

// --- my own includes
// type definitions
# include "const.hpp"

// --- dlib C++ includes
// command line parser
# include "dlib/cmd_line_parser.h"
// configuration file parser
# include "dlib/config_reader.h"
// string casts, concatenation, trimming (for config file parsing)
# include "dlib/string.h"
//misc_api_kernel_abstract.h"

// --- Standard C++ includes
// file operations using streams
# include <fstream>
// for seeding random numbers, std::time will be used
# include <ctime>

# include "logval.hpp"

typedef dlib::cmd_line_parser<char>::check_1a_c clp;
typedef dlib::config_reader::kernel_1a cr_type;

class simulation
{
    private:
        // == technical part
        std::string configuration_filename;
        // command line parser
        void    parse_commandline(int argc, char** argv);
        // configuration file parser
        void    parse_configuration_file();
        // reads a vector from a file
//        template <class T>
//             void read_vector_from_file(std::vector<T> &, std::string);
        void read_weights_from_file(std::vector<logval> &, std::vector<double> &, std::string);
        
        // == physical variables
        // lattice dimensions
        uint32_1d dimensions;
        // total number of spins in the lattice
        uint32 number_of_spins;
        // boundary conditions
        uint8 boundary_conditions;
        // interaction weight
        double kappa;
        // inverse temperature
        double beta;
        
        // == simulation parameters
        // number of sweeps for measurements
        int32 sweeps;
        // number of steps to conduct in a multicanonical simulation
        // to calculate the weights
        int32 muca_steps;
        // random number generator seed
        uint32 rng_seed;
        // initial microstate can be read in
        int32_1d initial_state;
        // weights if weights_file is given
        std::vector<logval> weights;
        std::vector<double> energies;
        bool weights_came_from_file;
        bool weights_over_time;

        // int32 number_of_bins;
        double E_min;
        double E_max;

        double wang_accuracy;
        int32 wang_iterations_test_flatness;

        int32 max_tunnel_events;

    public:
        // constructor
        simulation(int argc, char** argv);
        
        // getter methods
        void    get_dimensions(uint32_1d&);
        void    get_weights(std::vector<logval> &, std::vector<double> &);
        bool    weights_from_file();
        bool    get_weights_over_time();
        double  get_kappa();
        double  get_beta();
        uint32  get_number_of_spins();
        uint32  get_dimension();
        uint32  get_seed();
        int32   get_sweeps();
        int32   get_muca_steps();
        uint8   get_boundary_conditions();
        bool    get_initial_state(int32_1d&);
        int32   get_wang_iterations_test_flatness();
        double  get_wang_accuracy();
        int32   get_max_tunnel_events();

        // int32   get_number_of_bins();
        double  get_E_min();
        double  get_E_max();
};

/// Constructor for setting the simulation environment.\n
/// Parses commandline and configuration file
simulation::simulation(int argc, char** argv)
{
    try
    {
        // command line parsing
        parse_commandline(argc, argv);

        // configuration file reading
        parse_configuration_file();

   }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

// ***
// getter - methods
double  simulation::get_wang_accuracy() { return wang_accuracy; }
int32   simulation::get_wang_iterations_test_flatness() { return wang_iterations_test_flatness; }
int32   simulation::get_max_tunnel_events() { return max_tunnel_events; }
double  simulation::get_E_min() { return E_min; }
double  simulation::get_E_max() { return E_max; }
// int32   simulation::get_number_of_bins() { return number_of_bins; }
int32   simulation::get_sweeps()    { return sweeps; }
int32   simulation::get_muca_steps() { return muca_steps; }
double  simulation::get_kappa()     { return kappa; }
double  simulation::get_beta()      { return beta; }
uint32  simulation::get_seed()      { return rng_seed; }
uint32  simulation::get_dimension()
{
    return static_cast<uint32>(dimensions.size());
}

void simulation::get_dimensions(uint32_1d& return_vector)
{
    return_vector.assign
    (
        dimensions.begin(),
        dimensions.end()
    );
}

inline
void simulation::get_weights(std::vector<logval>& return_vector, std::vector<double> &energies_)
{
    return_vector.assign
    (
        weights.begin(),
        weights.end()
    );
    energies_.assign
    (
        energies.begin(),
        energies.end()
    );
}

bool simulation::weights_from_file() { return weights_came_from_file; }
bool simulation::get_weights_over_time() { return weights_over_time; }

bool simulation::get_initial_state(int32_1d& return_vector)
{
    return_vector.assign
    (
        initial_state.begin(),
        initial_state.end()
    );
    if (initial_state.size() > 0)
        return true;

    return false;
}

uint8 simulation::get_boundary_conditions()
{
    return boundary_conditions;
}

uint32 simulation::get_number_of_spins()
{
    return number_of_spins;
}

void simulation::parse_commandline(int argc, char** argv)
{
    clp parser;
    parser.add_option("h","Display this help message.", 0);
    
    // parse the commandline
    parser.parse(argc, argv);

    if (parser.option("h") or parser.number_of_arguments() != 1)
    {
        std::cout << "Usage: " << argv[0] << " configuration_file\n";
        parser.print_options(std::cout);
        std::cout << std::endl;
        exit(0);
    }
    
    configuration_filename = parser[0];
}   

void simulation::parse_configuration_file()
{
    cr_type cr(configuration_filename);
    // read content of configuration file

    if(cr["dimensions"] != "")
    {
        std::string full_string = cr["dimensions"];
        // split string at commas
        while(full_string.find(",",0) != std::string::npos)
        {
            // store the position of delimiter
            size_t pos = full_string.find(",",0);
            // copy the substring in front of the first comma
            std::string temp = full_string.substr(0,pos);
            // erase the substring from source string
            full_string.erase(0, pos + 1);
            // put the number into an array after deleting all whitespaces
            dimensions.push_back(dlib::string_cast<int>(dlib::trim(temp)));
        }
        // add the last dimension
        dimensions.push_back(
                dlib::string_cast<int>(dlib::trim(full_string))
                );
        number_of_spins = 1;
        for (std::size_t i = 0; i < dimensions.size(); i++)
        {
            number_of_spins *= dimensions.at(i);
        }
    } // if cr["dimensions"]
    
    // parse the boundary condition
    // note that the constants are defined in "const.hpp"
    if(cr["boundary_conditions"] != "")
    {
        std::string bc = dlib::trim(cr["boundary_conditions"]);
        if(bc == "periodic")
        {
            boundary_conditions = PERIODIC;
        }
        if(bc == "interface")
        {
            boundary_conditions = INTERFACE;
        }
        if(bc == "vanish")
        {
            boundary_conditions = VANISH;
        }
        if(bc == "helical")
        {
            boundary_conditions = HELICAL;
        }
    }

    if(cr["state_file"] != "")
    {
        std::ifstream state_file;
        state_file.open(cr["state_file"].c_str());
        if (!state_file)
        {
            std::cerr << "Error opening file " << cr["state_file"] << "\n";
            exit(1);
        }
        std::string state;
        while ( std::getline(state_file, state) )
        {
            while(state.find(",",0) != std::string::npos)
            {
                size_t pos = state.find(",",0);
                std::string temp = state.substr(0,pos);
                state.erase(0, pos + 1);
                initial_state.push_back(
                        dlib::string_cast<int32>(dlib::trim(temp)));
            }
            initial_state.push_back(
                    dlib::string_cast<int32>(dlib::trim(state)));
        }
        state_file.close();

        assert(initial_state.size() == number_of_spins);
    }
    
    // muca_steps = number_of_spins*10;
    muca_steps = number_of_spins*number_of_spins;
    if (cr["muca_steps"] != "")
    {
        muca_steps = dlib::string_cast<int32>(cr["muca_steps"]);
    }

    rng_seed = static_cast<uint32>(std::time(0));
    if(cr["random_numbers_seed"] != "")
    {
        rng_seed = dlib::string_cast<uint32>(cr["random_numbers_seed"]);
    } // if random numbers seed
    
//    kappa = 0.0;
//    if(cr["kappa"] != "")
//        kappa = dlib::string_cast<double>(cr["kappa"]);

    beta = 0.0;
    if(cr["beta"] != "")
        beta = dlib::string_cast<double>(cr["beta"]);
    
    sweeps = 1000000;
    if(cr["number_of_sweeps"] != "")
        sweeps = dlib::string_cast<int32>(cr["number_of_sweeps"]);
    
    weights_came_from_file = false;
    if(cr["weights_file"] != "")
    {
        weights_came_from_file = true;
        read_weights_from_file(weights, energies, cr["weights_file"]);
    }
    
    // TODO: minimal and maximal energy should be kappa - dependant

    int32 dim = get_dimension(); 
    E_min = - 1.0/4.0*dim*(dim-1)*number_of_spins;
    if(cr["E_min"] != "")
    {
        E_min = dlib::string_cast<double>(cr["E_min"]);
//        E_min *= static_cast<double>(number_of_spins);
        int32 Eminint = static_cast<int>(number_of_spins * E_min);
        if (std::abs(Eminint) % 2 == 1)
           Eminint += 1;
        E_min = static_cast<double>(Eminint);
    }
    E_max = 1.0/4.0*dim*(dim-1)*number_of_spins;
    if(cr["E_max"] != "")
    {
        E_max = dlib::string_cast<double>(cr["E_max"]);
        int32 Emaxint = static_cast<int>(number_of_spins * E_max);
        if (std::abs(Emaxint) % 2 == 1)
           Emaxint -= 1;
        E_max = static_cast<double>(Emaxint);
    }
    
    wang_iterations_test_flatness = 10000;
    if(cr["wang_iterations_test_flatness"] != "")
    {
        wang_iterations_test_flatness = 
            dlib::string_cast<int32>(cr["wang_iterations_test_flatness"]);
    }
    
    wang_accuracy = std::exp(pow(10.0, -8));
    if(cr["wang_accuracy"] != "")
    {
        wang_accuracy = 
            std::exp(std::pow(10.0, dlib::string_cast<int>(cr["wang_accuracy"])));
    }

    max_tunnel_events = 8;
    if(cr["max_tunnel_events"] != "")
    {
        max_tunnel_events = 
            dlib::string_cast<int32>(cr["max_tunnel_events"]);
    }

    weights_over_time = false;
    if(cr["weights_over_time"] != "")
    {
        int foo =
            dlib::string_cast<int32>(cr["weights_over_time"]);
        if (foo != 0)
            weights_over_time = true;
    }

/*
    number_of_bins = ceil(std::abs(E_min_muca)) + 1;
    if(cr["number_of_bins"] != "")
        number_of_bins = dlib::string_cast<int32>(cr["number_of_bins"]);
*/  
    //
    //TODO:     add other stuff here 
    //          (Monte Carlo Sweeps, Thermalization, ...)
}

void simulation::read_weights_from_file(std::vector<logval> &weights_, std::vector<double> &energies_, std::string filename)
{
    std::ifstream vector_file(filename.c_str());
    if ( !vector_file.is_open() )
    {
        std::cerr << "Error opening file '" << filename << "'. Aborting.\n";
        std::exit(1);
    }
    else
    {
        weights_.clear();
        std::string lineread;
        while(std::getline(vector_file, lineread))
        {   
            // skip comments
            if(lineread.at(0) != '#')
            {
                //std::cout << lineread << "\n";
                std::istringstream is(lineread);
                double E; 
                logval weight;
                is >> E >> weight;
                // std::cout << E << "\t" << weight << std::endl;
                weights_.push_back(weight);
                energies_.push_back(E);
            }
        }
    }
}

# endif // __SIMULATION_HPP__
