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
        
        // == physical variables
        // lattice dimensions
        uint32_1d dimensions;
        // total number of spins in the lattice
        uint32 number_of_spins;
        
        // == simulation parameters
        // number of sweeps for initial thermalisation
        uint32 thermalisation_sweeps;
        // random number generator seed
        uint32 rng_seed;
        // interval of betas to be measured
        double beta_min;
        double beta_max;
        // number of betas in the interval [beta_min, beta_max]
        uint32 nr_betas; 
        // thermalisation sweeps per beta
        uint32 tspb;
        // measurement sweeps per beta
        uint32 mspb;
        // probability of a spin being "up" in the initial configuration
        double p_init;

    public:
        // constructor
        simulation(int argc, char** argv);
        
        // getter methods
        void    get_dimensions(uint32_1d&);
        uint32  get_number_of_spins();
        uint32  get_dimension();
        uint32  get_seed();
        double  get_beta_min();
        double  get_beta_max();
        double  get_p_init();
        uint32  get_thermalisation_sweeps();
        uint32  get_nr_betas();
        inline uint32  get_tspb();
        inline uint32  get_mspb();
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
        std::ifstream config;
        config.open(configuration_filename.c_str());
        std::string s;
        while(getline(config, s))
            logger_ofstream << "# " << s << "\n";
   }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

// ***
// getter - methods
double  simulation::get_beta_min()  { return beta_min; }
double  simulation::get_beta_max()  { return beta_max; }
double  simulation::get_p_init()    { return p_init; }
uint32  simulation::get_nr_betas()  { return nr_betas; }
inline uint32  simulation::get_tspb()      { return tspb; }
inline uint32  simulation::get_mspb()      { return mspb; }
uint32  simulation::get_thermalisation_sweeps() { return thermalisation_sweeps; }
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
uint32 simulation::get_number_of_spins() { return number_of_spins; }

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
    
    rng_seed = static_cast<uint32>(std::time(0));
    if(cr["random_numbers_seed"] != "")
    {
        rng_seed = dlib::string_cast<uint32>(cr["random_numbers_seed"]);
    } // if random numbers seed

    // number of sweeps for initial thermalisation
    thermalisation_sweeps = 100000;
    if (cr["thermalisation_sweeps"] != "")
    {
        thermalisation_sweeps = 
            dlib::string_cast<uint32>(cr["thermalisation_sweeps"]);
    }
    // interval of betas to be measured
    beta_min = 1.25;
    if (cr["beta_min"] != "")
    {
        beta_min = 
            dlib::string_cast<double>(cr["beta_min"]);
    }
    beta_max = 1.39;
    if (cr["beta_max"] != "")
    {
        beta_max = 
            dlib::string_cast<double>(cr["beta_max"]);
    }
    // number of betas in the interval [beta_min, beta_max]
    nr_betas = 32; 
    if (cr["nr_betas"] != "")
    {
        nr_betas = 
            dlib::string_cast<uint32>(cr["nr_betas"]);
    }
    // thermalisation sweeps per beta
    tspb = 10000;
    if (cr["tspb"] != "")
    {
        tspb = 
            dlib::string_cast<uint32>(cr["tspb"]);
    }
    // measurement sweeps per beta
    mspb = 100000;
    if (cr["mspb"] != "")
    {
        mspb = 
            dlib::string_cast<uint32>(cr["mspb"]);
    }

    p_init = 0.5;
    if (cr["p_init"] != "")
    {
        p_init = 
            dlib::string_cast<double>(cr["p_init"]);
    }
    
}

# endif // __SIMULATION_HPP__
