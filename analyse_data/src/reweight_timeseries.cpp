# include <dlib/cmd_line_parser.h>
# include <iostream>
# include <iomanip>
# include <sstream>
# include "readfile.hpp"
# include "discrete_function.hpp"
# include "logval.hpp"
# include "gzstream.h"
# include "running_weighted.hpp"

using std::cout;
using std::cin;
using std::string;
using std::stringstream;

int main(int argc, char** argv){
    typedef dlib::cmd_line_parser<char>::check_1a_c clp_type;
    clp_type clp;
    
    clp.add_option("h", "help", 0);
    clp.add_option("w", "file with LOGARITHMIC (multicanonical/Wang-Landau) weights", 1);
    clp.add_option("sep", "separator (\"tab\", \"white\", or some string)", 1);
    clp.add_option("b0", "inverse simulation temperature (0.0)", 1);

    // TODO: one could also implement a file with beta-values instead which
    // would be more flexible
    clp.add_option("g", "assumes the input file is gzipped", 0);

    clp.parse(argc, argv);

    if (clp.option("h") || clp.number_of_arguments() < 1)
    {
        cout << "Usage: " << argv[0] << " filename [options]\n";
        clp.print_options(cout); 
        cout << std::endl;
        return EXIT_SUCCESS;
    }

    // define default values for the command line options:
    // separator
    const string sep_def(" ");
    double beta_sim = dlib::get_option(clp, "b0", 0.0);

    string filename      = clp[0];
    string sep           = dlib::get_option(clp, "sep", sep_def);
    if ( sep == "tab" ) sep = "\t";
    if ( sep == "white" ) sep = " ";

    discrete_function<double> logw(0,1,1);
    if ( clp.option("w") ){
        string weights_filename = clp.option("w").argument();
        logw = from_file<double>(weights_filename);
    }

    const size_t colE = 0; // column of energy
    const size_t colO = 0; // column of observable of interest

    data_file<double> ts;
    bool read_file = false;
    if ( clp.option("g") ){ // read gzstream
      stringstream sstr;
      igzstream infilegz(filename.c_str());
      sstr << infilegz.rdbuf();
      if (!sstr.good()) {
        std::cerr << "Error reading timeseries file.\n";
        return EXIT_FAILURE;
      }
      read_file = ts.read_data_2d(sstr, sep);
    }
    else {
      read_file = ts.read_data_2d(filename, sep);
    }

    if ( !read_file ){ 
      std::cerr << "Error reading timeseries file.\n";
      return EXIT_FAILURE;
    }

    cout << "# beta ... inverse temperature\n";
    cout << "# E ... internal energy\n";
    cout << "# C ... heat capacity ... beta**2*var(E)\n";
    cout << "# B ... Binder's energy cumulant\n";
    cout << "# beta E Cv B\n";
    cout << std::scientific << std::setprecision(15);

    // read inverse temperatures from stdin
    size_t nr_skipped = 0;
    while (!cin.eof()) {
      string line;
      getline(cin, line);

      if (cin.fail()) {
        // error
        break;
      }
      // ignore empty lines and comments
      try {
        const double beta = std::stod(line);
        const double dbeta = beta - beta_sim;
        statistics::running_weighted<logval> O1;
        statistics::running_weighted<logval> O2;
        statistics::running_weighted<logval> O4;
        for (size_t t = 0; t < ts.data[colE].size(); ++t){
          const double Et = ts.data[colE][t];
          const logval w = make_logval_from(1, -logw[Et] - dbeta*Et);
          const double meas = ts.data[colO][t];
          const logval obs1(meas);
          const logval obs2 = obs1*obs1;
          const logval obs4 = obs2*obs2;
          O1(obs1, w);
          O2(obs2, w);
          O4(obs4, w);
        }
        cout << beta 
          << " " << O1.mean()
          << " " << beta*beta*O1.variance_n()
          << " " << 1.0 - O4.mean()/(3.0*O2.mean()*O2.mean())
          << "\n";
      }
      catch ( ... ){
        nr_skipped++;
      }
    } // while we read from stdin
    if (nr_skipped){
      if ( nr_skipped == 1){
        std::cerr << "# Warning: skipped one ill-formatted line.\n";
      }
      else {
        std::cerr << "# Warning: skipped " << nr_skipped << 
          " ill-formatted lines.\n";
      }
    }

	return 0;
}
