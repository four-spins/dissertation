/******************************************************************************
 *
 * @file: readfile.hpp
 *
 * @date: 09/20/2012 01:06:11 PM (CEST)
 *
 * @author: Marco Müller <muelma@gmail.com>
 *
 ******************************************************************************/
# ifndef READFILE_HPP 
# define READFILE_HPP 

# include <iostream>
# include <fstream>
# include <cstdlib>
# include <cassert>
# include <sstream>
# include <iterator> // istream_iterator
# include <map>
# include <vector>
# include "string.hpp"

/// @brief Represents a data file with some meta-informations.
///
/// Meta information is contained in a header, where each line
/// is marked by a hash. The format can be used for a histogram or time series.
/// This class, however assumes that the file can be read into RAM.
///
/// example:
/// \verbatim
/// # simulation name
/// # git revision: GITHASH 
/// # time: TIMESPECIFIER
/// # additional stuff like parameters of the simulation,
/// # information to algorithms, seed, ...
/// # 
/// # parameter1 = 0.98321
/// # parameter2 = 42
/// #
/// # Energy Magnetisation MagnetisationSquared ... EnergyDensity
/// 3.134 1.559 2.265 ... 0.497
///  :::   :::   :::  ...  :::
///  :::   :::   :::  ...  :::
///  :::   :::   :::  ...  :::
/// \endverbatim
///
/// The last line tagged with the comment character (#) is assumed to hold 
/// the names of the observables. If this is not the case, variable names 
/// will become a string of the column number. 
/// Every time step is represented by one line and values of different 
/// observables are separated by spaces. Each line _must_ hold the same 
/// number of observables.
/// 
/// This class also works with files that contain no header, and consists only
/// of data points with type value_t.
///
template <typename value_t = double>
class data_file{
    private:
        typedef std::vector<value_t>    vec_1d_t; // one-dimensional vector type
        typedef std::vector<vec_1d_t>   vec_2d_t; // two-dimensional vector type

    
    public:
        /// the actual data
        vec_2d_t                    data;
        /// contains the header of the file ( ie. everything that is a comment )
        std::vector<std::string>    header;
        std::map<std::string, std::string> meta_data;
        std::vector<std::string>    observable_names;

        bool parse_header(  
            // parses the header (lines starting with comment-string at the
            // beginning of the file), and stores meta data i.e.
            // everything that was given in the format:
            // # parameter = value
            // will be stored in a hash-table meta_data, such that
            // meta_data["parameter"] == value_as_std::string
            //
            // returns false on error
            std::istream& in, 
            const std::string& comment, 
            const std::string& separator=" ");

        bool read_data_2d(
            // reads the data contained in a file with name filename
            // into vec_2d_t data, where data[0] is the first column
            // and it parses the header on the way
            const std::string& filename, 
            const std::string& separator=" ", 
            const std::string& comment="#" );

        bool read_data_2d( 
            // reads the data contained in the input stream istream
            // into vec_2d_t data, where data[0] is the first column and it
            // parses the header on the way
            std::istream& istream, 
            const std::string& separator=" ", 
            const std::string& comment="#" );

        bool read_data_2d_row_major(
            // reads the data contained in a file with name filename
            // into vec_2d_t data, where data[0] is the first row
            // and it parses the header on the way
            const std::string& filename, 
            const std::string& separator=" ", 
            const std::string& comment="#" );

        bool read_data_2d_row_major(
            // reads the data contained in the input stream istream
            // into vec_2d_t data, where data[0] is the first row and it
            // parses the header on the way
            std::istream& istream, 
            const std::string& separator=" ", 
            const std::string& comment="#");

};

template <typename value_t>
bool data_file<value_t>::parse_header( std::istream& in, const std::string& comment, const std::string& separator){
  std::streampos data_starts(0);
  std::string line("");
  while( std::getline( in, line ) ) {
    // look for comment character
    size_t pos_nonalpha = line.find_first_not_of(" \t\n");
    size_t pos = line.find(comment);
    if ( pos != line.npos && pos == pos_nonalpha ) {
      // line w/o comment
      std::string linec = string_helper::get_right_of(line, comment);
      header.push_back(linec/*.substr(pos, std::string::npos)*/);
      std::vector<std::string> split 
        = string_helper::split(linec, "=");
      if (split.size() == 2){
          meta_data[split[0]] = split[1];
      }
      // get the position of the cursor in the file
      // this points to the next line
      data_starts = in.tellg();
    } 
    else if ( pos == line.npos ){ // no comment character found
      // try deducing the observable names from the last line
      std::stringstream sstr;
      if( header.size() > 0 ){
        std::string obsnameline = header[header.size()-1];
        if ( ! obsnameline.empty() ) {
          sstr << obsnameline;
          std::string valstr("");
          while ( std::getline( sstr, valstr, separator[0] ) ) {
            if (valstr.substr(0,1) != comment){
              string_helper::trim_whitespaces(valstr);
              if (!valstr.empty())
                observable_names.push_back(valstr);
            }
          }
        }
    }
    in.seekg(data_starts);
    return true;
    }
    else { // some comment character found, but it is not the first
      // non-space-character
      std::cerr << "Error: Found comment character, but it is not"
          << " the first non-space character on the line.\n";
      return false;
    }
  }
  // if this return is conducted, the file only consists of a header
  return true;
}

template <typename value_t>
bool data_file<value_t>::read_data_2d( const std::string& filename, const std::string& separator, const std::string& comment) {
  std::ifstream infile(filename.c_str());
  if ( ! infile.is_open() ) { 
      std::cerr << "Error: " << filename << " could not be opened.\n";
      return false;
  }
  // read the whole filestream into a stringstream to fix 
  // file-ending problems under Windows
  std::stringstream sstream;
  sstream << infile.rdbuf();
  infile.close();
  bool successful = read_data_2d(sstream, separator, comment);
  return successful;
}

/// reads a data file into a two-dimensional vector
template <typename value_t>
bool data_file<value_t>::read_data_2d( std::istream& istream, const std::string& separator, const std::string& comment) {
  if ( separator.size() != 1 ) {
    std::cerr << "Error: bad separator.\n";
    return false;
  }

  // prepare members
  header.clear();
  data.clear();
  meta_data.clear();
  observable_names.clear();

  unsigned cols = 0; // counts the number of data columns
  unsigned lines = 0; // counts the number of lines of data read in
  if ( ! parse_header( istream, comment, separator ) ) {
    std::cerr << "Error: Could not parse header.\n";
    return false;
  }
  
  // the first line of data is treated separately, since
  // we want to know how many columns are present
  std::string line("");
  std::stringstream sstr;
  std::getline(istream, line);
  if ( ! line.empty() ) {
    sstr << line;
    std::string valstr("");
    while ( std::getline( sstr, valstr, separator[0] ) ) {
      value_t val;
      std::stringstream converter(valstr);
      if (converter >> val){
        data.push_back( vec_1d_t(1, val) );
        ++cols;
      }
    }
    ++lines;
  }

  if(observable_names.size() != cols){
    observable_names.clear();
    for ( unsigned col = 0; col < cols; ++col ){
      observable_names.push_back("$" + std::to_string(col));
    }
  }

  // The header is known and the number of columns 
  // was deduced from the first line of data.
  // The remaining lines of data will be parsed with more efficient
  // code that follows:

  // First, the file will be put into a buffer.
  // Therefore, we need to know how many characters are to be read.
  // 
  // Store the current position in the file
  std::streampos start_p = istream.tellg();
  // Jump to the end of the file
  istream.seekg(-1, std::ios::end);
  // Go backwards until a character is found that is not a whitespace
  while( std::isspace(istream.peek()) ) {
    istream.seekg(-1, std::ios::cur);
  }
  // calculate the number of characters that will be read in:
  // that is the current cursor position minus the starting position
  // plus one because we add the '\0' control character to the last
  // position in the buffer
  unsigned long length = istream.tellg() - start_p + 1;
  // create a buffer -- use std::vector for exception safety
  // note: you should compile with std=c++0x at least to ensure
  // that the underlying array in the vector is contiguous in 
  // memory
  //
  // reset cursor to the beginning of the file
  istream.clear();
  istream.seekg(start_p);
  // read the file into the buffer
  std::string buffer(std::istreambuf_iterator<char>(istream), {});
  if (!istream.good()) {
    std::cerr << "Error reading file into buffer\n";
  }
  // go through the buffer once for counting the lines and 
  // replace seperator character with spaces so that
  // istream_iterator and the C-api functions strto* 
  // will skip them implicitely
  //
  unsigned long nr_lines = 0;
  for (size_t p = 0; p < buffer.size(); ++p) {
    if (separator == buffer.substr(p, separator.size())) {
      for (size_t j = p; j < p + separator.size(); ++j) { buffer[j] = ' '; }
    }
    if ('\n' == buffer[p]) ++nr_lines; 
  }

  // resize each column to hold the number of lines ( plus the one
  // already read in when walking through the header )
  for ( unsigned col = 0; col < cols; ++col ){
    data.at(col).resize( nr_lines + 1 );
  }
  
  // specialize when underlying data type is some 
  // floating point number -- the C api functions strtod, strtold and 
  // strtof are the fastest way to parse input
  //  
  if ( 
  std::is_same<value_t, double>::value 
  or std::is_same<value_t, long double>::value 
  or std::is_same<value_t, float>::value){
    char* token_p = &buffer.front();
    char* end_p   = token_p + length;
    value_t val;
    while( token_p < end_p and cols > 0 ) {
      for ( unsigned col = 0; col < cols; ++col ){
        // read value_t at token_p and let it point behind it
        if( std::is_same<value_t, double>::value )
          val = strtod(token_p, &token_p);
        else if( std::is_same<value_t, long double>::value )
          val = strtold(token_p, &token_p);
        else if( std::is_same<value_t, float>::value )
          val = strtof(token_p, &token_p);
        // clean up when ill-formatted input is parsed
        if ( nr_lines < lines || errno != 0 ) {
          std::cerr << "Error: Ill-formatted file.\n";
          return false;
        }
        data[col][lines] = val;
      }
      ++lines;
    }
  }
  // use istream_iterator for other data types
  // use istream directly for all other data types
  else {
    istream.clear();
    istream.seekg(start_p);
    value_t val; 
    unsigned col = 0;
    while ( istream >> val ) {
      data[col][lines] = val;
      ++col;
      if ( col == cols ) { col = 0; ++lines; }
    }
  }
//  std::cout << "# Detected " << cols << " column(s) and " << nr_lines+1 << " row(s).\n";
//  std::cout << "# Read " << lines << " row(s).\n";
  return true;
}

template <typename value_t>        
/// reads a data file into a two-dimensional vector
bool data_file<value_t>::read_data_2d_row_major( const std::string& filename, 
    const std::string& separator, const std::string& comment) {
  std::ifstream infile(filename.c_str());
  if ( ! infile.is_open() ) { 
    std::cerr << "Error: " << filename << " could not be opened.\n";
    return false;
  }
  std::stringstream sstream;
  sstream << infile.rdbuf();
  infile.close();
  bool successful = read_data_2d_row_major(sstream, separator, comment);
  return successful;
}

template <typename value_t>
bool data_file<value_t>::read_data_2d_row_major( std::istream& istream, 
  const std::string& separator, const std::string& comment) {

  if ( separator.size() != 1 ) {
    std::cerr << "Error: bad separator.\n";
    return false;
  }

  // prepare members
  header.clear();
  data.clear();
  meta_data.clear();
  observable_names.clear();

  // helper that stores the first row
  vec_1d_t first_row;

  unsigned cols = 0; // counts the number of data columns
  unsigned lines = 0; // counts the number of lines of data read in
  if ( ! parse_header( istream, comment, separator ) ) {
    std::cerr << "Error: Could not parse header.\n";
    return false;
  }
  
  // the first line of data is treated separately, since
  // we want to know how many columns are present
  std::string line("");
  std::stringstream sstr;
  std::getline(istream, line);
  if ( ! line.empty() ) {
    sstr << line;
    std::string valstr("");
    //while ( std::getline( sstr, valstr, separator[0] ) ) {
    while ( std::getline( sstr, valstr, separator[0] ) ) {
      value_t val;
      std::stringstream converter(valstr);
      if (converter >> val){
        first_row.push_back(val);
        ++cols;
      }
    }
    data.push_back( first_row );
    ++lines;
  }

  if(observable_names.size() != cols){
    observable_names.clear();
    for ( unsigned col = 0; col < cols; ++col ){
      observable_names.push_back("$" + std::to_string(col));
    }
  }

  // The header is known and the number of columns 
  // was deduced from the first line of data.
  // The remaining lines of data will be parsed with more efficient
  // code that follows:

  // First, the file will be put into a buffer.
  // Therefore, we need to know how many characters are to be read.
  // 
  // Store the current position in the file
  std::streampos start_p = istream.tellg();
  // Jump to the end of the file
  istream.seekg(-1, std::ios::end);
  // Go backwards until a character is found that is not a whitespace
  while( std::isspace(istream.peek()) ) {
    istream.seekg(-1, std::ios::cur);
  }
  // calculate the number of characters that will be read in:
  // that is the current cursor position minus the starting position
  // plus one because we add the '\0' control character to the last
  // position in the buffer
  std::streamoff length = istream.tellg() - start_p + 1;
  // create a buffer -- use std::vector for exception safety
  // note: you should compile with std=c++0x at least to ensure
  // that the underlying array in the vector is contiguous in 
  // memory
  //
  // reset cursor to the beginning of the file
  istream.clear();
  istream.seekg(start_p);
  // read the file into the buffer
  std::string buffer(std::istreambuf_iterator<char>(istream), {});
  if (!istream.good()) {
    std::cerr << "Error reading file into buffer\n";
  }
  // go through the buffer once for counting the lines and 
  // replace seperator character with spaces so that
  // istream_iterator and the C-api functions strto* 
  // will skip them implicitely
  //
  unsigned long nr_lines = 0;
  for (size_t p = 0; p < buffer.size(); ++p) {
    if (separator == buffer.substr(p, separator.size())) {
      for (size_t j = p; j < p + separator.size(); ++j) { buffer[j] = ' '; }
    }
    if ('\n' == buffer[p]) ++nr_lines; 
  }

  // resize data to hold the number of lines ( plus the one
  // already read in when walking through the header )
  data.resize( nr_lines + 1, vec_1d_t(cols, 0) );
  
  // specialize when underlying data type is some 
  // floating point number -- the C api functions strtod, strtold and 
  // strtof are the fastest way to parse input
  //  
  if ( 
  std::is_same<value_t, double>::value 
  or std::is_same<value_t, long double>::value 
  or std::is_same<value_t, float>::value){
    char* token_p = &buffer.front();
    char* end_p   = token_p + length;
    value_t val;
    while( token_p < end_p and cols > 0 ) {
      for ( unsigned col = 0; col < cols; ++col ){
        // read value_t at token_p and let it point behind it
        if( std::is_same<value_t, double>::value )
          val = strtod(token_p, &token_p);
        else if( std::is_same<value_t, long double>::value )
          val = strtold(token_p, &token_p);
        else if( std::is_same<value_t, float>::value )
          val = strtof(token_p, &token_p);
        // clean up when ill-formatted input is parsed
        if ( nr_lines < lines || errno != 0 ) {
          std::cerr << "Error: Ill-formatted file.\n";
          return false;
        }
        data[lines][col] = val;
      }
      ++lines;
    }
  }
  // use istream directly for all other data types
  else {
    istream.clear();
    istream.seekg(start_p);
    value_t val; 
    unsigned col = 0;
    while ( istream >> val ) {
      data[lines][col] = val;
      ++col;
      if ( col == cols ) { col = 0; ++lines; }
    }
  }
//  std::cout << "# Detected " << cols << " column(s) and " << nr_lines+1 << " row(s).\n";
//  std::cout << "# Read " << lines << " row(s).\n";
  return true;
}

# endif // READFILE_HPP
