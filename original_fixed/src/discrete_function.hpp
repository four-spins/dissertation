/******************************************************************************
 *
 * @file: discrete_function.hpp
 *
 * @date: 12/06/2012 01:39:53 PM (CET)
 *
 * @author: Marco Müller <muelma@gmail.com>
 *
 ******************************************************************************/
#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip> // stream manipulators
#include <sstream> // string streams used for output
#include <limits>
#include <cmath>

/// class to manage a functional dependency y = f(x) with automatic 
/// discretisation in x
///
/// The naming of the class is as follows:\n
/// x is called the abscissa\n
/// y is called the ordinate\n
/// The interval [abscissa_min, abscissa_max] is divided in a number of
/// nr_bins subintervals of a given length where each subinterval 
/// holds exactly one value for the ordinate.\n
/// Subintervals are also referred to as "bins".
template<class T>
class discrete_function {

  public:
    typedef typename std::vector<T>::iterator               iterator;
    typedef typename std::vector<T>::reverse_iterator       reverse_iterator;
    typedef typename std::vector<T>::const_iterator         const_iterator;
    typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
    typedef T                                               value_type;

    /// constructor defining
    /// interval limits and length of one subinterval
    discrete_function(double abscissa_min_, double abscissa_max_, double bin_width_);

    // members defining the discretisation:
    double abscissa_min; ///< minimum of the interval
    double abscissa_max; ///< maximum of the interval
    double bin_width;    ///< length of a subinterval
    size_t nr_bins;
    size_t normalization;

    // vector storing the ordinates
    std::vector<T> data;

    /// access by index 
    T& operator[]( size_t index );
    const T& operator[]( size_t index ) const;
    /// access by value for the abscissa
    T& operator[]( double abscissa_val );
    const T& operator[]( double abscissa_val ) const;

    /// returns iterator to begin of underlying data
    iterator begin()               { return data.begin(); }
    /// returns const iterator to begin of underlying data
    const_iterator begin()  const  { return data.begin(); }
    /// returns const iterator to begin of underlying data
    const_iterator cbegin() const  { return data.cbegin(); }
    /// returns iterator to end of underlying data
    iterator end()                 { return data.end(); }
    /// returns const iterator to end of underlying data
    const_iterator end()  const    { return data.end(); }
    /// returns const iterator to end of underlying data
    const_iterator cend() const    { return data.cend(); }
    /// returns reverse iterator to begin of underlying data
    reverse_iterator rbegin()              { return data.rbegin(); }
    /// returns const reverse iterator to begin of underlying data
    const_reverse_iterator rbegin()  const { return data.rbegin(); }
    /// returns const reverse iterator to begin of underlying data
    const_reverse_iterator crbegin() const { return data.crbegin(); }
    /// returns reverse iterator to end of underlying data
    reverse_iterator rend()                { return data.rend(); }
    /// returns const reverse iterator to end of underlying data
    const_reverse_iterator rend()  const   { return data.rend(); }
    /// returns const reverse iterator to end of underlying data
    const_reverse_iterator crend() const   { return data.crend(); }

    /// get internal vector storing ordinates
    std::vector<T>& get_data()      { return data; }
    /// get lower limit of the interval
    double get_abscissa_min() const { return abscissa_min; }
    /// get upper limit of the interval
    double get_abscissa_max() const { return abscissa_max; }
    /// get size of a subinterval
    double get_bin_width() const    { return bin_width; }
    /// return normalization
    size_t get_norm() const         { return normalization; }

    /// get number of subintervals
    size_t size() const;
    /// get index for a value of the abscissa
    size_t get_index( double abscissa_val ) const;
    /// get the abscissa for an index interpolated to be the 
    /// mean value of the subinterval
    double get_abscissa(size_t index) const;
    /// get a string representation
    std::string to_string() const;

    //template<class T2>
    //friend std::ostream& operator<<(std::ostream& os, discrete_function<T2>& func);

}; // class discrete_function

// constructor definition
template<class T>
discrete_function<T>
::discrete_function(double abscissa_min_, double abscissa_max_, double bin_width_)
: abscissa_min(abscissa_min_)
  , abscissa_max(abscissa_max_)
  , bin_width(bin_width_) 
  , normalization( std::numeric_limits<double>::quiet_NaN() ){
  const double diff = abscissa_max - abscissa_min;

  // Error:
  // either the lower bound of the interval is greater than
  // the upper bound or they are equal
  assert( diff > 0 && "Error: creating discrete_function failed due to improper definition of the interval.");
  nr_bins = static_cast<size_t>( (diff / bin_width) );

  // Error:
  // length of a subinterval is bigger than the length of the 
  // whole interval
  assert( nr_bins > 0 && "Error: discrete_function needs at least one subinterval.");
  data.resize(nr_bins);
}

// number of subintervals
template<class T>
size_t discrete_function<T>
::size() const { return nr_bins; }

// access by index
template<class T>
const T& discrete_function<T>
::operator[](size_t index) const {
  assert( index < nr_bins );
  return data[index];
}

template<class T>
T& discrete_function<T>
::operator[](size_t index) {
  assert( index < nr_bins );
  return data[index];
}

// access by value
template<class T>
inline const T& discrete_function<T>
::operator[](double abscissa) const {
  return data[get_index(abscissa)];
}

template<class T>
inline T& discrete_function<T>
::operator[](double abscissa) {
  return data[get_index(abscissa)];
}

template<class T>
inline size_t discrete_function<T>
::get_index( double abscissa_val ) const {
  size_t index = 0;
  index = static_cast<size_t>((abscissa_val - abscissa_min)/bin_width);
  assert( index < nr_bins );
  return index;
}

template<class T> 
double discrete_function<T>
::get_abscissa(size_t index) const{
  double tmp =static_cast<double>(abscissa_min) + static_cast<double>(bin_width) * (index + 0.5);
  return tmp;
}

template<class T>
std::string discrete_function<T>
::to_string() const {
  std::ostringstream os;
  os  << "# discrete function:\n"
    << "# abscissa_min = " << std::scientific << std::setprecision(16) << abscissa_min << "\n"
    << "# abscissa_max = " << std::scientific << std::setprecision(16) << abscissa_max << "\n"
    << "# nr_bins = " << nr_bins << "\n"
    << "# -----\n" << "# Every line represents a right-open subinterval of the abscissa with its ordinate\n"
    << "#" << std::setw(29) << "x_mean" << std::setw(31) << "f[x_min, x_max)\n";

  // TODO: how to write type information ... logval vs. double ...
  os.setf(std::ios::floatfield);
  os.precision(16);

  for ( size_t i = 0; i < nr_bins; i++ ) {
    double mean = this->get_abscissa(i) ;
    os << std::setw(30) << mean << " "
      << std::setw(29) << this->operator[](i) << " "
      << "\n";
  }
  return os.str();
}

template<class T>
std::ostream& operator<<(std::ostream& os, discrete_function<T>& func){
  os << func.to_string();
  return os;
}
