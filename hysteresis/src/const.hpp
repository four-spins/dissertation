/**
* @brief header defining types and constants
*
* @author Marco Mueller
* @date 11.04.2011
*
*/
# ifndef __CONST_HPP__
# define __CONST_HPP__

# include <vector>
# include "dlib/uintn.h"
# include "dlib/string.h"
# include <map>
# include <fstream>
// logfile for output
std::ofstream logger_ofstream;

struct double_one
{
    double value;
    double_one(): value(1.0){} // default value 100
};

// dlib types that are checked at compile time
using dlib::uint32;
using dlib::int32;
using dlib::uint8;
using dlib::uint64;

// typedefs for readability
typedef std::vector<uint32> uint32_1d; 
    ///< one dimensional vector of unsigned ints
typedef std::vector<uint32_1d> uint32_2d;
    ///< two dimensional vector of unsigned ints (e.g. for neighbour table)
typedef std::vector<uint32_2d> uint32_3d;
    ///< 3-dimensional vector of unsigned ints (e.g for plaquette table) 

typedef std::vector<int32> int32_1d;
    ///< one dimensional vector of ints
typedef std::vector<int32_1d> int32_2d;
    ///< two dimensional vector of ints

// constants for boundary conditions (for readability, consistency)
const uint8 PERIODIC = 0;
const uint8 INTERFACE = 1;
const uint8 VANISH = 2;
const uint8 HELICAL = 3;

template <class T> inline
std::string join(T& vec)
{
    std::string temp = "";
    for(size_t i = 0; i < vec.size() - 1; i++)
        temp += dlib::cast_to_string(vec.at(i)) + ", ";
    temp += dlib::cast_to_string(vec.at(vec.size() - 1));
    return temp;

}

template <class T> inline
std::string joinm(std::map<double, T>& map)
{
    std::string temp = "";
    typename std::map<double, T>::iterator it;
    it = map.begin();
    while(it != map.end())
    {   
        T foo = it->second;
        temp += dlib::cast_to_string(foo.value()) + ", ";
        it ++;
    }
    return temp;
}

inline std::string joinm(std::map<double, uint64>& map)
{
    std::string temp = "";
    std::map<double, uint64>::const_iterator it;
    it = map.begin();
    while(it != map.end())
    {   
        uint64 foo = it->second;
        temp += dlib::cast_to_string(foo) + ", ";
        it ++;
    }
    return temp;
}

# endif
