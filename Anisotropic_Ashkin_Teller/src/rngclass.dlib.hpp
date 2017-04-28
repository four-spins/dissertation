#ifndef __RNGCLASS_DLIB_HPP__
#define __RNGCLASS_DLIB_HPP__

//#include "dlib/all/source.cpp"
//#define NO_MAKEFILE

// random number generator
#include "dlib/rand.h"
typedef dlib::rand::float_1a random_number_generator;

class rngclass
{

	private: 
		// random number generator initialization
        unsigned seed;	
        random_number_generator rng;

    public:
	
	rngclass()
    : seed(static_cast<int>(time(NULL)))
	{
		// calculate some seed out of the system's time, 
		// if there is no seed the environment variables given
		// seed the generator
        rng.set_seed(dlib::cast_to_string<int>(seed));
	}
	
	rngclass(int seed_)
    : seed(seed_)
	{
        rng.set_seed(dlib::cast_to_string<int>(seed));
	}

    void set_seed(int seed_)
    {
		// calculate some seed out of the system's time, 
		// seed the generator
        seed = seed_;
        rng.set_seed(dlib::cast_to_string<int>(seed));
    }

	~rngclass()
	{
	}
	
	/// draws a random number out of (0,1)
	inline double rand()
	{
		return rng.get_random_double();
	}

	const char * get_type()
	{
		return "dlib mersenne twister";
	}

    unsigned get_seed()
    {
        return seed;
    }
	
};
#endif //__RNGCLASS_DLIB_HPP__
