/**
* @brief class wrapper for random numbers by gsl
*
* @author Marco Mueller
* @date 17.08.2010
* @note guest student programme on parallel computing at JSC
*
*/
# pragma once
// random number generator using gsl
#include "gsl/gsl_rng.h"
#include <ctime>

class rngclass
{

	private: 
		gsl_rng * rng;  /* global generator */
		// random number generator initialization
		const gsl_rng_type * T;
        unsigned seed;	
	
    public:
	
	rngclass()
	{
		// get environment variables
		gsl_rng_env_setup();
		T = gsl_rng_default;
		// allocate the actual random number generator
		rng = gsl_rng_alloc (T);
		// calculate some seed out of the system's time, 
		// if there is no seed the environment variables given
		seed = gsl_rng_default_seed ? gsl_rng_default_seed : static_cast<unsigned>(time(NULL));
		// seed the generator
		gsl_rng_set(rng,seed);
	}

	rngclass(unsigned int seed_)
    : seed(seed_)
	{
		// get environment variables
		gsl_rng_env_setup();
		T = gsl_rng_default;
		// allocate the actual random number generator
		rng = gsl_rng_alloc (T);
		// calculate some seed out of the system's time, 
		// if there is no seed the environment variables given
		seed = gsl_rng_default_seed ? gsl_rng_default_seed : seed;
		// seed the generator
		gsl_rng_set(rng,seed);
	}

	~rngclass()
	{
		gsl_rng_free (rng);
	}
	
	/// draws a random number out of (0,1)
	double rand()
	{
		return gsl_rng_uniform_pos(rng);
	}

	const char * get_type()
	{
		return gsl_rng_name(rng);
	}

    unsigned get_seed()
    {
        return seed;
    }

};
