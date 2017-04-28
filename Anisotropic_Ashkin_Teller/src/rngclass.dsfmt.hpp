#ifndef __RNGCLASS_DSFMT_HPP__
#define __RNGCLASS_DSFMT_HPP__

//#include "dlib/all/source.cpp"
//#define NO_MAKEFILE

// random number generator
#include "dSFMT.h"

class rngclass
{

	private: 
		// random number generator initialization
        unsigned seed;	
        dsfmt_t rng;

    public:
	
	rngclass()
    : seed(static_cast<int>(time(NULL)))
	{
		// calculate some seed out of the system's time, 
		// if there is no seed the environment variables given
		// seed the generator
        dsfmt_init_gen_rand(&rng, seed);
	}
	
	rngclass(int seed_)
    : seed(seed_)
	{
        dsfmt_init_gen_rand(&rng, seed);
	}

    void set_seed(int seed_)
    {
		// calculate some seed out of the system's time, 
		// seed the generator
        seed = seed_;
        dsfmt_init_gen_rand(&rng, seed);
    }

	~rngclass()
	{
	}
	
	/// draws a random number out of (0,1)
	inline double rand()
	{
        return dsfmt_genrand_close_open(&rng);
	}

	const char * get_type()
	{
		return "dsfmt";
	}

    unsigned get_seed()
    {
        return seed;
    }
	
};
#endif //__RNGCLASS_DSFMT_HPP__
