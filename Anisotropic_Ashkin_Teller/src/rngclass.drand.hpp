/**
* @brief class wrapper for random numbers using Marsaglia RNG
*
* @author Marco Mueller
* @date 26.08.2010
* @note guest student programme on parallel computing at JSC
*
*/
#ifndef __RNGCLASS_DRAND_HPP__
#define __RNGCLASS_DRAND_HPP__

# include<cstdlib>

class rngclass
{

	private: 
		// random number generator initialization
        unsigned seed;	

    public:
	
	rngclass()
    : seed(static_cast<int>(time(NULL)))
	{
		// calculate some seed out of the system's time, 
		// if there is no seed the environment variables given
		// seed the generator
		srand48(seed);
	}
	
	rngclass(int seed_)
    : seed(seed_)
	{
        srand48(seed);
	}

    void set_seed(int seed_)
    {
		// calculate some seed out of the system's time, 
		// seed the generator
        seed = seed_;
        srand48(seed);
    }

	~rngclass()
	{
	}
	
	/// draws a random number out of (0,1)
	inline double rand()
	{
		return drand48();
	}

	const char * get_type()
	{
		return "DRAND48 RNG";
	}

    unsigned get_seed()
    {
        return seed;
    }
	
};
#endif //__RNGCLASS_DRAND_HPP__
