#ifndef RANDOM_NUMBER_GENERATOR_H
#define RANDOM_NUMBER_GENERATOR_H

#include <cstdio>
#include <cmath>


using namespace std;

typedef unsigned long uint32;

#define MT_BUFSZ       (624)     // length of state vector
// RANGE_LT1 puts a uint32 in the range [0,1) (0<=x<1)
#define RANGE_LT1 2.3283064370807974e-10


// This is the ``Mersenne Twister'' random number generator MT19937, which
// generates pseudorandom integers uniformly distributed in 0..(2^32 - 1)
// starting from any odd seed in 0..(2^32 - 1).  This version is a recode
// by Shawn Cokus (Cokus@math.washington.edu) on March 8, 1998 of a version by
// Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
// July-August 1997).
class MT_RNG 	
{
	private:
    	uint32   state[MT_BUFSZ+1];  // state vector + 1 extra to not violate ANSI C
    	uint32   *next;          	// next random value is computed from here
    	int      left;      		// can *next++ this many times before reloading

    	uint32 reloadMT();

	public:
		MT_RNG(const uint32 seed = 615460891);
    	void Seed(const uint32 seed);
    	
		uint32 lRand();
    	inline double dRand() { // double in  [0,1)-interval
			return ((double)lRand() * RANGE_LT1); 	
		}
		inline double operator()() {
			return dRand();
		}
};

class nonrep_mt_rng {
	private:
		uint32* seed_cache;
		double* cache;
		int seed_cache_fill;
		int cache_fill;
		MT_RNG rng;
		FILE* seedsource;
		
		void fill_cache() {
			if (seed_cache_fill == 0) {
				if (seedsource) {
					fread(reinterpret_cast<void*>(seed_cache), sizeof(uint32), 1024, seedsource);
				} else {
					for (int i=0; i<1024;i++)
						seed_cache[i] += uint32(i);
				}
				seed_cache_fill = 1024;
			}
			
			rng.Seed(seed_cache[--seed_cache_fill]);

			for (int i=0; i<4096;i++)
				cache[i] = rng();
				
			cache_fill = 4096;
		}
		
	public:
		nonrep_mt_rng() 
			: seed_cache(0), cache(0), seed_cache_fill(0), cache_fill(0),
				seedsource(0) 
		{
			seedsource = fopen("/dev/urandom", "r");
			seed_cache = new uint32[1024];
			cache = new double[4096];
		}
		~nonrep_mt_rng() {
			delete[] cache;
			delete[] seed_cache;
		}
		inline double operator()() {
			if (!cache_fill)
				fill_cache();
				
			return cache[--cache_fill];
		}
};


template<class RNG>
class Random {
	protected:
		
	
		template<class T> 
		void swap_randomly(T& v, const long pos) {
 			const long index = randindex(v.size()-pos);
 			typename T::value_type value = v[pos+index];
 			v[pos+index] = v[pos];		// swap elements 
 			v[pos] = value;
		} 
	
	public:
		static RNG rng;
		unsigned long randindex(const unsigned long N) {	// return random index in range [0 ... N)
			double x = rng();
			return  (unsigned long) floor(x * ((double)N)); 
		} 		
		bool get_bool() {
			return (rng() > 0.5);
		}
		double operator()() {
			return rng();
		}
		template<class T>
		T operator()(const T lower_bound, const T upper_bound) {
			return lower_bound + ((T)(rng() * (upper_bound - lower_bound)));
		}
		template<class T>		// randomly permute elements of in container of type T
		void permute(T& sequence) {
			for (long i=0; i<sequence.size(); i++) 
				swap_randomly(sequence, i);	
		}
};
	
template<class RNG>
RNG Random<RNG>::rng;
	
	
#undef RANGE_LT1 
#undef MT_BUFSZ

#endif
