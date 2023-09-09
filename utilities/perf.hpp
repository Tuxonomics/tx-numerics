// perf.hpp

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
	uint64_t oldstate = rng->state;
	// Advance internal state
	rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
	// Calculate output function (XSH RR), uses old state for max ILP
	uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
	uint32_t rot = oldstate >> 59u;
	return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// + more

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

bool pcg32_random_r_probability(pcg32_random_t* rng, float probablity)
{
	uint64_t oldstate = rng->state;
	// Advance internal state
	rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
	// Calculate output function (XSH RR), uses old state for max ILP
	uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
	uint32_t rot = oldstate >> 59u;
	uint32_t rval = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
	return (rval & 16777215) < (probablity*16777215);
}

static uint32_t big_int_val = 1<<30;
static uint32_t big_int_mask = big_int_val-1;
uint32_t pcg32_random_r_range(pcg32_random_t* rng, uint32_t minVal, uint32_t maxVal)
{
	uint64_t oldstate = rng->state;
	// Advance internal state
	rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
	// Calculate output function (XSH RR), uses old state for max ILP
	uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
	uint32_t rot = oldstate >> 59u;
	uint32_t rval = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
	float inRange = ((float)(rval&big_int_mask)) / big_int_val;
	return minVal + (maxVal-minVal) * inRange;
}

float pcg32_random_r_rangef(pcg32_random_t* rng, float minVal, float maxVal)
{
	uint64_t oldstate = rng->state;
	// Advance internal state
	rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
	// Calculate output function (XSH RR), uses old state for max ILP
	uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
	uint32_t rot = oldstate >> 59u;
	uint32_t rval = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
	float inRange = ((float)(rval&big_int_mask)) / big_int_val;
	return minVal + (maxVal-minVal) * inRange;
}

// utilities based on
// https://github.com/raspofabs/dodbooksourcecode/blob/master/common.h

#include <string.h>
#include <algorithm>
#include <cmath>
#include <chrono>


const float TRIAL_TIMEOUT = 20.0f * 10000.0f;

class Timer {
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<ms_>
            (clock_::now() - beg_).count();
    }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::milli> ms_;
    std::chrono::time_point<clock_> beg_;
};


struct Stat {
	float average;
	float standardDeviation;
	float ninetyFivePercentMin;
	float ninetyFivePercentMax;
	float lowest, highest;
	void Calculate( float *series, int length ) {
		lowest = highest = series[0];
		float sum = 0;
		for( int i = 0; i < length; ++i ) {
			float val = series[i];
			lowest = lowest < val ? lowest : val;
			highest = highest > val ? highest : val;
			sum += val;
		}
		average = sum / length;
		float sds = 0;
		for( int i = 0; i < length; ++i ) {
			float d = series[i] - average;
			sds += d*d;
		}
		standardDeviation = sqrtf( sds / (length-1) );
		ninetyFivePercentMin = average - standardDeviation*2;
		ninetyFivePercentMax = average + standardDeviation*2;
	}
	void Magnify( float multiplier ) {
		average *= multiplier;
		standardDeviation *= multiplier;
		ninetyFivePercentMin *= multiplier;
		ninetyFivePercentMax *= multiplier;
		lowest *= multiplier;
		highest *= multiplier;
	}
};


volatile int touchVal;

struct CacheClearer {
	// a massive array used to read and write to, to flush any data in the cache.
	char *horribleBuffer = 0;
	static const int HORRIBLE_BUFFER_SIZE = 12 * 1024 * 1024;
	CacheClearer() {
		horribleBuffer = (char*)malloc( HORRIBLE_BUFFER_SIZE );
		if( horribleBuffer ) {
			memset( horribleBuffer, touchVal, HORRIBLE_BUFFER_SIZE );
		}
	}
	~CacheClearer() {
		if( horribleBuffer ) {
			touchVal = horribleBuffer[(touchVal*touchVal)%HORRIBLE_BUFFER_SIZE] + 1;
			free( horribleBuffer );
		}
		horribleBuffer = 0;
	}
	void ClearCaches() {
		for( int i = 0; i < HORRIBLE_BUFFER_SIZE; ++i ) {
			horribleBuffer[i] += 1;
		}
	}
};


typedef void (*TestFunc)();
struct Test {
	TestFunc func;
	char name[64];
	Test( TestFunc f, const char* n ) : func(f) { strcpy( name, n ); }
};

void Shuffle( pcg32_random_t *pcg, int *array, int length ) {
	for( int i = 0; i < length-1; ++i ) {
		int source = pcg32_random_r_range( pcg, i, length );
		int t = array[i];
		array[i] = array[source];
		array[source] = t;
	}
#if 0
	printf( "New order = " );
	for( int i = 0; i < length; ++i ) {
		printf( "%i ", array[i] );
	}
	printf( "\n" );
#endif
}

template<typename TestStruct, int numTests>
void RunTests( TestStruct (&testArray)[numTests] ) {
	CacheClearer cc;

	const int TRIALS = 128;

	struct TimingData {
		int testID;
		float timing[TRIALS];
		Stat s;
	};
	TimingData timingData[numTests];
	int testToDo[numTests];
	for( int i = 0; i < numTests; ++i ) {
		testToDo[i] = i;
	}
	pcg32_random_t pcg;
	pcg32_srandom_r( &pcg, 1234, 5678 );
	Shuffle( &pcg, testToDo, numTests );
	Timer trialTimer;
	int trial = 0;
	while( trial < TRIALS && trialTimer.elapsed() < TRIAL_TIMEOUT ) {
		Shuffle( &pcg, testToDo, numTests );
		for( int selector = 0; selector < numTests; ++selector ) {
			int tid = testToDo[selector];

			cc.ClearCaches();
			auto &test = testArray[tid];
			auto &td = timingData[tid];
			td.testID = tid;
			// warm the engines
			Timer t;
			test.func();
			td.timing[trial] = t.elapsed();
		}
		++trial;
	}
	for( auto &td : timingData ) {
		td.s.Calculate( td.timing, trial );
	}
	printf( "Managed %i trials in %fms\n\n", trial, trialTimer.elapsed() );

	std::sort(
			std::begin(timingData),
			std::end(timingData),
			[]( const TimingData &a, const TimingData &b ){
				return a.s.average > b.s.average;
			} );

	const char *timesuffix = "ms";
	if( timingData[0].s.average < 20 ) {
		timesuffix = "us";
		for( auto &td : timingData ) {
			td.s.Magnify( 1000.0f );
		}
	}

	for( int tid = 0; tid < numTests; ++tid ) {
		auto &td = timingData[tid];
		auto &test = testArray[td.testID];
		printf( "Average% 9.2f%s (sd %.2f%s) for test [%s]\n",
				td.s.average, timesuffix,
				td.s.standardDeviation, timesuffix,
				test.name );
	}
#ifdef TARGET
#define TOSTR(X) #X
#define STR(x) TOSTR(x)
	char outfilename[256];
	sprintf( outfilename, "testreport_%s.txt", STR(TARGET) );
	if( FILE *fp = fopen( outfilename, "wt" ) ) {
		for( int tid = 0; tid < numTests; ++tid ) {
			auto &td = timingData[tid];
			auto &test = testArray[td.testID];
			fprintf( fp, "Average% 9.2f%s (sd %.2f%s) for test [%s]\n",
					td.s.average, timesuffix,
					td.s.standardDeviation, timesuffix,
					test.name );
		}
		fclose( fp );
		//printf( "Wrote report for " STR(TARGET) " to file %s\n", outfilename );
	} else {
		printf( "unable to open file %s for writing\n", outfilename );
	}
	sprintf( outfilename, "testdata_%s.csv", STR(TARGET) );
	if( FILE *fp = fopen( outfilename, "wt" ) ) {
		fprintf( fp, "trialseq,testtime,testname\n" );
		for( int tid = 0; tid < numTests; ++tid ) {
			auto &td = timingData[tid];
			auto &test = testArray[td.testID];
			for( int i = 0; i < trial; ++i ) {
				fprintf( fp, "%i,%f,%s\n", i, td.timing[i], test.name );
			}
		}
		fclose( fp );
		//printf( "Wrote data for " STR(TARGET) " to file %s\n", outfilename );
	} else {
		printf( "unable to open file %s for writing\n", outfilename );
	}
#endif
}