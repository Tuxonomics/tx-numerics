// rng_perf.cpp

#include "rng.hpp"
#include "../utilities/perf.hpp"


using namespace std;
using namespace tx_rng;

// #define BUFFER_SIZE (4096 << 10)
#define BUFFER_SIZE (4096)
f64 buffer[BUFFER_SIZE];

Generator<u32> rng1;
Generator<u32> rng2;

void init( void )
{
    // Xoshiro256ss *state1 = (Xoshiro256ss *) malloc( sizeof(*state1) );
    // rng1 = xoshiro256ss_init( state1, 3747 );
    Xoshiro128p *state1 = (Xoshiro128p *) malloc( sizeof(*state1) );
    rng1 = xoshiro128p_init( state1, 3747 );

    Xoshiro128pv *state2 = (Xoshiro128pv *) malloc( sizeof(*state2) );
    rng2 = xoshiro128pv_init( state2, 3747 );
}


void test_uniform_batches_1( void )
{
    uniform_n( rng1, buffer, BUFFER_SIZE );
}


void test_uniform_scalar_1( void )
{
    for ( u32 i = 0; i < BUFFER_SIZE; ++i ) {
        buffer[i] = uniform( rng1 );
    }
}


void test_uniform_batches_2( void )
{
    uniform_n( rng2, buffer, BUFFER_SIZE );
}


void test_uniform_scalar_2( void )
{
    for ( u32 i = 0; i < BUFFER_SIZE; ++i ) {
        buffer[i] = uniform( rng2 );
    }
}



int main( int argn, const char *argv[] )
{
    init();

    Test tests[] = {
        (Test){ test_uniform_batches_1, "test_uniform_batches_1" },
        (Test){ test_uniform_scalar_1, "test_uniform_scalar_1" },

        (Test){ test_uniform_batches_2, "test_uniform_batches_2" },
        (Test){ test_uniform_scalar_2, "test_uniform_scalar_2" },
    };

    RunTests( tests );

    f64 sum = 0.0;
    for ( size_t i = 0; i < BUFFER_SIZE; ++i ) {
        sum += buffer[i];
    }

    return (int)sum;
}

