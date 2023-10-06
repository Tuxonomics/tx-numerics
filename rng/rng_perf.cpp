// rng_perf.cpp

#include "rng.hpp"
#include "../utilities/perf.hpp"


using namespace std;
using namespace tx_rng;

#define BUFFER_SIZE (4096 << 10)
f64 buffer[BUFFER_SIZE];


Generator rng;

void init( void )
{
    Xoshiro256ss *state = (Xoshiro256ss *) malloc( sizeof(*state) );
    rng = xoshiro256ss_init( state, 3747 );
}


void test_uniform_batches( void )
{
    uniform_n( rng, buffer, BUFFER_SIZE );
}


void test_uniform_scalar( void )
{
    for ( u32 i = 0; i < BUFFER_SIZE; ++i ) {
        buffer[i] = uniform( rng );
    }
}



int main( int argn, const char *argv[] )
{
    init();

    Test tests[] = {
        (Test){ test_uniform_batches, "test_uniform_batches" },
        (Test){ test_uniform_scalar, "test_uniform_scalar" },
    };

    RunTests( tests );

    return 0;
}

