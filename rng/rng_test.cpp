// rng_test.cpp

#include "../utilities/tx_tests.h"
#include "rng.hpp"


using namespace tx_rng;


#define ASSERT TX_ASSERT
#define ASSERT_MSG TX_ASSERT_MSG
#define ASSERT_MSG_VA TX_ASSERT_MSG

#define APPROX TX_APPROX
#define APPROX_ARRAY TX_APPROX_ARRAY


template <typename T>
T array_mean( T *array, size_t length )
{
    T sum = 0.0;
    for ( size_t i=1; i<length; ++i ) {
        sum += array[i];
    }
    return sum / (T) length;
}


template <typename T>
f64 array_variance( T *array, size_t length )
{
    T mean = array_mean( array, length );
    T tmp, sum  = 0.0;
    for ( size_t i=1; i<length; ++i ) {
        tmp = ( array[i] - mean );
        sum += tmp * tmp;
    }
    return sum / (T) (length - 1);
}


// ------------------------------

void test_xorshift1024s( void )
{
    Xorshift1024s x;
    xorshift1024s_seed( &x, 37473 );

    ASSERT( xorshift1024s_next( &x ) > 0 );

    xorshift1024s_jump( &x);

    ASSERT( xorshift1024s_next( &x ) > 0 );
}


void test_xoshiro256ss( void )
{
    Xoshiro256ss x;
    xoshiro256ss_seed( &x, 37473 );

    ASSERT( xoshiro256ss_next( &x ) > 0 );

    xoshiro256ss_jump( &x);

    ASSERT( xoshiro256ss_next( &x ) > 0 );
}


void test_xoshiro128p( void )
{
    Xoshiro128p x;
    xoshiro128p_seed( &x, 37473 );

    ASSERT( xoshiro128p_next( &x ) > 0 );

    xoshiro128p_jump( &x);

    ASSERT( xoshiro128p_next( &x ) > 0 );
}


void test_xoshiro128pv( void )
{
    Xoshiro128pv x;
    xoshiro128pv_seed( &x, 37473 );

    ASSERT( xoshiro128pv_next( &x ) > 0 );

    u32 array[255];
    xoshiro128pv_nextn( &x, array, 255 );

    f64 mean = 0.0;
    for ( size_t i = 0; i < 255; ++i ) {
        mean += to_f64(array[i]);
        ASSERT( array[i] > 0 );
    }
    mean /= 255;
    printf("%f\n", mean);

    xoshiro128pv_jump( &x);

    ASSERT( xoshiro128pv_next( &x ) > 0 );
}


void test_generator( void )
{
    u64 seed = seed_value();

    Xoshiro256ss state;
    Generator rng = xoshiro256ss_init( &state, seed );

    jump( rng );

    ASSERT( next( rng ) > 0 );
}


void test_uniform( void )
{
    Xoshiro256ss state;
    Generator rng = xoshiro256ss_init( &state, 3747 );

    ASSERT( uniform( rng ) <= 1.0 );

    f64 buffer[256];
    uniform_n( rng, buffer, 256 );

    for ( u32 i = 0; i < 256; ++i ) {
        ASSERT( 0.0 <= buffer[i] && buffer[i] <= 1.0 );
    }

    ASSERT( uniform_pdf( 1.0, 3.0 ) == 0.5 );
    ASSERT( uniform_lpdf( 1.0, 3.0 ) == std::log(0.5) );
    ASSERT( uniform_cdf( 2.0, 1.0, 3.0 ) == 0.5 );
}


void test_normal( void )
{
#define N 1000

    f64 array[N];

    Xorshift1024s state;
    Generator rng = xorshift1024s_init( &state, 3747 );

    for ( u32 i=0; i<N; ++i ) {
        array[i] = normal( rng );
    }

    f64 mean = array_mean( array, N );
    f64 var  = array_variance( array, N );

    ASSERT( APPROX( mean, 0,   0.01 ) );
    ASSERT( APPROX( var,  1.0, 0.1 ) );

    f64 pdf  = normal_pdf( 1.0, 1.0, 1.0 );
    f64 lpdf = normal_lpdf( 1.0, 1.0, 1.0 );
    f64 cdf  = normal_cdf( 1.0, 1.0, 1.0 );

    ASSERT( APPROX( pdf,  0.3989423,  1E-7 ) );
    ASSERT( APPROX( lpdf, std::log( pdf ), 1E-14 ) );
    ASSERT( APPROX( cdf,  0.5,        1E-14 ) );
#undef N
}


TX_TEST_LIST = {
    TX_ADD_TEST(test_xorshift1024s),
    TX_ADD_TEST(test_xoshiro256ss),
    TX_ADD_TEST(test_xoshiro128p),
    TX_ADD_TEST(test_xoshiro128pv),
    TX_ADD_TEST(test_generator),
    TX_ADD_TEST(test_uniform),
    TX_ADD_TEST(test_normal),
};


TX_TESTS_MAIN

