// rng_test.cpp

#include "lin_alg.hpp"
#include "../utilities/perf.hpp"

#define N (32)

double m1[N*N];
double m2[N*N];
double m3[N*N];


void init( void )
{
    srand(9874);

    for ( int i = 0; i < N*N; i++ )
    {
        m1[i] = (double)rand() / (double)RAND_MAX;
    }
    for ( int i = 0; i < N*N; i++ )
    {
        m2[i] = (double)rand() / (double)RAND_MAX;
    }
}

void test_mm_naive( void )
{
    printf("hello\n");
}


int main( int argn, const char *argv[] )
{
    init();

    Test tests[] = {
        (Test){ test_mm_naive, "test_mm_naive" },
    };

    RunTests( tests );

    double sum = 0.0;
    for ( size_t i = 0; i < N*N; ++i ) {
        sum += m3[i];
    }

    return (int)sum;
}
