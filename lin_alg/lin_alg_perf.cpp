// rng_test.cpp

#include <stdlib.h>
#include "lin_alg.hpp"
#include "../utilities/perf.hpp"
#include <Accelerate/Accelerate.h>

#define N (512)

double m1[N*N];
double m2[N*N];
double m3[N*N];
double test[N*N];


void init( void )
{
    printf("Running matrix size %dx%d.\n", N, N);

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
    mm_naive(test, m1, m2, N, N, N);
}

void test_mm_linear( void )
{
    mm_linear(test, m1, m2, N, N, N);
}

void test_mm_linear_x4( void )
{
    mm_linear_x4(test, m1, m2, N, N, N);
}

void test_mm_linear_x4_ptr( void )
{
    mm_linear_x4_ptr(test, m1, m2, N, N, N);
}

void test_benchmark_dgemm( void )
{
    cblas_dgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        N, N, N, 1.0, m1, N, m2, N, 0.0, test, N
    );
}

int main( int argn, const char *argv[] )
{
    init();

    Test tests[] = {
        (Test){ test_mm_naive, "test_mm_naive" },
        (Test){ test_mm_linear, "test_mm_linear" },
        (Test){ test_mm_linear_x4, "test_mm_linear_x4" },
        (Test){ test_mm_linear_x4_ptr, "test_mm_linear_x4_ptr" },
        (Test){ test_benchmark_dgemm, "test_benchmark_dgemm" },
    };

    RunTests( tests );

    double sum = 0.0;
    for ( size_t i = 0; i < N*N; ++i ) {
        sum += m3[i];
    }

    return (int)sum;
}


