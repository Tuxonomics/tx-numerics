// rng_test.cpp

#include "../utilities/tx_tests.h"
#include "lin_alg.hpp"


#define ASSERT TX_ASSERT
#define ASSERT_MSG TX_ASSERT_MSG
#define ASSERT_MSG_VA TX_ASSERT_MSG

#define APPROX TX_APPROX
#define APPROX_ARRAY TX_APPROX_ARRAY


#define N (32)

double m1[N*N];
double m2[N*N];
double m3[N*N];
double m4[N*N];


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
    init();

    printf("hello\n");
}


TX_TEST_LIST = {
    TX_ADD_TEST(test_mm_naive),
};


TX_TESTS_MAIN
