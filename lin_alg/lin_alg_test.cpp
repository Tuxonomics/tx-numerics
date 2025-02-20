// rng_test.cpp

#include "../utilities/tx_tests.h"
#include "lin_alg.hpp"
#include <Accelerate/Accelerate.h>


#define ASSERT TX_ASSERT
#define ASSERT_MSG TX_ASSERT_MSG
#define ASSERT_MSG_VA TX_ASSERT_MSG

#define APPROX TX_APPROX
#define APPROX_ARRAY TX_APPROX_ARRAY

void print_matrix(double *mat, int n, int m, const char *name)
{
    printf("%s:\n", name);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            printf("%6.2f ", mat[i + j * n]);
        }
        printf("\n");
    }
}

#define n 5
#define m 7
#define p 6

double a[n*m] = { 0.28757752, 0.78830514, 0.40897692, 0.88301740, 0.94046728, 0.04555650, 0.52810549,
0.89241904, 0.55143501, 0.45661474, 0.95683335, 0.45333416, 0.67757064, 0.57263340,
0.10292468, 0.89982497, 0.24608773, 0.04205953, 0.32792072, 0.95450365, 0.88953932,
0.69280341, 0.64050681, 0.99426978, 0.65570580, 0.70853047, 0.54406602, 0.59414202,
0.28915974, 0.14711365, 0.96302423, 0.90229905, 0.69070528, 0.79546742, 0.02461368, };

double b[m*p] = {0.4777959711, 0.7584595375, 0.2164079358, 0.3181810076, 0.2316257854,
0.1428000224, 0.4145463358, 0.4137243263, 0.3688454509, 0.1524447477,
0.1388060634, 0.2330340995, 0.4659624503, 0.2659726404, 0.8578277153,
0.0458311667, 0.4422000742, 0.7989248456, 0.1218992600, 0.5609479838,
0.2065313896, 0.1275316502, 0.7533078643, 0.8950453592, 0.3744627759,
0.6651151946, 0.0948406609, 0.3839696378, 0.2743836446, 0.8146400389,
0.4485163414, 0.8100643530, 0.8123895095, 0.7943423211, 0.4398316876,
0.7544751586, 0.6292211316, 0.7101824014, 0.0006247733, 0.4753165741,
0.2201188852, 0.3798165377,
};

double c[n*p];

double test[n*p] = {1.371766, 1.565810, 1.551816, 1.669752, 1.3047436,
1.200126, 1.279145, 1.217315, 1.279541, 0.9335935,
2.095562, 1.673502, 1.278969, 1.745646, 1.8033127,
2.292969, 1.855114, 2.094199, 2.157503, 1.3729823,
2.983127, 2.441049, 2.473277, 2.601277, 2.1097713,
1.870270, 2.040924, 2.048884, 2.258442, 1.4239615,
};


void test_mm_naive( void )
{
    mm_naive(c, a, b, n, m, p);
    mm_naive(c, a, b, n, m, p);

    print_matrix(c, n, p, "c");
    print_matrix(test, n, p, "test");

    ASSERT(APPROX_ARRAY(c, test, 1e-6, n*p, double));
    memset(c, 0, sizeof(*c)*n*p);
}

void test_mm_linear( void )
{
    mm_linear(c, a, b, n, m, p);
    mm_linear(c, a, b, n, m, p);

    print_matrix(c, n, p, "c");
    print_matrix(test, n, p, "test");

    ASSERT(APPROX_ARRAY(c, test, 1e-6, n*p, double));
    memset(c, 0, sizeof(*c)*n*p);
}

void test_mm_linear_x4( void )
{
    mm_linear_x4(c, a, b, n, m, p);
    mm_linear_x4(c, a, b, n, m, p);

    print_matrix(c, n, p, "c");
    print_matrix(test, n, p, "test");

    ASSERT(APPROX_ARRAY(c, test, 1e-6, n*p, double));
    memset(c, 0, sizeof(*c)*n*p);
}

void test_mm_linear_x4_ptr( void )
{
    mm_linear_x4_ptr(c, a, b, n, m, p);
    mm_linear_x4_ptr(c, a, b, n, m, p);

    print_matrix(c, n, p, "c");
    print_matrix(test, n, p, "test");

    ASSERT(APPROX_ARRAY(c, test, 1e-6, n*p, double));
    memset(c, 0, sizeof(*c)*n*p);
}

void test_dgemm( void )
{
    cblas_dgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        n, p, m, 1.0, a, n, b, m, 0.0, c, n
    );

    print_matrix(c, n, p, "c");
    print_matrix(test, n, p, "test");

    ASSERT(APPROX_ARRAY(c, test, 1e-6, n*p, double));
    memset(c, 0, sizeof(*c)*n*p);
}


TX_TEST_LIST = {
    TX_ADD_TEST(test_mm_naive),
    TX_ADD_TEST(test_mm_linear),
    TX_ADD_TEST(test_mm_linear_x4),
    TX_ADD_TEST(test_mm_linear_x4_ptr),
    TX_ADD_TEST(test_dgemm),
};


TX_TESTS_MAIN
