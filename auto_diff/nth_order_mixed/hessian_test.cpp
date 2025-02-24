// hessian_test.cpp

#include "../../utilities/tx_tests.hpp"
#include "hessian.hpp"

using namespace std;

#define ASSERT TX_ASSERT
#define ASSERT_MSG TX_ASSERT_MSG
#define ASSERT_MSG_VA TX_ASSERT_MSG

#define APPROX TX_APPROX
#define APPROX_ARRAY TX_APPROX_ARRAY


template <typename T>
T test_hessian_func_1( T *x, size_t n )
{
    T res = x[0]*x[0]*x[1] + std::exp(x[1]);
    return res;
}

template <typename T>
void test_hessian_1( T test_precision )
{
#define n 2
    _TAPE_init<Fwd<T>>();

    T x[n] = { T(4.0), T(1.0) };

    T fval;
    T g[n];
    T h[n*n];

    T fval_test;
    T g_test[n];
    T h_test[n*n];

    mixed_hessian( &fval, h, g, x, n, test_hessian_func_1<Rev<Fwd<T>>> );

    fwd_hessian( &fval_test, h_test, g_test, x, n, test_hessian_func_1<Fwd<Fwd<T>>> );

    _TAPE_deinit<Fwd<T>>();

    ASSERT( TX_APPROX( fval/fval_test, T(1.0), test_precision) );
    ASSERT( TX_APPROX_ARRAY(g, g_test, test_precision, n, T) );
    ASSERT( TX_APPROX_ARRAY(h, h_test, test_precision, n, T) );

#undef n
}

void test_hessian_1_double( void )
{
    test_hessian_1( 1e-10 );
}

void test_hessian_1_float( void )
{
    test_hessian_1( 1e-10f );
}


template <typename T>
struct test_hessian_func_2 {
    T a;

    test_hessian_func_2(T a) : a(a) {}

    template <typename S>
    S operator () ( S *x, size_t n ) {
        S res =
            x[0]*x[0]*x[1]
            + std::exp(x[1])
            + std::tanh(x[2]) * std::log(x[3])
            + std::pow(x[3], a)
        ;
        return res;
    }
};

template <typename T>
void test_hessian_2( T test_precision )
{
#define n 4
    _TAPE_init<Fwd<T>>();

    T x[n] = { T(4.0), T(1.0), T(0.1), T(3.0) };

    T fval;
    T g[n];
    T h[n*n];

    T fval_test;
    T g_test[n];
    T h_test[n*n];

    test_hessian_func_2<T> test_fun( T(3.5) );

    mixed_hessian( &fval, h, g, x, n, test_fun );

    fwd_hessian( &fval_test, h_test, g_test, x, n, test_fun );

    _TAPE_deinit<Fwd<T>>();

    ASSERT( TX_APPROX( fval/fval_test, T(1.0), test_precision) );
    ASSERT( TX_APPROX_ARRAY(g, g_test, test_precision, n, T) );
    ASSERT( TX_APPROX_ARRAY(h, h_test, test_precision, n, T) );

#undef n
}

void test_hessian_2_double( void )
{
    test_hessian_2( 1e-10 );
}

void test_hessian_2_float( void )
{
    test_hessian_2( 1e-6f );
}


TX_TEST_LIST = {
    TX_ADD_TEST(test_hessian_1_double),
    TX_ADD_TEST(test_hessian_1_float),

    TX_ADD_TEST(test_hessian_2_double),
    TX_ADD_TEST(test_hessian_2_float),
};


TX_TESTS_MAIN

