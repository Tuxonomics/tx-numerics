// fwd_test.cpp

#include "fwd.hpp"
#include "../finite_differences/finite_differences.hpp"
#include "../../utilities/tx_tests.h"


using namespace std;

#define ASSERT TX_ASSERT
#define ASSERT_MSG TX_ASSERT_MSG
#define ASSERT_MSG_VA TX_ASSERT_MSG


template <typename T>
bool approx( T a, T b, T eps )
{
    if ( abs(a - b) > eps ) {
        return false;
    }
    else {
        return true;
    }
}


template <typename T>
T _eps( void )
{
    return T(NAN);
}

template <>
double _eps( void )
{
    return 1e-8;
}

template <>
float _eps( void )
{
    return 1e-4f;
}


void test_fwd( void )
{
    Fwd<double> x; x.val = 5.0; x.dot = 1.0;

    double y = double(x);

    ASSERT( x.val == y );

    Fwd<double> z = x;

    ASSERT( x.val == z.val );
}


template <typename T>
void test_fwd_add( void )
{
    Fwd<T> x(T(1.0), T(1.0));
    Fwd<T> y(T(2.0));

    ASSERT( x < y && x <= y );

    Fwd<T> z = x + y;

    ASSERT( z.val == x.val + y.val );
    ASSERT( z.dot == x.dot );

    Fwd<T> zz = z + T(4.0);

    ASSERT( zz.val == z.val + T(4.0) );
    ASSERT( zz.dot == z.dot );

    Fwd<T> zzz = T(4.0) + z;

    ASSERT( zzz.val == z.val + T(4.0) );
    ASSERT( zzz.dot == z.dot );
}

void test_fwd_sub( void )
{
    Fwd<double> x(2.0, 1.0);
    Fwd<double> y(1.0);

    Fwd<double> z = x - y;

    ASSERT( z.val == x.val - y.val );
    ASSERT( z.dot == x.dot );

    Fwd<double> zz = z - 4.0;

    ASSERT( zz.val == z.val - 4.0 );
    ASSERT( zz.dot == z.dot );

    Fwd<double> zzz = 4.0 - z;

    ASSERT( zzz.val == 4.0 - z.val );
    ASSERT( zzz.dot == -z.dot );
}

void test_fwd_mul( void )
{
    Fwd<double> x(2.0, 1.0);
    Fwd<double> y(1.0);

    Fwd<double> z = x * y;

    ASSERT( z.val == x.val * y.val );
    ASSERT( z.dot == y.val * x.dot );

    Fwd<double> zz = z * 4.0;

    ASSERT( zz.val == z.val * 4.0 );
    ASSERT( zz.dot == z.dot * 4.0 );

    Fwd<double> zzz = 4.0 * z;

    ASSERT( zzz.val == 4.0 * z.val );
    ASSERT( zzz.dot == 4.0 * z.dot );
}

void test_fwd_div( void )
{
    Fwd<double> x(2.0, 1.0);
    Fwd<double> y(3.0);

    ASSERT( y > x && y >= x );

    Fwd<double> z = x / y;

    ASSERT( z.val == x.val / y.val );
    ASSERT( z.dot == x.dot / y.val );

    Fwd<double> zz = z / 4.0;

    ASSERT( zz.val == z.val / 4.0 );
    ASSERT( zz.dot == z.dot / 4.0 );

    Fwd<double> zzz = 4.0 / z;

    ASSERT( zzz.val == 4.0 / z.val );
    ASSERT( zzz.dot == - (4.0*z.dot) / (z.val*z.val) );
}

void test_incr_decr( void )
{
    Fwd<double> x(2.0, 1.0);

    x += 3.0;

    ASSERT( x.val == 5.0 );
    ASSERT( x.dot == 1.0 );

    x -= 3.0;

    ASSERT( x.val == 2.0 );
    ASSERT( x.dot == 1.0 );
}

void test_abs( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = abs(x);

    ASSERT( z.val == x.val );
    ASSERT( z.dot == x.dot );


    z = abs(-x);

    ASSERT( z.val == x.val );
    ASSERT( z.dot == x.dot );
}

void test_sqrt( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = sqrt(x);

    ASSERT( z.val == sqrt(x.val) );
    ASSERT( z.dot == 1.0 / (2.0 * sqrt(x.val)) );
}

void test_pow( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = pow(x, 2.0);

    ASSERT( z.val == pow(x.val, 2.0) );
    ASSERT( z.dot == 2.0 * x.val );
}

void test_sin( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = sin(x);

    ASSERT( z.val == sin(x.val) );
    ASSERT( z.dot == cos(x.val) );
}

void test_cos( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = cos(x);

    ASSERT( z.val == cos(x.val) );
    ASSERT( z.dot == -sin(x.val) );
}

void test_tan( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z  = tan(x);
    Fwd<double> zz = sin(x) / cos(x);

    ASSERT( approx(z.val, zz.val, 1e-10) );
    ASSERT( approx(z.dot, zz.dot, 1e-10) );
}

void test_atan( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = atan(x);

    ASSERT( approx( z.val, atan(x.val), 1e-10) );
    ASSERT( approx( z.dot, 1.0 / (1.0 + (x.val*x.val)), 1e-10) );
}


template <typename T>
void test_exp( void )
{
    Fwd<T> x(T(2.0), T(1.0));

    Fwd<T> z = exp(x);

    ASSERT( approx(z.val, exp(x.val), T(1e-10)) );
    ASSERT( approx(z.dot, exp(x.val), T(1e-10)) );
}

void test_log( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = log(x);

    ASSERT( approx(z.val, log(x.val), 1e-10) );
    ASSERT( approx(z.dot, 1.0 / x.val, 1e-10) );
}

void test_logabs( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = logabs(-x);

    ASSERT( approx(z.val, log(abs(x.val)), 1e-10) );
    ASSERT( approx(z.dot, 1.0 / x.val, 1e-10) );
}

void test_sinh( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = sinh(x);

    ASSERT( approx(z.val, sinh(x.val), 1e-10) );
    ASSERT( approx(z.dot, cosh(x.val), 1e-10) );
}

void test_cosh( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = cosh(x);

    ASSERT( approx(z.val, cosh(x.val), 1e-10) );
    ASSERT( approx(z.dot, sinh(x.val), 1e-10) );
}

void test_tanh( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = tanh(x);

    ASSERT( approx(z.val, tanh(x.val), 1e-10) );
    ASSERT( approx(z.dot, 1.0 - tanh(x.val)*tanh(x.val), 1e-10) );
}

void test_atanh( void )
{
    Fwd<double> x(0.5, 1.0);

    Fwd<double> z = atanh(x);

    ASSERT( approx(z.val, atanh(x.val), 1e-10) );
    ASSERT( approx(z.dot, 1.0 / (1.0 - x.val*x.val), 1e-10) );
}


template <typename T>
T test_grad_func_1( T *x, size_t n )
{
    T res = x[0]*x[0]*x[1] + std::exp(x[1]);
    return res;
}


template <typename T>
void test_grad_1( T test_precision )
{
#define n 2

   T x[n] = { T(4.0), T(1.0) };
   T g1[n];
   T g2[n];

   T eps[n] = { _eps<T>(), _eps<T>() };

   T fval;

   fd_grad( &fval, g1, x, eps, n, test_grad_func_1<T> );
   fwd_gradient( &fval, g2, x, n, test_grad_func_1<Fwd<T>> );

   for ( size_t i=0; i<n; ++i ) {
       ASSERT( approx(g1[i]/g2[i], T(1.0), test_precision) );
   }

#undef n
}

void test_grad_1_double( void )
{
    test_grad_1( 1e-8 );
}

void test_grad_1_float( void )
{
    test_grad_1( 1e-2f );
}



template <typename T>
void test_hess_1( T test_precision )
{
#define n 2

    T x[n] = { T(4.0), T(1.0) };
    T g1[n];
    T g2[n];
    T h1[n*n];
    T h2[n*n];

    T eps[n] = { _eps<T>(), _eps<T>() };

    T fval;

    fwd_gradient( &fval, g1, x, n, test_grad_func_1<Fwd<T>> );
    fwd_hessian( &fval, h2, g2, x, n, test_grad_func_1< Fwd<Fwd<T>> > );

    for ( size_t i=0; i<n; ++i ) {
        ASSERT( approx(g1[i]/g2[i], T(1.0), test_precision) );
    }

#undef n
}

void test_hess_1_double( void )
{
    test_hess_1( 1e-8 );
}

void test_hess_1_float( void )
{
    test_hess_1( 1e-8f );
}



TX_TEST_LIST = {
    TX_ADD_TEST(test_fwd),
    TX_ADD_TEST(test_fwd_add<double>),
    TX_ADD_TEST(test_fwd_add<float>),
    TX_ADD_TEST(test_fwd_sub),
    TX_ADD_TEST(test_fwd_mul),
    TX_ADD_TEST(test_fwd_div),
    TX_ADD_TEST(test_incr_decr),
    TX_ADD_TEST(test_abs),
    TX_ADD_TEST(test_sqrt),
    TX_ADD_TEST(test_pow),
    TX_ADD_TEST(test_sin),
    TX_ADD_TEST(test_cos),
    TX_ADD_TEST(test_tan),
    TX_ADD_TEST(test_atan),
    TX_ADD_TEST(test_exp<double>),
    TX_ADD_TEST(test_exp<float>),
    TX_ADD_TEST(test_log),
    TX_ADD_TEST(test_logabs),
    TX_ADD_TEST(test_sinh),
    TX_ADD_TEST(test_cosh),
    TX_ADD_TEST(test_tanh),
    TX_ADD_TEST(test_atanh),

    TX_ADD_TEST(test_grad_1_double),
    TX_ADD_TEST(test_grad_1_float),

    TX_ADD_TEST(test_hess_1_double),
    TX_ADD_TEST(test_hess_1_float),
};


TX_TESTS_MAIN

