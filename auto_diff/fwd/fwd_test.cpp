//
// For testing this file individually, compile with
// clang++ -o fwd_test -std=c++17 -DFWD_AD_MAIN fwd_test.cpp


#include "fwd.hpp"
#include "../finite_differences/finite_differences.hpp"


using namespace std;


bool approx( double a, double b, double eps )
{
    if ( abs(a - b) > eps ) {
        return false;
    }
    else {
        return true;
    }
}


void test_fwd( void )
{
    Fwd<double> x; x.val = 5.0; x.dot = 1.0;

    double y = double(x);

    ASSERT( x.val == y );

    Fwd<double> z = x;

    ASSERT( x.val == z.val );
}

void test_fwd_add( void )
{
    Fwd<double> x(1.0, 1.0);
    Fwd<double> y(2.0);

    ASSERT( x < y && x <= y );

    Fwd<double> z = x + y;

    ASSERT( z.val == x.val + y.val );
    ASSERT( z.dot == x.dot );

    Fwd<double> zz = z + 4.0;

    ASSERT( zz.val == z.val + 4.0 );
    ASSERT( zz.dot == z.dot );

    Fwd<double> zzz = 4.0 + z;

    ASSERT( zzz.val == z.val + 4.0 );
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

void test_exp( void )
{
    Fwd<double> x(2.0, 1.0);

    Fwd<double> z = exp(x);

    ASSERT( approx(z.val, exp(x.val), 1e-10) );
    ASSERT( approx(z.dot, exp(x.val), 1e-10) );
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

void test_grad_1( void )
{
#define n 2

    double x[n] = { 4.0, 1.0 };
    double g1[n];
    double g2[n];

    double eps[n] = { 1e-8, 1e-8 };

    fd_grad( NULL, g1, x, eps, n, test_grad_func_1<double> );
    fwd_gradient( NULL, g2, x, n, test_grad_func_1<Fwd<double>> );

    for ( size_t i=0; i<n; ++i ) {
        ASSERT( approx(g1[i], g2[i], 1e-6) );
    }

#undef n
}

#include <stdio.h>

void test_hess_1( void )
{
#define n 2

    double x[n] = { 4.0, 1.0 };
    double g1[n];
    double g2[n];
    double h1[n*n];
    double h2[n*n];

    double eps[n] = { 1e-8, 1e-8 };

    fwd_gradient( NULL, g1, x, n, test_grad_func_1<Fwd<double>> );
    fwd_hessian( NULL, h2, g2, x, n, test_grad_func_1< Fwd<Fwd<double>> > );

    for ( size_t i=0; i<n; ++i ) {
        ASSERT( approx(g1[i], g2[i], 1e-6) );
    }

    // for ( size_t i = 0; i < n*n; i++ ) {
    //     printf("h[%zu] = %.4f\n", i, h2[i]);
    // }

#undef n
}


test tests[N_TESTS] = {
    ADD_TEST(test_fwd),
    ADD_TEST(test_fwd_add),
    ADD_TEST(test_fwd_sub),
    ADD_TEST(test_fwd_mul),
    ADD_TEST(test_fwd_div),
    ADD_TEST(test_incr_decr),
    ADD_TEST(test_abs),
    ADD_TEST(test_sqrt),
    ADD_TEST(test_pow),
    ADD_TEST(test_sin),
    ADD_TEST(test_cos),
    ADD_TEST(test_tan),
    ADD_TEST(test_atan),
    ADD_TEST(test_exp),
    ADD_TEST(test_log),
    ADD_TEST(test_logabs),
    ADD_TEST(test_sinh),
    ADD_TEST(test_cosh),
    ADD_TEST(test_tanh),
    ADD_TEST(test_atanh),
    ADD_TEST(test_grad_1),
    ADD_TEST(test_hess_1),
};

