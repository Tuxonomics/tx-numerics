//
// For testing this file individually, compile with
// clang++ -o fwd_test -std=c++17 -DFWD_AD_MAIN fwd_test.cpp


#include "fwd.hpp"
#include "../finite_differences/finite_differences.hpp"


#if !defined(ASSERT)
    #include <assert.h>
    #define ASSERT assert
#endif


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
    Fwd<double> x(1.0);

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
    T res = x[0]*x[0] + std::exp(x[1]);
    return res;
}

void test_grad_1( void )
{
#define n 2

    double x[n] = { 4.0, 1.0 };
    double g1[n];
    double g2[n];

    double eps[n] = { 1e-8, 1e-8 };

    finite_differences( g1, x, eps, n, test_grad_func_1<double> );
    fwd_gradient( NULL, g2, x, n, test_grad_func_1<Fwd<double>> );

    for ( size_t i=0; i<n; ++i ) {
        ASSERT( approx(g1[i], g2[i], 1e-6) );
    }

#undef n
}


// add: Fwd<T>( val + rhs.val, dot + rhs.dot );
// mul: Fwd<T>( val * rhs.val, (val * rhs.dot) + (dot * rhs.val) )

template <typename T, typename S>
auto multivariate( T x, S y )
{
    return x*x*x; // + x*y + x;
}


template <typename T>
Fwd<T> multivariate_grad_target( Fwd<T> *x, size_t n )
{
    return multivariate( x[0], x[1] );
}


#if FWD_AD_MAIN

#include "stdio.h"


int main( int argn, const char *argv[] )
{
    test_fwd();
    test_fwd_add();
    test_fwd_sub();
    test_fwd_mul();
    test_fwd_div();
    test_incr_decr();
    test_sqrt();
    test_pow();
    test_sin();
    test_cos();
    test_tan();
    test_atan();
    test_exp();
    test_log();
    test_logabs();
    test_sinh();
    test_cosh();
    test_tanh();
    test_atanh();
    test_grad_1();

    {
        Fwd<double> x(3.0, 1.0);
        Fwd<double> y(5.0);

        Fwd<double> z = multivariate(x, y);

        printf("z.val = %.4f, z.dot = %.4f\n\n", z.val, z.dot);


        printf("Computing the gradient\n");
#define n 2

        double g[n];
        double vals[n] = {3.0, 5.0};

        fwd_gradient(NULL, g, vals, n, multivariate_grad_target<double>);

        printf("d.f/d.x = %.4f \t d.f/d.y = %.4f\n", g[0], g[1]);
    }

    {
        Fwd<Fwd<double>> x = {0};
        x.val.val = 4.0; x.val.dot = 1.0;
        x.dot.val = 1.0; x.dot.dot = 0.0;

        Fwd<Fwd<double>> y = {0};
        y.val.val = 5.0;


        Fwd<Fwd<double>> z = multivariate(x, y);

        printf("z.val.val = %.4f, z.val.dot = %.4f\nz.dot.val = %.4f, z.dot.dot = %.4f\n\n",
            z.val.val, z.val.dot, z.dot.val, z.dot.dot);

    }
    return 0;
}
#endif

