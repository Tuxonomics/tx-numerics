// fwd_perf.cpp


#include "fwd.hpp"
#include "../../utilities/perf.hpp"


using namespace std;


double _SUM_d;
float  _SUM_f;


void test_exp_d( void )
{
    _SUM_d = 0.0;

    Fwd<double> x(2.0, 1.0);
    Fwd<double> z = exp(x);

    _SUM_d += z.dot;
}


void test_exp_f( void )
{
    _SUM_f = 0.0;

    Fwd<float> x(2.0, 1.0);
    Fwd<float> z = exp(x);

    _SUM_f += z.dot;
}



template <typename T>
T test_grad_func_1( T *x, size_t n )
{
    T res = x[0]*x[0]*x[1] + std::exp(x[1]);
    return res;
}


template <typename T>
T test_hess_1( void )
{
#define n 2

    T x[n] = { T(4.0), T(1.0) };
    T g1[n];
    T g2[n];
    T h1[n*n];
    T h2[n*n];

    T fval;

    fwd_gradient( &fval, g1, x, n, test_grad_func_1<Fwd<T>> );
    fwd_hessian( &fval, h2, g2, x, n, test_grad_func_1< Fwd<Fwd<T>> > );


    return fval;

#undef n
}

void test_hess_1_double( void )
{
    _SUM_d = test_hess_1<double>();
}

void test_hess_1_float( void )
{
    _SUM_f = test_hess_1<float>();
}



int main( int argn, const char *argv[] )
{
    Test tests[] = {
        (Test){ test_exp_d, "exp( Fwd<double> )" },
        (Test){ test_exp_f, "exp( Fwd<float> )" },

        (Test){ test_hess_1_double, "Hessian (double)" },
        (Test){ test_hess_1_float,  "Hessian (float)" },
    };

    RunTests( tests );

    return 0;
}

