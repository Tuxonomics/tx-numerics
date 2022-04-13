// fwd_draft.cpp
// clang++ -o fwd_draft -std=c++17 fwd_draft.cpp

#include "stdio.h"
#include <stdlib.h>  // malloc, free
#include <cmath>


struct Fwd {
    double val;
    double dot;


    Fwd( void ) {};

    Fwd( double val ) : val( val ) { dot = 0.0; };

    Fwd( const Fwd &x ) : val( x.val ), dot( x.dot ) {};

    Fwd( double val, double dot ) : val( val ), dot( dot ) {};


    explicit operator double& ()       { return val; }
    explicit operator double  () const { return val; }


    // unary operator overloads

    Fwd& operator=( double new_val );
    Fwd  operator-();

    Fwd  operator+() const
    {
        return *this;
    }


    // binary operator overloads

    Fwd operator+( Fwd rhs );
    Fwd operator+( double rhs );

    Fwd operator*( Fwd rhs );
    Fwd operator*( double rhs );
};


// ------------------------ //
// unary operator overloads //
// ------------------------ //

Fwd& Fwd::operator=( double new_val )
{
    val = new_val;
    dot = 0.0;
    return *this;
}

Fwd Fwd::operator-()
{
    return Fwd( -val, -dot );
}


// ------------------------- //
// binary operator overloads //
// ------------------------- //

Fwd Fwd::operator+( Fwd rhs )
{
    return Fwd( val + rhs.val, dot + rhs.dot );
}

Fwd Fwd::operator+( double rhs )
{
    return Fwd( val + rhs, dot );
}

Fwd operator+( double lhs, Fwd rhs )
{
    return rhs + lhs;
}


Fwd Fwd::operator*( Fwd rhs )
{
    return Fwd( val * rhs.val, (val * rhs.dot) + (dot * rhs.val) );
}

Fwd Fwd::operator*( double rhs )
{
    return Fwd( val * rhs, dot * rhs );
}

Fwd operator*( double lhs, Fwd rhs )
{
    return rhs * lhs;
}



// ----------------------------- //
// elementary functions overload //
// ----------------------------- //

namespace std {

Fwd exp( Fwd x )
{
    double tmp = std::exp(x.val);
    return Fwd( tmp, tmp * x.dot );
}

Fwd log( Fwd x )
{
    double tmp = std::log(x.val);
    double dot = tmp != tmp ? tmp : (x.dot / x.val);
    return Fwd( tmp, dot );
}

}



// prototype for f: Fwd f( Fwd *x, size_t n )
template <typename F>
void fwd_gradient_no_alloc( double *g, double *x, Fwd *x_cpy, size_t n, F f )
{
    Fwd tmp = { 0 };

    for ( size_t i=0; i<n; ++i ) {
        x_cpy[i] = Fwd(x[i]);
    }

    for ( size_t i=0; i<n; ++i ) {
        x_cpy[i].dot = 1.0;

        tmp = f( x_cpy, n );

        g[i] = tmp.dot;

        x_cpy[i].dot = 0.0;
    }
}


template <typename F>
void fwd_gradient( double *g, double *x, size_t n, F f )
{
    Fwd *x_cpy = (Fwd*) malloc(n * sizeof(Fwd));

    fwd_gradient_no_alloc( g, x, x_cpy, n, f );

    free(x_cpy);
}




// -------- //
// Examples //
// -------- //


Fwd univariate( Fwd x )
{
	return x*x;
}

template <typename T, typename S>
auto multivariate( T x, S y )
{
	return std::log(std::exp(x*y + x));
}


Fwd multivariate_grad_target( Fwd *x, size_t n )
{
	return multivariate( x[0], x[1] );
}


int main( int argn, const char *argv[] )
{
	printf("Univariate example\n");

	Fwd x(5.0, 1.0);

	Fwd z = univariate(x);

	printf("z.val = %.4f, z.dot = %.4f\n\n", z.val, z.dot);


	printf("Multivariate example\n");

	x.val = 3.0;
	x.dot = 1.0;

	Fwd y(5.0);

	z = multivariate(x, y);

	printf("z.val = %.4f, z.dot = %.4f\n\n", z.val, z.dot);

	z = multivariate(x, 5.0);

	printf("z.val = %.4f, z.dot = %.4f\n\n", z.val, z.dot);


	printf("Computing the gradient\n");
#define n 2

	double g[n];
	double vals[n] = {3.0, 5.0};

	fwd_gradient(g, vals, n, multivariate_grad_target);

	printf("d.f/d.x = %.4f \t d.f/d.y = %.4f\n", g[0], g[1]);

}


