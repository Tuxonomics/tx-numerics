// fwd.hpp
//
// A simple scalar-based implementation of forward mode automatic differentiation.
// This is for testing and applications where a small number of parameters
// can be expected.


#include <stdlib.h>  // malloc, free
#include <cmath>


// --------------------- //
// Forward-mode Variable //
// --------------------- //

template <typename T>
struct Fwd {
    T val;
    T dot;


    Fwd( void ) {};

    Fwd( double val ) : val( val ) { dot = T(0.0); };

    Fwd( Fwd<double> x ) : val( x.val ), dot( x.dot ) {};

    Fwd( Fwd<Fwd<double>> x ) : val( x.val ), dot( x.dot ) {};

    Fwd( T val, T dot ) : val( val ), dot( dot ) {};


    explicit operator T& ()       { return val; }
    explicit operator T  () const { return val; }


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


    Fwd operator-( Fwd rhs );
    Fwd operator-( double rhs );


    Fwd operator*( Fwd rhs );
    Fwd operator*( double rhs );


    Fwd operator/( Fwd rhs );
    Fwd operator/( double rhs );


    Fwd& operator+=( Fwd rhs );
    Fwd& operator+=( double rhs );

    Fwd& operator-=( Fwd rhs );
    Fwd& operator-=( double rhs );

    Fwd& operator*=( Fwd rhs );
    Fwd& operator*=( double rhs );

    Fwd& operator/=( Fwd rhs );
    Fwd& operator/=( double rhs );
};


// ------------------------ //
// unary operator overloads //
// ------------------------ //

template <typename T>
Fwd<T>& Fwd<T>::operator=( double new_val )
{
    val = T(new_val);
    dot = T(0.0);
    return *this;
}

template <typename T>
Fwd<T> Fwd<T>::operator-()
{
    return Fwd( -val, -dot );
}


// ------------------------- //
// binary operator overloads //
// ------------------------- //

template <typename T>
Fwd<T> Fwd<T>::operator+( Fwd<T> rhs )
{
    return Fwd<T>( val + rhs.val, dot + rhs.dot );
}

template <typename T>
Fwd<T> Fwd<T>::operator+( double rhs )
{
    return Fwd<T>( val + rhs, dot );
}

template <typename T>
Fwd<T> operator+( double lhs, Fwd<T> rhs )
{
    return rhs + lhs;
}


template <typename T>
Fwd<T> Fwd<T>::operator-( Fwd<T> rhs )
{
    return Fwd<T>( val - rhs.val, dot - rhs.dot );
}

template <typename T>
Fwd<T> Fwd<T>::operator-( double rhs )
{
    return Fwd<T>( val - rhs, dot );
}

template <typename T>
Fwd<T> operator-( double lhs, Fwd<T> rhs )
{
    return Fwd<T>( lhs - rhs.val, -rhs.dot );
}


template <typename T>
Fwd<T> Fwd<T>::operator*( Fwd<T> rhs )
{
    return Fwd<T>( val * rhs.val, (val * rhs.dot) + (dot * rhs.val) );
}

template <typename T>
Fwd<T> Fwd<T>::operator*( double rhs )
{
    return Fwd<T>( val * rhs, dot * rhs );
}

template <typename T>
Fwd<T> operator*( double lhs, Fwd<T> rhs )
{
    return rhs * lhs;
}


template <typename T>
Fwd<T> Fwd<T>::operator/( Fwd<T> rhs )
{
    return Fwd<T>( val / rhs.val, ((dot*rhs.val) - (val*rhs.dot)) / (rhs.val*rhs.val));
}

template <typename T>
Fwd<T> Fwd<T>::operator/( double rhs )
{
    return Fwd<T>( val / rhs, dot / rhs );
}

template <typename T>
Fwd<T> operator/( double lhs, Fwd<T> rhs )
{
    return Fwd<T>( lhs / rhs.val, - (lhs*rhs.dot) / (rhs.val*rhs.val) );
}


template <typename T>
Fwd<T>& Fwd<T>::operator+=( Fwd<T> rhs )
{
    *this = *this + rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator+=( double rhs )
{
    *this = *this + rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator-=( Fwd<T> rhs )
{
    *this = *this - rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator-=( double rhs )
{
    *this = *this - rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator*=( Fwd<T> rhs )
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator*=( double rhs )
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator/=( Fwd<T> rhs )
{
    *this = *this / rhs;
    return *this;
}

template <typename T>
Fwd<T>& Fwd<T>::operator/=( double rhs )
{
    *this = *this / rhs;
    return *this;
}


// ----------------------------- //
// comparison operator overloads //
// ----------------------------- //


template <typename T>
bool operator==( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val == rhs.val;
}

template <typename T>
bool operator==( Fwd<T> lhs, double rhs )
{
    return lhs.val == rhs;
}

template <typename T>
bool operator==( double lhs, Fwd<T> rhs )
{
    return lhs == rhs.val;
}


template <typename T>
bool operator!=( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val != rhs.val;
}

template <typename T>
bool operator!=( Fwd<T> lhs, double rhs )
{
    return lhs.val != rhs;
}

template <typename T>
bool operator!=( double lhs, Fwd<T> rhs )
{
    return lhs != rhs.val;
}


template <typename T>
bool operator<( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val < rhs.val;
}

template <typename T>
bool operator<( Fwd<T> lhs, double rhs )
{
    return lhs.val < rhs;
}

template <typename T>
bool operator<( double lhs, Fwd<T> rhs )
{
    return lhs < rhs.val;
}


template <typename T>
bool operator>( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val > rhs.val;
}

template <typename T>
bool operator>( Fwd<T> lhs, double rhs )
{
    return lhs.val > rhs;
}

template <typename T>
bool operator>( double lhs, Fwd<T> rhs )
{
    return lhs > rhs.val;
}


template <typename T>
bool operator<=( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val <= rhs.val;
}

template <typename T>
bool operator<=( Fwd<T> lhs, double rhs )
{
    return lhs.val <= rhs;
}

template <typename T>
bool operator<=( double lhs, Fwd<T> rhs )
{
    return lhs <= rhs.val;
}


template <typename T>
bool operator>=( Fwd<T> lhs, Fwd<T> rhs )
{
    return lhs.val >= rhs.val;
}

template <typename T>
bool operator>=( Fwd<T> lhs, double rhs )
{
    return lhs.val >= rhs;
}

template <typename T>
bool operator>=( double lhs, Fwd<T> rhs )
{
    return lhs >= rhs.val;
}



// ----------------------------- //
// elementary functions overload //
// ----------------------------- //


namespace std {

template <typename T>
Fwd<T> abs( Fwd<T> x )
{
    if ( x.val > 0.0 ) {
        return x;
    }
    else if ( x.val < 0.0 ) {
        return -x;
    }
    else if ( x.val == 0.0 ) {
        return Fwd(T(0.0), T(0.0));
    }
    else { // x.val = inf, -inf, nan, ...
        return Fwd(abs(x.val), T(NAN));
    }
}

template <typename T>
Fwd<T> sqrt( Fwd<T> x )
{
    T tmp = std::sqrt( x.val );
    return Fwd<T>( tmp, (T(0.5)*x.dot) / tmp );
}

template <typename T>
Fwd<T> pow( Fwd<T> x, double a )
{
    T tmp = std::pow( x.val, a - 1.0 );
    return Fwd<T>( x.val * tmp, T(a) * tmp * x.dot );
}

template <typename T>
Fwd<T> sin( Fwd<T> x )
{
    return Fwd<T>( std::sin(x.val), std::cos(x.val)*x.dot );
}

template <typename T>
Fwd<T> cos( Fwd<T> x )
{
    return Fwd<T>( std::cos(x.val), -std::sin(x.val)*x.dot );
}

template <typename T>
Fwd<T> tan( Fwd<T> x )
{
    T tmp = std::cos(x.val);
    return Fwd<T>( std::tan(x.val), x.dot / (tmp*tmp) );
}

template <typename T>
Fwd<T> atan( Fwd<T> x )
{
    return Fwd<T>( std::atan(x.val), x.dot / (T(1.0) + (x.val*x.val)) );
}

template <typename T>
Fwd<T> exp( Fwd<T> x )
{
    T tmp = std::exp(x.val);
    return Fwd<T>( tmp, tmp * x.dot );
}

template <typename T>
Fwd<T> log( Fwd<T> x )
{
    T tmp = std::log(x.val);
    T dot = tmp != tmp ? tmp : (x.dot / x.val);
    return Fwd<T>( tmp, dot );
}

template <typename T>
Fwd<T> logabs( Fwd<T> x )
{
    return Fwd<T>( std::log(std::abs(x.val)), x.dot / x.val );
}

template <typename T>
Fwd<T> sinh( Fwd<T> x )
{
    return Fwd<T>( std::sinh(x.val), x.dot * std::cosh(x.val) );
}

template <typename T>
Fwd<T> cosh( Fwd<T> x )
{
    return Fwd<T>( std::cosh(x.val), x.dot * std::sinh(x.val) );
}

template <typename T>
Fwd<T> tanh( Fwd<T> x )
{
    T tmp = std::tanh(x.val);
    return Fwd<T>( tmp, x.dot * (T(1.0) - (tmp * tmp)) );
}

template <typename T>
Fwd<T> atanh( Fwd<T> x )
{
    return Fwd<T>( std::atanh(x.val), T(1.0) / (T(1.0) - (x.val*x.val)) );
}

} // end namespace std



// Gradient from function / functor fulfilling prototype
// f: Fwd<double> f( Fwd<double> x[], size_t n )


template <typename F>
void fwd_gradient_no_alloc( double *f_val, double g[], double x[], Fwd<double> x_cpy[], size_t n, F f )
{
    Fwd<double> tmp = { 0 };
    int f_val_not_set = 1;

    for ( size_t i=0; i<n; ++i ) {
        x_cpy[i] = Fwd<double> (x[i]);
    }

    for ( size_t i=0; i<n; ++i ) {
        x_cpy[i].dot = 1.0;

        tmp = f( x_cpy, n );

        if (f_val_not_set && f_val != NULL) {
            *f_val = tmp.val;
            f_val_not_set = 0;
        }

        g[i] = tmp.dot;

        x_cpy[i].dot = 0.0;
    }
}


template <typename F>
void fwd_gradient( double *f_val, double g[], double x[], size_t n, F f )
{
#define FWD_GRAD_STACK_SIZE 50

    Fwd<double> buff[FWD_GRAD_STACK_SIZE];
    Fwd<double> *x_cpy = buff;

    if ( n >= FWD_GRAD_STACK_SIZE ) {
        x_cpy = (Fwd<double>*) malloc(n * sizeof(*x_cpy));
    }

    fwd_gradient_no_alloc( f_val, g, x, x_cpy, n, f );

    if ( n >= FWD_GRAD_STACK_SIZE ) {
        free(x_cpy);
    }

#undef FWD_GRAD_STACK_SIZE
}



// Hessian from function / functor fulfilling prototype
// f: Fwd<Fwd<double>> f( Fwd<Fwd<double>> x[], size_t n )
// h is column-major nxn matrix

template <typename F>
void fwd_hessian_no_alloc( double *f_val, double h[], double g[], double x[], Fwd<Fwd<double>> x_cpy[], size_t n, F f )
{
    Fwd<Fwd<double>> tmp = { 0 };
    int f_val_not_set = 1;

    for ( size_t i=0; i<n; ++i ) {
        x_cpy[i] = Fwd<Fwd<double>>( Fwd<double>(x[i]) );
    }

    for ( size_t i=0; i<n; ++i ) {
        x_cpy[i].val.dot = 1.0;

        for ( size_t j=i; j<n; ++j ) {
            x_cpy[j].dot.val = 1.0;

            tmp = f( x_cpy, n );

            if ( f_val_not_set && f_val != NULL ) {
                *f_val = tmp.val.val;
                f_val_not_set = 0;
            }

            if ( j == i && g != NULL ) {
                g[i] = tmp.val.dot;
            }

            h[i + j*n] = tmp.dot.dot;
            h[j + i*n] = tmp.dot.dot;

            x_cpy[j].dot.val = 0.0;
        }

        x_cpy[i].val.dot = 0.0;
    }
}


template <typename F>
void fwd_hessian( double *f_val, double h[], double g[], double x[], size_t n, F f )
{
#define FWD_HESS_STACK_SIZE 50

    Fwd<Fwd<double>> buff[FWD_HESS_STACK_SIZE];
    Fwd<Fwd<double>> *x_cpy = buff;

    if ( n >= FWD_HESS_STACK_SIZE ) {
        x_cpy = (Fwd<Fwd<double>>*) malloc(n * sizeof(*x_cpy));
    }

    fwd_hessian_no_alloc( f_val, h, g, x, x_cpy, n, f );

    if ( n >= FWD_HESS_STACK_SIZE ) {
        free(x_cpy);
    }

#undef FWD_HESS_STACK_SIZE
}



