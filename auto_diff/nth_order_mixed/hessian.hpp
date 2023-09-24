// hessian.cpp

#ifndef TX_HESSIAN_HPP
#define TX_HESSIAN_HPP

#include "../fwd/fwd.hpp"
#include "../rev/rev.hpp"


Tape<Fwd<double>> _TAPE_fd = { 0 };

template <>
Tape<Fwd<double>> *get_tape( void )
{
    return &_TAPE_fd;
}


// Hessian of function / functor fulfilling prototype
// f: Rev<Fwd<T>> f( Rev<Fwd<T>> x[], size_t n )
// h is column-major nxn matrix
template <typename T, typename F>
void mixed_hessian_no_alloc( T *f_val, T h[], T g[], T x[], Rev<Fwd<T>> x_cpy[], size_t n, F f )
{
	int f_val_not_set = 1;

	Tape_Position tp0 = _TAPE_mark<Fwd<T>>();

    for ( size_t i = 0; i < n; ++i ) {
        x_cpy[i] = Rev<Fwd<T>>( Fwd<T>(x[i]) );
    }

    Tape_Position tp1 = _TAPE_mark<Fwd<T>>();

    for ( size_t j = 0; j < n; ++j ) {
        if (j > 0 ) {
             _TAPE_reset_adjoints_until<Fwd<T>>( tp0 );
        }

        x_cpy[j].val.dot = 1.0;

        Rev<Fwd<T>> z = test_hessian_func_1( x_cpy, n );

        backprop_to_mark<Fwd<T>>( z, tp0 );

        for ( size_t i = 0; i < n; ++i ) {
            Fwd<T> adj = rev_adjoint( x_cpy[i] );
            if ( j == 0 ) {
                g[i] = adj.val;
            }

            if ( i >= j ) {
                h[i + j*n] = adj.dot;
                if ( i != j ) {
                    h[j + i*n] = adj.dot;
                }
            }

        }

        if ( f_val && f_val_not_set ) {
            *f_val = z.val.val;
            f_val_not_set = 0;
        }

        x_cpy[j].val.dot = 0.0;

        _TAPE_rewind<Fwd<T>>( tp1 );
    }

    _TAPE_rewind<Fwd<T>>( tp0 );
}

template <typename T, typename F>
void mixed_hessian_no_alloc( double *f_val, double h[], double g[], double x[], Rev<Fwd<double>> x_cpy[], size_t n, F f )
{
    mixed_hessian_no_alloc( f_val, h, g, x, x_cpy, n, f );
}

template <typename T, typename F>
void mixed_hessian_no_alloc( float *f_val, float h[], float g[], float x[], Rev<Fwd<float>> x_cpy[], size_t n, F f )
{
    mixed_hessian_no_alloc( f_val, h, g, x, x_cpy, n, f );
}


template <typename T, typename F>
void mixed_hessian( T *f_val, T h[], T g[], T x[], size_t n, F f )
{
#define MXD_HESS_STACK_SIZE 2500

    Rev<Fwd<T>> buff[MXD_HESS_STACK_SIZE];
    Rev<Fwd<T>> *x_cpy = buff;

    if ( n*n >= MXD_HESS_STACK_SIZE ) {
        x_cpy = (Rev<Fwd<T>>*) malloc(n * sizeof(*x_cpy));
    }

    mixed_hessian_no_alloc( f_val, h, g, x, x_cpy, n, f );

    if ( n*n >= MXD_HESS_STACK_SIZE ) {
        free(x_cpy);
    }

#undef MXD_HESS_STACK_SIZE
}

template <typename T, typename F>
void mixed_hessian( double *f_val, double h[], double g[], double x[], size_t n, F f )
{
    mixed_hessian( f_val, h, g, x, n, f );
}

template <typename T, typename F>
void mixed_hessian( float *f_val, float h[], float g[], float x[], size_t n, F f )
{
    mixed_hessian( f_val, h, g, x, n, f );
}

#endif

