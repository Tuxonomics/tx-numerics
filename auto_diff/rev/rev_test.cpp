// rev_t_test.cpp

#include "../../utilities/tx_tests.h"
#include "rev.hpp"
#include "../finite_differences/finite_differences.hpp"


using namespace std;

#define ASSERT TX_ASSERT
#define ASSERT_MSG TX_ASSERT_MSG
#define ASSERT_MSG_VA TX_ASSERT_MSG

#define APPROX TX_APPROX


template <typename T>
T _eps( void )
{
    return T(NAN);
}

template <>
double _eps( void )
{
    return 1e-10;
}

template <>
float _eps( void )
{
    return 1e-4f;
}


#ifndef EPS
    #define EPS 1E-10
#endif


void test_memory_blocks( void )
{
#define T char

    size_t n_blocks = 1024;

    Memory_Blocks<T> mb = mb_make<T>( n_blocks );

    T *ar1 = mb_alloc( &mb, n_blocks/2 );
    T *ar2 = mb_alloc( &mb, n_blocks/2 );
    T *ar3 = mb_alloc( &mb, n_blocks/2 );

    (void) ar2;
    (void) ar3;

    ASSERT( mb.raw.len == 2 );

    Block_Position bp = mb_mark( &mb );

    T *ar4 = mb_alloc( &mb, n_blocks/2 );
    T *ar5 = mb_alloc( &mb, n_blocks/2 );

    ASSERT( mb.raw.len == 3 );

    mb_rewind( &mb, bp );

    T *ar6 = mb_alloc( &mb, n_blocks/2 );
    T *ar7 = mb_alloc( &mb, n_blocks/2 );

    ASSERT( mb.raw.len == 3 );
    ASSERT( ar4 == ar6 );
    ASSERT( ar5 == ar7 );

//    mb_print(mb, NULL);


    mb_clear( &mb );

//    mb_print(mb, NULL);

    T *ar8 = mb_alloc( &mb, n_blocks/2 );

    ASSERT( ar8 == ar1 );
    ASSERT( mb.raw.len == 3 );


    mb_free( &mb );

//    mb_print(mb, NULL);


#undef T
}


void test_node( void )
{
#define T char
    Node<T> node = nd_make<T>( 0 );

    ASSERT( node.adjoint == 0 );

    nd_propagate_one( node );
#undef T
}


void test_tape( void )
{
#define T double

    size_t nNodes = 2;

    Tape<T> t = tp_make<T>( nNodes, 2*nNodes );

    Node<T> *n0 = tp_alloc_node( &t, 0 );
    Node<T> *n1 = tp_alloc_node( &t, 0 );
    Node<T> *n2 = tp_alloc_node( &t, 0 );


    ASSERT( n0->adjoint == 0 );
    ASSERT( n1->adjoint == 0 );
    ASSERT( n2->adjoint == 0 );


    n0->adjoint = 4.0;
    n1->adjoint = -1.0;
    n2->adjoint = 10.0;

    tp_reset_adjoints( &t );


    ASSERT( n0->adjoint == 0 );
    ASSERT( n1->adjoint == 0 );
    ASSERT( n2->adjoint == 0 );


    tp_clear( &t );

    Node<T> *n3 = tp_alloc_node( &t, 0 );

    ASSERT( n3 == n0 );


    tp_free( &t );

#undef T
}


template <typename T>
void test_rev_alloc( void )
{
    _TAPE_init<T>();

    Rev<T> x = 1;

    _TAPE_deinit<T>();
}



template <typename T>
void test_basic_ops( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(1.0));

    ASSERT( x.val == T(1.0) && rev_adjoint(x) == T(0.0) );
    ASSERT( x.node->n == 0 );

    _TAPE_clear<T>();

    x = T(2.0);

    Rev<T> z = T(2.0) * x;

    backprop( z );

    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_reset_adjoints<T>();

    ASSERT( rev_adjoint(z) == T(0.0) );
    ASSERT( rev_adjoint(x) == T(0.0) );

    backprop( z );

    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_add( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(2.0) + x;

    backprop( z );

    ASSERT( APPROX(z.val, T(4.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x + T(2.0);

    backprop( z );

    ASSERT( APPROX(z.val, T(4.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(3.0));

    z = y + x;

    backprop( z );

    ASSERT( APPROX(z.val, T(5.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(y), T(1.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_subtract( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(3.0) - x;

    backprop( z );

    ASSERT( APPROX(z.val, T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(-1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x - T(3.0);

    backprop( z );

    ASSERT( APPROX(z.val, T(-1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(3.0));

    z = y - x;

    backprop( z );

    ASSERT( APPROX(z.val, T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(-1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(y), T(1.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_multiply( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(2.0) * x;

    backprop( z );

    ASSERT( APPROX( z.val, T(4.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x * T(2.0);

    backprop( z );

    ASSERT( APPROX( z.val, T(4.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(3.0));

    z = y * x;

    backprop( z );

    ASSERT( APPROX( z.val, T(6.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(3.0), EPS) );
    ASSERT( APPROX(rev_adjoint(y), T(2.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_divide( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(4.0) / x;

    backprop( z );

    ASSERT( APPROX( z.val, T(2.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(-1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x / T(4.0);

    backprop( z );

    ASSERT( APPROX( z.val, T(0.5), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(0.25), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(4.0));

    z = y / x;

    backprop( z );

    ASSERT( APPROX( z.val, T(2.0), EPS) );
    ASSERT( APPROX(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(x), T(-1.0), EPS) );
    ASSERT( APPROX(rev_adjoint(y), T(0.5), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_tape_operations( void )
{
    _TAPE_init<T>();

    T xval = T(3.0);
    T yval = T(2.0);

    Tape_Position tp0 = _TAPE_mark<T>();

    Rev<T> x(xval);
    Rev<T> y(yval);

    Rev<T> z = x*y;

    T **adj_ptr_y = &z.node->adj_ptrs[1];

    Tape_Position tp1 = _TAPE_mark<T>();

    int n = 4;

    for ( int i = 0; i < n; ++i ) {
        Rev<T> v = T(2.0)*z; // error
        Rev<T> w = T(2.0)*v;

        backprop_to_mark<T>( w, tp1 );

        _TAPE_rewind<T>( tp1 );
    }

    backprop_between<T>( tp1, tp0 );

    ASSERT( APPROX( rev_adjoint(x), T(4.0)*n*yval, EPS ) );
    ASSERT( APPROX( rev_adjoint(y), T(4.0)*n*xval, EPS ) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_fabs( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = fabs(x);

    backprop( z );

    ASSERT( APPROX(z.val,          x.val, T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), 1.0,   T(1e-10)) );

    _TAPE_reset_adjoints<T>();


    z = fabs(-x);

    backprop( z );

    ASSERT( APPROX(z.val,          x.val, T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), 1.0,   T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_sqrt( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = sqrt(x);

    backprop( z );

    ASSERT( APPROX(z.val,          sqrt(x.val),                     T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0) / (T(2.0) * sqrt(x.val)), T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_pow_1( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    T a = T(2.0);

    Rev<T> z = pow(x, a);

    backprop( z );

    ASSERT( APPROX(z.val,          pow(x.val, a),              T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), a * pow(x.val, a - T(1.0)), T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_sin( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = sin(x);

    backprop( z );

    ASSERT( APPROX(z.val,          sin(x.val), T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), cos(x.val), T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_cos( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = cos(x);

    backprop( z );

    ASSERT( APPROX(z.val,          cos(x.val),  T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), -sin(x.val), T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_tan( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = tan(x);

    backprop( z );

    ASSERT( APPROX(z.val,          tan(x.val),                         T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0) / (cos(x.val) * cos(x.val)), T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_atan( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = atan(x);

    backprop( z );

    ASSERT( APPROX(z.val,          atan(x.val),                       T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0) / (T(1.0) + (x.val*x.val)), T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_exp( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = exp(x);

    backprop( z );

    ASSERT( APPROX(z.val,          exp(x.val), T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), exp(x.val), T(1e-10)) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_log( void )
{
    _TAPE_init<T>();

    Rev<T> x( T(2.0) );

    Rev<T> z = log(x);

    backprop( z );

    ASSERT( APPROX(z.val,          log(x.val),     T(1e-10)) );
    ASSERT( APPROX(rev_adjoint(x), T(1.0) / x.val, T(1e-10)) );

    _TAPE_deinit<T>();
}


// template <typename T>
// void test_logabs( void )
// {
//     Fwd<T> x(T(2.0), T(1.0));

//     Fwd<T> z = logabs(-x);

//     ASSERT( APPROX(z.val, log(abs(x.val)), T(1e-10)) );
//     ASSERT( APPROX(z.dot, T(1.0) / x.val,  T(1e-10)) );
// }

// template <typename T>
// void test_sinh( void )
// {
//     Fwd<T> x(T(2.0), T(1.0));

//     Fwd<T> z = sinh(x);

//     ASSERT( APPROX(z.val, sinh(x.val), T(1e-10)) );
//     ASSERT( APPROX(z.dot, cosh(x.val), T(1e-10)) );
// }

// template <typename T>
// void test_cosh( void )
// {
//     Fwd<T> x(T(2.0), T(1.0));

//     Fwd<T> z = cosh(x);

//     ASSERT( APPROX(z.val, cosh(x.val), T(1e-10)) );
//     ASSERT( APPROX(z.dot, sinh(x.val), T(1e-10)) );
// }

// template <typename T>
// void test_tanh( void )
// {
//     Fwd<T> x(T(2.0), T(1.0));

//     Fwd<T> z = tanh(x);

//     ASSERT( APPROX(z.val, tanh(x.val), T(1e-10)) );
//     ASSERT( APPROX(z.dot, 1.0 - tanh(x.val)*tanh(x.val), T(1e-10)) );
// }

// template <typename T>
// void test_atanh( void )
// {
//     Fwd<T> x(T(0.5), T(1.0));

//     Fwd<T> z = atanh(x);

//     ASSERT( APPROX(z.val, atanh(x.val), T(1e-10)) );
//     ASSERT( APPROX(z.dot, T(1.0) / (T(1.0) - x.val*x.val), T(1e-10)) );
// }


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
    _TAPE_init<T>();

    T x[n] = { T(4.0), T(1.0) };
    T g1[n];
    T g2[n];

    T eps[n] = { _eps<T>(), _eps<T>() };

    T fval;

    fd_grad( &fval, g1, x, eps, n, test_grad_func_1<T> );
    rev_gradient( &fval, g2, x, n, test_grad_func_1<Rev<T>> );

    for ( size_t i=0; i<n; ++i ) {
       ASSERT( APPROX(g1[i]/g2[i], T(1.0), test_precision) );
    }

    _TAPE_deinit<T>();

#undef n
}

void test_grad_1_double( void )
{
    test_grad_1( 1e-6 );
}

void test_grad_1_float( void )
{
    test_grad_1( 1e-2f );
}


TX_TEST_LIST = {
    TX_ADD_TEST(test_memory_blocks),
    TX_ADD_TEST(test_node),
    TX_ADD_TEST(test_tape),

    TX_ADD_TEST(test_rev_alloc<double>),
    TX_ADD_TEST(test_rev_alloc<float>),

    TX_ADD_TEST(test_basic_ops<double>),
    TX_ADD_TEST(test_basic_ops<float>),

    TX_ADD_TEST(test_add<double>),
    TX_ADD_TEST(test_add<float>),
    TX_ADD_TEST(test_subtract<double>),
    TX_ADD_TEST(test_subtract<float>),
    TX_ADD_TEST(test_multiply<double>),
    TX_ADD_TEST(test_multiply<float>),
    TX_ADD_TEST(test_divide<double>),
    TX_ADD_TEST(test_divide<float>),

    TX_ADD_TEST(test_fabs<double>),
    TX_ADD_TEST(test_fabs<float>),

    TX_ADD_TEST(test_sqrt<double>),
    TX_ADD_TEST(test_sqrt<float>),

    TX_ADD_TEST(test_pow_1<double>),
    TX_ADD_TEST(test_pow_1<float>),

    TX_ADD_TEST(test_sin<double>),
    TX_ADD_TEST(test_sin<float>),

    TX_ADD_TEST(test_cos<double>),
    TX_ADD_TEST(test_cos<float>),

    TX_ADD_TEST(test_tan<double>),
    TX_ADD_TEST(test_tan<float>),

    TX_ADD_TEST(test_atan<double>),
    TX_ADD_TEST(test_atan<float>),

    TX_ADD_TEST(test_exp<double>),
    TX_ADD_TEST(test_exp<float>),

    TX_ADD_TEST(test_log<double>),
    TX_ADD_TEST(test_log<float>),

    TX_ADD_TEST(test_tape_operations<double>),
    TX_ADD_TEST(test_tape_operations<float>),

    TX_ADD_TEST(test_grad_1_double),
    TX_ADD_TEST(test_grad_1_float),
};


TX_TESTS_MAIN

