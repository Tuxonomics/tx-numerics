// rev_t_test.cpp

#include "rev_t.hpp"
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

    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_reset_adjoints<T>();

    ASSERT( rev_adjoint(z) == T(0.0) );
    ASSERT( rev_adjoint(x) == T(0.0) );

    backprop( z );

    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_add( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(2.0) + x;

    backprop( z );

    ASSERT( approx(z.val, T(4.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x + T(2.0);

    backprop( z );

    ASSERT( approx(z.val, T(4.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(3.0));

    z = y + x;

    backprop( z );

    ASSERT( approx(z.val, T(5.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(y), T(1.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_subtract( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(3.0) - x;

    backprop( z );

    ASSERT( approx(z.val, T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(-1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x - T(3.0);

    backprop( z );

    ASSERT( approx(z.val, T(-1.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(3.0));

    z = y - x;

    backprop( z );

    ASSERT( approx(z.val, T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(-1.0), EPS) );
    ASSERT( approx(rev_adjoint(y), T(1.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_multiply( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(2.0) * x;

    backprop( z );

    ASSERT( approx( z.val, T(4.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x * T(2.0);

    backprop( z );

    ASSERT( approx( z.val, T(4.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(2.0), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(3.0));

    z = y * x;

    backprop( z );

    ASSERT( approx( z.val, T(6.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(3.0), EPS) );
    ASSERT( approx(rev_adjoint(y), T(2.0), EPS) );

    _TAPE_deinit<T>();
}


template <typename T>
void test_divide( void )
{
    _TAPE_init<T>();

    Rev<T> x(T(2.0));

    Rev<T> z = T(4.0) / x;

    backprop( z );

    ASSERT( approx( z.val, T(2.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(-1.0), EPS) );

    _TAPE_reset_adjoints<T>();


    z = x / T(4.0);

    backprop( z );

    ASSERT( approx( z.val, T(0.5), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(0.25), EPS) );

    _TAPE_reset_adjoints<T>();


    Rev<T> y(T(4.0));

    z = y / x;

    backprop( z );

    ASSERT( approx( z.val, T(2.0), EPS) );
    ASSERT( approx(rev_adjoint(z), T(1.0), EPS) );
    ASSERT( approx(rev_adjoint(x), T(-1.0), EPS) );
    ASSERT( approx(rev_adjoint(y), T(0.5), EPS) );

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

    ASSERT( approx( rev_adjoint(x), T(4.0)*n*yval, EPS ) );
    ASSERT( approx( rev_adjoint(y), T(4.0)*n*xval, EPS ) );

    _TAPE_deinit<T>();
}


test tests[N_TESTS] = {
    ADD_TEST(test_memory_blocks),
    ADD_TEST(test_node),
    ADD_TEST(test_tape),

    ADD_TEST(test_rev_alloc<double>),
    ADD_TEST(test_rev_alloc<float>),

    ADD_TEST(test_basic_ops<double>),
    ADD_TEST(test_basic_ops<float>),

    ADD_TEST(test_add<double>),
    ADD_TEST(test_add<float>),
    ADD_TEST(test_subtract<double>),
    ADD_TEST(test_subtract<float>),
    ADD_TEST(test_multiply<double>),
    ADD_TEST(test_multiply<float>),
    ADD_TEST(test_divide<double>),
    ADD_TEST(test_divide<float>),

    ADD_TEST(test_tape_operations<double>),
    ADD_TEST(test_tape_operations<float>),
};
