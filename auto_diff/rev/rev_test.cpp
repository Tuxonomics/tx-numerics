// rev_test.cpp

#include "rev.hpp"
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

    size_t nBlocks = 1024;

    Memory_Blocks mb = mb_make( sizeof(T), nBlocks );

    T *ar1 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );
    T *ar2 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );
    T *ar3 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );

    (void) ar2;
    (void) ar3;

    ASSERT( mb.raw.len == 2 );

    Block_Position bp = mb_mark( &mb );

    T *ar4 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );
    T *ar5 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );

    ASSERT( mb.raw.len == 3 );

    mb_rewind( &mb, bp );

    T *ar6 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );
    T *ar7 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );

    ASSERT( mb.raw.len == 3 );
    ASSERT( ar4 == ar6 );
    ASSERT( ar5 == ar7 );

//    mb_print(mb, NULL);


    mb_clear( &mb );

//    mb_print(mb, NULL);

    T *ar8 = (T*) mb_alloc( &mb, nBlocks/2, sizeof(T) );

    ASSERT( ar8 == ar1 );
    ASSERT( mb.raw.len == 3 );


    mb_free( &mb );

//    mb_print(mb, NULL);


#undef TYPE
}


void test_node( void )
{
    Node node = nd_make( 0 );

    ASSERT( node.adjoint == 0 );

    nd_propagate_one( node );
}


void test_tape( void )
{
    size_t nNodes = 2;

    Tape t = tp_make( nNodes, 2*nNodes );

    Node *n0 = (Node*) tp_alloc_node( &t, 0 );
    Node *n1 = (Node*) tp_alloc_node( &t, 0 );
    Node *n2 = (Node*) tp_alloc_node( &t, 0 );


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

    Node *n3 = (Node*) tp_alloc_node( &t, 0 );

    ASSERT( n3 == n0 );


    tp_free( &t );
}


void test_rev_alloc( void )
{
    _TAPE_init();

    Rev x = 1;
    // rev_print(x);

    _TAPE_deinit();
}


void test_basic_ops( void )
{
    _TAPE_init();

    Rev x(1.0);

    ASSERT( x.val == 1.0 && rev_adjoint(x) == 0.0 );
    ASSERT( x.node->n == 0 );

    _TAPE_clear();

    x = 2.0;

    Rev z = 2.0 * x;

    backprop( z );

    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 2.0, EPS) );

    _TAPE_reset_adjoints();

    ASSERT( rev_adjoint(z) == 0.0 );
    ASSERT( rev_adjoint(x) == 0.0 );

    backprop( z );

    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 2.0, EPS) );

    _TAPE_deinit();
}


void test_add( void )
{
    _TAPE_init();

    Rev x(2.0);

    Rev z = 2.0 + x;

    backprop( z );

    ASSERT( approx(z.val, 4.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 1.0, EPS) );

    _TAPE_reset_adjoints();


    z = x + 2.0;

    backprop( z );

    ASSERT( approx(z.val, 4.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 1.0, EPS) );

    _TAPE_reset_adjoints();


    Rev y(3.0);

    z = y + x;

    backprop( z );

    ASSERT( approx(z.val, 5.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(y), 1.0, EPS) );

    _TAPE_deinit();
}


void test_subtract( void )
{
    _TAPE_init();

    Rev x(2.0);

    Rev z = 3.0 - x;

    backprop( z );

    ASSERT( approx(z.val, 1.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), -1.0, EPS) );

    _TAPE_reset_adjoints();


    z = x - 3.0;

    backprop( z );

    ASSERT( approx(z.val, -1.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 1.0, EPS) );

    _TAPE_reset_adjoints();


    Rev y(3.0);

    z = y - x;

    backprop( z );

    ASSERT( approx(z.val, 1.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), -1.0, EPS) );
    ASSERT( approx(rev_adjoint(y), 1.0, EPS) );

    _TAPE_deinit();
}


void test_multiply( void )
{
    _TAPE_init();

    Rev x(2.0);

    Rev z = 2.0 * x;

    backprop( z );

    ASSERT( approx( z.val, 4.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 2.0, EPS) );

    _TAPE_reset_adjoints();


    z = x * 2.0;

    backprop( z );

    ASSERT( approx( z.val, 4.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 2.0, EPS) );

    _TAPE_reset_adjoints();


    Rev y(3.0);

    z = y * x;

    backprop( z );

    ASSERT( approx( z.val, 6.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 3.0, EPS) );
    ASSERT( approx(rev_adjoint(y), 2.0, EPS) );

    _TAPE_deinit();
}


void test_divide( void )
{
    _TAPE_init();

    Rev x(2.0);

    Rev z = 4.0 / x;

    backprop( z );

    ASSERT( approx( z.val, 2.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), -1.0, EPS) );

    _TAPE_reset_adjoints();


    z = x / 4.0;

    backprop( z );

    ASSERT( approx( z.val, 0.5, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), 0.25, EPS) );

    _TAPE_reset_adjoints();


    Rev y(4.0);

    z = y / x;

    backprop( z );

    ASSERT( approx( z.val, 2.0, EPS) );
    ASSERT( approx(rev_adjoint(z), 1.0, EPS) );
    ASSERT( approx(rev_adjoint(x), -1.0, EPS) );
    ASSERT( approx(rev_adjoint(y), 0.5, EPS) );

    _TAPE_deinit();
}


void test_tape_operations( void )
{

    _TAPE_init();

    double xval = 3.0;
    double yval = 2.0;

    Tape_Position tp0 = _TAPE_mark();

    Rev x(xval);
    rev_print(x, "x");
    Rev y(yval);
    rev_print(y, "y");

    Rev z = x*y;
    rev_print(z, "z");
    
    double **adj_ptr_y = &z.node->adj_ptrs[1];
    printf("adj_ptr_y = %p\n", (void*) *adj_ptr_y);

    Tape_Position tp1 = _TAPE_mark();

    int n = 4;

    for ( int i = 0; i < n; ++i ) {
        Rev v = 2.0*z; // error
        
        printf("\n\n-------------------\n[%i]\n", i);
        rev_print(v, "v");
        
        Rev w = 2.0*v;
        rev_print(w, "w");

        backprop_to_mark( w, tp1 );

        printf("----\n");
        rev_print(v, "v");
        rev_print(w, "w");
        rev_print(z, "z");

        _TAPE_rewind( tp1 );
    }

    backprop_between( tp1, tp0 );
    
    printf("\n\n-------------------\n");
    rev_print(x, "x");
    rev_print(y, "y");
    rev_print(z, "z");

    ASSERT( approx( rev_adjoint(x), 4.0*n*yval, EPS ) );
    ASSERT( approx( rev_adjoint(y), 4.0*n*xval, EPS ) );

    _TAPE_deinit();
}


test tests[N_TESTS] = {
    ADD_TEST(test_memory_blocks),
    ADD_TEST(test_node),
    ADD_TEST(test_tape),
    ADD_TEST(test_rev_alloc),
    ADD_TEST(test_basic_ops),
    ADD_TEST(test_add),
    ADD_TEST(test_subtract),
    ADD_TEST(test_multiply),
    ADD_TEST(test_divide),
    ADD_TEST(test_tape_operations),
};
