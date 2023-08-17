// rev_test.cpp

#include "rev.hpp"
#include "../finite_differences/finite_differences.hpp"


using namespace std;


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
    rev_print(x);

    _TAPE_deinit();
}


test tests[N_TESTS] = {
    ADD_TEST(test_memory_blocks),
    ADD_TEST(test_node),
    ADD_TEST(test_tape),
    ADD_TEST(test_rev_alloc),
};
