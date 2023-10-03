// rev_t.hpp
//
// A scalar-based implementation of reverse mode automatic differentiation.
//
// Main reference:
// Antoine Savine, "Modern Computational Finance - AAD and Parallel Simuations", 2018
//


#ifndef TX_REV_HPP
#define TX_REV_HPP


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <assert.h>


// -------------------- //
// Forward Declarations //
// -------------------- //


// "Modern Computational Finance - AAD and Parallel Simuations", p. 361

template <typename T>
struct Node {
    size_t     n;
    T         *derivs;     // derivatives of node wrt. parent nodes
    T        **adj_ptrs;   // pointers to adjoints of parent nodes

    T          adjoint;    // derivative of node wrt to variable that "started" the reverse pass
};

struct Tape_Position;


// ---------------------- //
// Memory Blocks for Tape //
// ---------------------- //

// Similar to the blocklist in
// "Modern Computational Finance - AAD and Parallel Simuations", p. 363

template <typename T>
struct Memory_Blocks {
    struct {
        T       **data;
        size_t    len;
        size_t    cap;
    } raw;

    size_t block_size;

    size_t cur_block;
    size_t cur_block_pos;
};


template <typename T>
void mb_print( Memory_Blocks<T> mb, const char *name )
{
    printf("(%s) = {\n", name);
    printf("\traw.len       = %zu\n", mb.raw.len);
    printf("\traw.cap       = %zu\n\n", mb.raw.cap);
    printf("\tblock_size    = %zu\n\n", mb.block_size);
    printf("\tcur_block     = %zu\n", mb.cur_block);
    printf("\tcur_block_pos = %zu\n\n", mb.cur_block_pos);
    printf("}\n");
}


// TODO: zero out rest of current block when new block is allocated

template <typename T>
T *mb_alloc(
    Memory_Blocks<T> *mblocks,
    size_t count
) {
    assert( count <= mblocks->block_size );

    if ( (mblocks->cur_block_pos + count) > mblocks->block_size ) {
        // alloc additional block
        if ( mblocks->cur_block >= mblocks->raw.len - 1 ) {
            T *new_block = (T *) malloc( mblocks->block_size * sizeof(T) );

            if ( mblocks->raw.cap == 0 || (! mblocks->raw.data) ) {
                mblocks->raw.data = (T **) malloc( sizeof(*mblocks->raw.data) );
            }
            else if ( mblocks->raw.len >= mblocks->raw.cap - 1 ) {
                size_t new_cap;
                if ( mblocks->raw.cap == 1 ) {
                    new_cap = 2;
                }
                else {
                    new_cap = (size_t) (1.618 * mblocks->raw.cap);
                }

                mblocks->raw.data = (T**) realloc(mblocks->raw.data, new_cap * sizeof(*mblocks->raw.data));
                mblocks->raw.cap  = new_cap;
            }

            mblocks->raw.data[mblocks->raw.len] = new_block;
            mblocks->raw.len++;
        }

        mblocks->cur_block    += 1;
        mblocks->cur_block_pos = 0;
    }

    T *ptr = (T *) &(mblocks->raw.data[mblocks->cur_block][mblocks->cur_block_pos]);
    mblocks->cur_block_pos += count;

    return ptr;
}


template <typename T>
void mb_reset( Memory_Blocks<T> *mblocks )
{
    mblocks->cur_block     = 0;
    mblocks->cur_block_pos = 0;
}


template <typename T>
Memory_Blocks<T> mb_make( size_t block_size )
{
    Memory_Blocks<T> mblocks;

    mblocks.block_size  = block_size;

    mblocks.raw.data    = (T **) malloc( sizeof(*mblocks.raw.data) );
    T *ptr              = (T *)  malloc( mblocks.block_size * sizeof(*ptr) );
    mblocks.raw.data[0] = ptr;
    mblocks.raw.len     = 1;
    mblocks.raw.cap     = 1;

    mblocks.cur_block     = 0;
    mblocks.cur_block_pos = 0;

    return mblocks;
}


template <typename T>
void mb_free( Memory_Blocks<T> *mblocks )
{
    for ( size_t i = 0; i < mblocks->raw.len; ++i ) {
        free( mblocks->raw.data[i] );
    }

    free( mblocks->raw.data );

    memset( mblocks, 0, sizeof(*mblocks) );
}


template <typename T>
void mb_clear( Memory_Blocks<T> *mblocks )
{
    mblocks->cur_block     = 0;
    mblocks->cur_block_pos = 0;
}


typedef struct Block_Position {
    size_t block;
    size_t element;
} Block_Position;


// next open spot in the memory block
template <typename T>
Block_Position mb_mark( Memory_Blocks<T> *mblocks )
{
    Block_Position pos;
    pos.block   = mblocks->cur_block;
    pos.element = mblocks->cur_block_pos;
    return pos;
}


template <typename T>
void mb_rewind( Memory_Blocks<T> *mblocks, Block_Position p )
{
    mblocks->cur_block     = p.block;
    mblocks->cur_block_pos = p.element;
}


// ---- //
// Node //
// ---- //

template <typename T>
Node<T> nd_make( size_t n )
{
    Node<T> node;

    node.n = n;

    node.derivs   = NULL;
    node.adj_ptrs = NULL;
    node.adjoint  = 0.0;

    return node;
}


template <typename T>
void nd_propagate_one( Node<T> node )
{
    if ( ! node.n || (node.adjoint == T(0.0)) )
        return;

    for ( size_t i=0; i<node.n; ++i ) {
        *(node.adj_ptrs[i]) += node.derivs[i] * node.adjoint;
    }
}



// ---- //
// Tape //
// ---- //

// Similar to the tape in
// "Modern Computational Finance - AAD and Parallel Simuations", p. 376

template <typename T>
struct Tape {
    Memory_Blocks<Node<T>>  nodes;
    Memory_Blocks<T>        derivs;
    Memory_Blocks<T*>       arg_ptrs;

    char _pad[64];
};


template <typename T>
void tp_print( Tape<T> t, const char * name )
{
    printf("(%s) = {\n", name);
    mb_print(t.nodes, "nodes");
    mb_print(t.derivs, "derivs");
    mb_print(t.arg_ptrs, "argPtrs");
    printf("}\n");
}


template <typename T>
Tape<T> tp_make( size_t num_nodes, size_t num_derivs )
{
    Tape<T> t;

    t.nodes    = mb_make<Node<T>>( num_nodes );
    t.derivs   = mb_make<T>( num_derivs );
    t.arg_ptrs = mb_make<T*>( num_derivs );

    return t;
}


template <typename T>
void tp_free( Tape<T> *t )
{
    mb_free( &(t->nodes)   );
    mb_free( &(t->derivs)  );
    mb_free( &(t->arg_ptrs) );
}


template <typename T>
Node<T> * tp_alloc_node( Tape<T> *t, size_t n )
{
    Node<T> *node = mb_alloc( &(t->nodes), 1 );

    if ( n > 0 ) {
        *node = nd_make<T>( n );

        node->derivs   = mb_alloc( &(t->derivs),   n );
        node->adj_ptrs = mb_alloc( &(t->arg_ptrs), n );
    }
    else {
        memset( node, 0, sizeof(*node) );
    }

    return node;
}


template <typename T>
void tp_clear( Tape<T> *t )
{
    mb_clear( &(t->nodes) );
    mb_clear( &(t->derivs) );
    mb_clear( &(t->arg_ptrs) );
}


struct Tape_Position {
    Block_Position nodes;
    Block_Position derivs;
    Block_Position arg_ptrs;
};


template <typename T>
Tape_Position tp_mark( Tape<T> *t )
{
    Tape_Position tp;
    tp.nodes    = mb_mark( &(t->nodes) );
    tp.derivs   = mb_mark( &(t->derivs) );
    tp.arg_ptrs = mb_mark( &(t->arg_ptrs) );
    return tp;
}


template <typename T>
void tp_rewind( Tape<T> *t, Tape_Position tp )
{
    mb_rewind( &(t->nodes),    tp.nodes );
    mb_rewind( &(t->derivs),   tp.derivs );
    mb_rewind( &(t->arg_ptrs), tp.arg_ptrs );
}


template <typename T>
void tp_reset_adjoints_until( Tape<T> *t, Tape_Position tp )
{
    if ( t->nodes.block_size == 0 )
        return;

    size_t n_block    = t->nodes.cur_block;
    size_t last_pos   = t->nodes.cur_block_pos;
    size_t block_size = t->nodes.block_size;

    size_t block_idx = tp.nodes.block;
    size_t node_idx  = tp.nodes.element;

    Node<T>  *node;
    Node<T> **data = (Node<T> **) t->nodes.raw.data;

    for ( ; block_idx < n_block; ++block_idx ) {
        for ( node_idx = 0; node_idx < block_size; node_idx++ ) {
            node = data[block_idx] + node_idx;
            node->adjoint = T(0);
        }
    }

    for ( node_idx = 0; node_idx < last_pos; node_idx++ ) {
        node = data[block_idx] + node_idx;
        node->adjoint = T(0);
    }
}


template <typename T>
void tp_reset_adjoints( Tape<T> *t )
{
    if ( t->nodes.block_size == 0 )
        return;

    size_t n_block    = t->nodes.cur_block;
    size_t last_pos   = t->nodes.cur_block_pos;
    size_t block_size = t->nodes.block_size;

    size_t block_idx = 0;
    size_t node_idx  = 0;

    Node<T>  *node;
    Node<T> **data = (Node<T> **) t->nodes.raw.data;

    for ( ; block_idx < n_block; ++block_idx ) {
        for ( node_idx = 0; node_idx < block_size; node_idx++ ) {
            node = data[block_idx] + node_idx;
            node->adjoint = T(0);
        }
    }

    for ( node_idx = 0; node_idx < last_pos; node_idx++ ) {
        node = data[block_idx] + node_idx;
        node->adjoint = T(0);
    }
}



// TODO: use TLS
Tape<float>    _TAPE_f = { 0 };
Tape<double>   _TAPE_d = { 0 };

const
size_t _NUM_NODES = 1e4;


template <typename T>
Tape<T> *get_tape( void )
{
    return NULL;
}

template <>
Tape<float> *get_tape( void )
{
    return &_TAPE_f;
}

template <>
Tape<double> *get_tape( void )
{
    return &_TAPE_d;
}



template <typename T>
void _TAPE_init( void )
{
    Tape<T> *t = get_tape<T>();
    *t = tp_make<T>( _NUM_NODES, 3*_NUM_NODES );
}


template <typename T>
void _TAPE_clear( void )
{
    Tape<T> *t = get_tape<T>();
    tp_clear( t );
}


template <typename T>
void acquire_tape( const Tape<T> *new_tape )
{
    Tape<T> *t = get_tape<T>();
    *t = *new_tape;
}


template <typename T>
void release_tape( Tape<T> *old )
{
    Tape<T> *t = get_tape<T>();
    old = t;
#if DEBUG
    memset(t, 0, sizeof(Tape<T>));
#endif
}


template <typename T>
Tape_Position _TAPE_mark( void )
{
    Tape<T> *t = get_tape<T>();
    return tp_mark( t );
}


template <typename T>
void _TAPE_rewind( Tape_Position tp )
{
    Tape<T> *t = get_tape<T>();
    tp_rewind( t, tp );
}


template <typename T>
void _TAPE_deinit( void )
{
    Tape<T> *t = get_tape<T>();
    tp_free( t );
}


template <typename T>
void _TAPE_reset_adjoints_until( Tape_Position tp )
{
    Tape<T> *t = get_tape<T>();
    tp_reset_adjoints_until( t, tp );
}


template <typename T>
void _TAPE_reset_adjoints( void )
{
    Tape<T> *t = get_tape<T>();
    tp_reset_adjoints( t );
}


// --------------------- //
// Reverse-mode Variable //
// --------------------- //

// Equivalent to the Number type in
// "Modern Computational Finance - AAD and Parallel Simuations", p. 379


template <typename T>
struct Rev {
    T        val;
    Node<T> *node;


    Rev( void ) {};


    // produces a leaf node on tape
    Rev( T val ) : val( val )
    {
        Tape<T> *t = get_tape<T>();
        node = tp_alloc_node( t, 0 );
    }


    // unary operation
    Rev( Node<T> *arg, T val ) : val( val )
    {
        Tape<T> *t = get_tape<T>();

        node = tp_alloc_node( t, 1 );
        node->adj_ptrs[0] = &(arg->adjoint);
    }


    // binary operation
    Rev( Node<T> *lhs, Node<T> *rhs, T val ) : val( val )
    {
        Tape<T> *t = get_tape<T>();

        node = tp_alloc_node( t, 2 );
        node->adj_ptrs[0] = &(lhs->adjoint);
        node->adj_ptrs[1] = &(rhs->adjoint);
    }


    // ternary operation
    Rev( Node<T> *n1, Node<T> *n2, Node<T> *n3, T val ) : val( val )
    {
        Tape<T> *t = get_tape<T>();

        node = tp_alloc_node( t, 3 );
        node->adj_ptrs[0] = &(n1->adjoint);
        node->adj_ptrs[1] = &(n2->adjoint);
        node->adj_ptrs[2] = &(n3->adjoint);
    }


    // explicit to avoid silent casting / conversion
    explicit operator T& ()       { return val; }
    explicit operator T  () const { return val; }


    // unary operator overloads

    Rev& operator=( T newVal );

    Rev operator-() const
    {
        return T(0.0) - *this;
    }

    Rev operator+() const
    {
        return *this;
    }


    // binary operator overloads

    Rev operator+( Rev rhs );
    Rev operator+( T   rhs );

    Rev operator-( Rev rhs );
    Rev operator-( T   rhs );

    Rev operator*( Rev rhs );
    Rev operator*( T   rhs );


    Rev operator/( Rev rhs );
    Rev operator/( T   rhs );


    Rev& operator+=( Rev rhs );
    Rev& operator+=( T   rhs );

    Rev& operator-=( Rev rhs );
    Rev& operator-=( T   rhs );

    Rev& operator*=( Rev rhs );
    Rev& operator*=( T   rhs );

    Rev& operator/=( Rev rhs );
    Rev& operator/=( T   rhs );
};


template <typename T>
void rev_print( Rev<T> v, const char *name = NULL) {
    printf("%s = {\n", name);
    printf("\tval  = %.4f\n", v.val);
    printf("\tnode = {\n");
    printf("\t\tn        = %zu\n", v.node->n);
    printf("\t\tderivs   = ["); for (size_t i = 0; i < v.node->n; i++) { printf("%.4f, ", v.node->derivs[i]); } printf("]\n");
    printf("\t\tadj_ptrs = ["); for (size_t i = 0; i < v.node->n; i++) { printf("%p, ", (void*) v.node->adj_ptrs[i]); } printf("]\n");
    printf("\t\tadj      = %.4f (%p)\n", v.node->adjoint, (void*) &v.node->adjoint);
    printf("\t}");
    printf("}\n");
}


template <typename T>
Rev<T>& Rev<T>::operator=( T newVal )
{
    Tape<T> *t = get_tape<T>();

    val  = newVal;
    node = tp_alloc_node( t, 0 );
    return *this;
}


template <typename T>
inline
Rev<T> rev_make( T val, size_t n )
{
    Tape<T> *t = get_tape<T>();

    Rev<T> v;

    v.val  = val;
    v.node = tp_alloc_node( t, n );

    return v;
}


template <typename T>
inline
Rev<T> rev_make_unary( Node<T> *arg, T val )
{
    Rev<T> v;

    v.val  = val;
    v.node = tp_alloc_node( get_tape<T>(), 1 );
    v.node->adj_ptrs[0] = &(arg->adjoint);

    return v;
}


template <typename T>
inline
Rev<T> rev_make_binary( Node<T> *lhs, Node<T> *rhs, T val )
{
    Rev<T> v;

    v.val  = val;
    v.node = tp_alloc_node( get_tape<T>(), 2 );
    v.node->adj_ptrs[0] = &(lhs->adjoint);
    v.node->adj_ptrs[1] = &(rhs->adjoint);

    return v;
}


// to produce a leaf node on tape
template <typename T>
inline
void rev_put_on_tape( Rev<T> *v )
{
    v->node = tp_alloc_node( get_tape<T>(), 0 );
}


template <typename T>
inline
void rev_put_on_tape_n( Rev<T> **a, size_t n_a )
{
    for ( size_t i = 0; i < n_a; ++i ) {
        rev_put_on_tape( a[i] );
    }
}


// access to adjoint
template <typename T>
inline
T rev_adjoint( Rev<T> v )
{
    return v.node->adjoint;
}


// accessing local derivatives with bounds check
template <typename T>
inline
T& rev_deriv( Rev<T> v )
{
    assert( v.node->n == 1 );
    return v.node->derivs[0];
}


template <typename T>
inline
T& rev_deriv_l( Rev<T> v )
{
    assert( v.node->n == 2 );
    return v.node->derivs[0];
}


template <typename T>
inline
T& rev_deriv_r( Rev<T> v )
{
    assert( v.node->n == 2 );
    return v.node->derivs[1];
}



template <typename T>
inline
T& rev_deriv_n( Rev<T> v, size_t n )
{
    assert( n < v.node->n );
    return v.node->derivs[n];
}


// accessing child node adjoints with bounds check
template <typename T>
inline
T * rev_adj( Rev<T> v )
{
    assert( v.node->n == 1 );
    return v.node->adj_ptrs[0];
}


template <typename T>
inline
T * rev_adj_l( Rev<T> v )
{
    assert( v.node->n == 2 );
    return v.node->adj_ptrs[0];
}


template <typename T>
inline
T * rev_adj_r( Rev<T> v )
{
    assert( v.node->n == 2 );
    return v.node->adj_ptrs[1];
}


template <typename T>
inline
T * rev_adj_n( Rev<T> v, size_t n )
{
    assert( n < v.node->n );
    return v.node->adj_ptrs[0];
}


// core function where back propagation of takes place
// "Modern Computational Finance - AAD and Parallel Simuations", p. 389

template <typename T>
void backprop_between( Tape_Position start, Tape_Position end )
{
    Tape<T> *tape = get_tape<T>();

#define NODE_PROPAGATION \
    node = (Node<T> *) (data[block_idx] + node_idx); \
    nd_propagate_one( *node );

#define BACKPROP_LOOP \
    for ( int64_t node_idx = block_size - 1; node_idx >= 0; node_idx-- ) { \
        NODE_PROPAGATION \
    }

    if ( tape->nodes.block_size == 0 )
        return;

    int64_t block    = (int64_t)end.nodes.block;
    int64_t position = (int64_t)end.nodes.element;

    int64_t block_idx = (int64_t)start.nodes.block;
    int64_t last_pos  = (int64_t)start.nodes.element;

    int64_t block_size = (int64_t)tape->nodes.block_size;

    int64_t node_idx  = 0;

    Node<T>   *node = NULL;
    Node<T>  **data = (Node<T> **) tape->nodes.raw.data;

    int64_t first_block_pos = 0;

    if ( block_idx == block ) {
        first_block_pos = position;
    }
    for ( node_idx = last_pos - 1; node_idx >= first_block_pos; node_idx -- ) {
        NODE_PROPAGATION
    }

    block_idx -= 1;

    for ( ; block_idx > block; --block_idx ) {
        BACKPROP_LOOP
    }

    if ( block_idx == block ) {
        BACKPROP_LOOP
    }

#undef NODE_PROPAGATION
#undef BACKPROP_LOOP
}


template <typename T>
inline
void backprop_until( Rev<T> v, Tape_Position tp )
{
    v.node->adjoint = 1.0;
    backprop_between<T>( _TAPE_mark<T>(), tp );
}


template <typename T>
inline
void backprop_to_mark( Rev<T> v, Tape_Position tp ) {
    backprop_until<T>( v, tp );
}


template <typename T>
inline
void backprop( Rev<T> v ) {
    backprop_until<T>( v, { 0 } );
}



// Binary operations
// Equivalent to
// "Modern Computational Finance - AAD and Parallel Simuations", p. 393

template <typename T>
inline
Rev<T> Rev<T>::operator+( Rev<T> rhs )
{
    T resVal = val + rhs.val;

    Rev<T> res( node, rhs.node, resVal );

    rev_deriv_l( res ) = 1.0;
    rev_deriv_r( res ) = 1.0;

    return res;
}

template <typename T>
inline
Rev<T> Rev<T>::operator+( T rhs )
{
    T resVal = val + rhs;

    Rev<T> res( node, resVal );

    rev_deriv( res ) = 1.0;

    return res;
}


template <typename T>
inline
Rev<T> operator+( T lhs, Rev<T> rhs )
{
    return rhs + lhs;
}


template <typename T>
inline
Rev<T> Rev<T>::operator-( Rev<T> rhs )
{
    T resVal = val - rhs.val;

    Rev<T> res( node, rhs.node, resVal );

    rev_deriv_l( res ) = 1.0;
    rev_deriv_r( res ) = -1.0;

    return res;
}

template <typename T>
inline
Rev<T> Rev<T>::operator-( T rhs )
{
    T resVal = val - rhs;

    Rev<T> res( node, resVal );

    rev_deriv( res ) = 1.0;

    return res;
}

template <typename T>
inline
Rev<T> operator-( T lhs, Rev<T> rhs )
{
    T resVal = lhs - rhs.val;

    Rev<T> res( rhs.node, resVal );

    rev_deriv( res ) = -1.0;

    return res;
}


template <typename T>
inline
Rev<T> Rev<T>::operator*( Rev<T> rhs )
{
    T resVal = val * rhs.val;

    Rev<T> res( node, rhs.node, resVal );

    rev_deriv_l( res ) = rhs.val;
    rev_deriv_r( res ) = val;

    // assert( std::isfinite(rev_deriv_l(res)) );
    // assert( std::isfinite(rev_deriv_r(res)) );

    return res;
}

template <typename T>
inline
Rev<T> Rev<T>::operator*( T rhs )
{
    T resVal = val * rhs;

    Rev<T> res( node, resVal );

    rev_deriv( res ) = rhs;

    // assert( std::isfinite(rev_deriv(res)) );

    return res;
}

template <typename T>
inline
Rev<T> operator*( T lhs, Rev<T> rhs )
{
    return rhs * lhs;
}


template <typename T>
inline
Rev<T> Rev<T>::operator/( Rev<T> rhs )
{
    T resVal = val / rhs.val;

    Rev<T> res( node, rhs.node, resVal );

    T irhs = 1.0 / rhs.val;

    rev_deriv_l( res ) = irhs;
    rev_deriv_r( res ) = -val * irhs * irhs;

    // assert( std::isfinite(rev_deriv_l(res)) );
    // assert( std::isfinite(rev_deriv_r(res)) );

    return res;
}

template <typename T>
inline
Rev<T> Rev<T>::operator/( T rhs )
{
    T resVal = val / rhs;

    Rev<T> res( node, resVal );

    rev_deriv( res ) = 1.0 / rhs;

    // assert( std::isfinite(rev_deriv(res)) );

    return res;
}

template <typename T>
inline
Rev<T> operator/( T lhs, Rev<T> rhs )
{
    T resVal = lhs / rhs.val;

    Rev<T> res( rhs.node, resVal );

    rev_deriv( res ) = -lhs / (rhs.val * rhs.val);

    // assert( std::isfinite(rev_deriv(res)) );

    return res;
}


template <typename T>
Rev<T>& Rev<T>::operator+=( Rev<T> rhs )
{
    *this = *this + rhs;
    return *this;
}

template <typename T>
Rev<T>& Rev<T>::operator+=( T rhs )
{
    *this = *this + rhs;
    return *this;
}

template <typename T>
Rev<T>& Rev<T>::operator-=( Rev<T> rhs )
{
    *this = *this - rhs;
    return *this;
}

template <typename T>
Rev<T>& Rev<T>::operator-=( T rhs )
{
    *this = *this - rhs;
    return *this;
}

template <typename T>
Rev<T>& Rev<T>::operator*=( Rev<T> rhs )
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
Rev<T>& Rev<T>::operator*=( T rhs )
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
Rev<T>& Rev<T>::operator/=( Rev<T> rhs )
{
    *this = *this / rhs;
    return *this;
}

template <typename T>
Rev<T>& Rev<T>::operator/=( T rhs )
{
    *this = *this / rhs;
    return *this;
}


template <typename T>
bool operator==( Rev<T> lhs, Rev<T> rhs )
{
    return lhs.val == rhs.val;
}

template <typename T>
bool operator==( Rev<T> lhs, T rhs )
{
    return lhs.val == rhs;
}

template <typename T>
bool operator==( T lhs, Rev<T> rhs )
{
    return lhs == rhs.val;
}


template <typename T>
bool operator!=( Rev<T> lhs, Rev<T> rhs )
{
    return lhs.val != rhs.val;
}

template <typename T>
bool operator!=( Rev<T> lhs, T rhs )
{
    return lhs.val != rhs;
}

template <typename T>
bool operator!=( T lhs, Rev<T> rhs )
{
    return lhs != rhs.val;
}


template <typename T>
bool operator<( Rev<T> lhs, Rev<T> rhs )
{
    return lhs.val < rhs.val;
}

template <typename T>
bool operator<( Rev<T> lhs, T rhs )
{
    return lhs.val < rhs;
}

template <typename T>
bool operator<( T lhs, Rev<T> rhs )
{
    return lhs < rhs.val;
}


template <typename T>
bool operator>( Rev<T> lhs, Rev<T> rhs )
{
    return lhs.val > rhs.val;
}

template <typename T>
bool operator>( Rev<T> lhs, T rhs )
{
    return lhs.val > rhs;
}

template <typename T>
bool operator>( T lhs, Rev<T> rhs )
{
    return lhs > rhs.val;
}


template <typename T>
bool operator<=( Rev<T> lhs, Rev<T> rhs )
{
    return lhs.val <= rhs.val;
}

template <typename T>
bool operator<=( Rev<T> lhs, T rhs )
{
    return lhs.val <= rhs;
}

template <typename T>
bool operator<=( T lhs, Rev<T> rhs )
{
    return lhs <= rhs.val;
}


template <typename T>
bool operator>=( Rev<T> lhs, Rev<T> rhs )
{
    return lhs.val >= rhs.val;
}

template <typename T>
bool operator>=( Rev<T> lhs, T rhs )
{
    return lhs.val >= rhs;
}

template <typename T>
bool operator>=( T lhs, Rev<T> rhs )
{
    return lhs >= rhs.val;
}



// ----------------------------- //
// elementary functions overload //
// ----------------------------- //


namespace std {

template <typename T>
Rev<T> log( Rev<T> x );

template <typename T>
Rev<T> pow( Rev<T> x, T a );

template <typename T>
Rev<T> pow( T x, Rev<T> a );

// template <typename T>
// Rev<T> pow( Rev<T> x, Rev<T> a );


template <typename T>
Rev<T> fabs( Rev<T> x )
{
    if ( x.val > T(0.0) ) {
        return x;
    }
    else if ( x.val < T(0.0) ) {
        return -x;
    }
    else if ( x.val == T(0.0) ) {
        return Rev<T>( T(0.0) );
    }
    else { // x.val = inf, -inf, nan, ...
        return Rev<T>( x.val );
    }
}


template <typename T>
Rev<T> sqrt( Rev<T> x )
{
    T y_val = std::sqrt( x.val );

    Rev<T> y( x.node, y_val );

    rev_deriv( y ) = T(0.5) / y_val;

    return y;
}


template <typename T>
Rev<T> pow( Rev<T> x, T a )
{
    T tmp = std::pow( x.val, a - T(1.0) );

    Rev<T> y( x.node, tmp * x.val );

    rev_deriv( y ) = a * tmp;

    return y;
}


template <typename T>
Rev<T> pow( T x, Rev<T> a )
{
    T val = std::pow( x, a.val );

    Rev<T> y( a.node, val );

    rev_deriv( y ) = std::log(x) * val;

    return y;
}


// template <typename T>
// Rev<T> pow( Rev<T> x, Rev<T> a )
// {

// }


template <typename T>
Rev<T> sin( Rev<T> x )
{
    T val = std::sin( x.val );

    Rev<T> y( x.node, val );

    rev_deriv( y ) = std::cos(x.val);

    return y;
}


template <typename T>
Rev<T> cos( Rev<T> x )
{
    T val = std::cos( x.val );

    Rev<T> y( x.node, val );

    rev_deriv( y ) = -std::sin(x.val);

    return y;
}


template <typename T>
Rev<T> tan( Rev<T> x )
{
    T tmp = std::cos( x.val );

    Rev<T> y( x.node, std::tan(x.val) );

    rev_deriv( y ) = T(1.0) / (tmp * tmp);

    return y;
}


template <typename T>
Rev<T> atan( Rev<T> x )
{
    T tmp = std::cos( x.val );

    Rev<T> y( x.node, std::atan(x.val) );

    rev_deriv( y ) = T(1.0) / (T(1.0) + (x.val * x.val));

    return y;
}


template <typename T>
Rev<T> exp( Rev<T> x )
{
    T y_val = std::exp( x.val );

    Rev<T> y(x.node, y_val);

    rev_deriv( y ) = y_val;

    return y;
}


template <typename T>
Rev<T> log( Rev<T> x )
{
    T y_val = std::log( x.val );

    Rev<T> y(x.node, y_val);

    rev_deriv( y ) = y_val != y_val ? y_val : (T(1.0) / x.val);

    return y;
}


template <typename T>
Rev<T> logabs( Rev<T> x )
{
    T y_val = std::log( std::fabs(x.val) );

    Rev<T> y(x.node, y_val);

    rev_deriv( y ) = T(1.0) / x.val;

    return y;
}


template <typename T>
Rev<T> sinh( Rev<T> x )
{
    T y_val = std::sinh( x.val );

    Rev<T> y(x.node, y_val);

    rev_deriv( y ) = std::cosh(x.val);

    return y;
}


template <typename T>
Rev<T> cosh( Rev<T> x )
{
    T y_val = std::cosh( x.val );

    Rev<T> y(x.node, y_val);

    rev_deriv( y ) = std::sinh(x.val);

    return y;
}


template <typename T>
Rev<T> tanh( Rev<T> x )
{
    T y_val = std::tanh( x.val );

    Rev<T> y(x.node, y_val);

    rev_deriv( y ) = T(1.0) - (y_val * y_val);

    return y;
}


template <typename T>
Rev<T> atanh( Rev<T> x )
{
    T y_val = std::atanh( x.val );

    Rev<T> y(x.node, y_val);

    rev_deriv( y ) = T(1.0) / (T(1.0) - (x.val * x.val));

    return y;
}


} // end namespace std


// Gradient of function / functor fulfilling prototype
// f: Rev<T> f( Rev<T> x[], size_t n )

template <typename T, typename F>
void rev_gradient_no_alloc( T *f_val, T g[], T x[], Rev<T> x_cpy[], size_t n, F f )
{
    Tape_Position tp0 = _TAPE_mark<T>();

    for ( size_t i = 0; i < n; ++i ) {
        x_cpy[i] = Rev<T>( x[i] );
    }

    Rev<T> z = f( x_cpy, n );

    backprop_to_mark<T>( z, tp0 );

    for ( size_t i = 0; i < n; ++i ) {
        g[i] = rev_adjoint( x_cpy[i] );
    }

    if ( f_val ) {
        *f_val = z.val;
    }

    _TAPE_rewind<T>( tp0 );
}

template <typename T, typename F>
void rev_gradient_no_alloc( double *f_val, double g[], double x[], Rev<double> x_cpy[], size_t n, F f )
{
    rev_gradient_no_alloc( f_val, g, x, x_cpy, n, f );
}

template <typename T, typename F>
void rev_gradient_no_alloc( float *f_val, float g[], float x[], Rev<float> x_cpy[], size_t n, F f )
{
    rev_gradient_no_alloc( f_val, g, x, x_cpy, n, f );
}


template <typename T, typename F>
void rev_gradient( T *f_val, T g[], T x[], size_t n, F f )
{
#define REV_GRAD_STACK_SIZE 50

    Rev<T> buff[REV_GRAD_STACK_SIZE];
    Rev<T> *x_cpy = buff;

    if ( n > REV_GRAD_STACK_SIZE ) {
        x_cpy = (Rev<T>*) malloc(n * sizeof(*x_cpy));
    }

    rev_gradient_no_alloc( f_val, g, x, x_cpy, n, f );

    if ( n >= REV_GRAD_STACK_SIZE ) {
        free(x_cpy);
    }

#undef REV_GRAD_STACK_SIZE
}

template <typename T, typename F>
void rev_gradient( double *f_val, double g[], double x[], size_t n, F f )
{
    rev_gradient( f_val, g, x, n, f );
}

template <typename T, typename F>
void rev_gradient( float *f_val, float g[], float x[], size_t n, F f )
{
    rev_gradient( f_val, g, x, n, f );
}

#endif
