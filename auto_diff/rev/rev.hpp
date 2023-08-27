// rev.hpp
//
// A scalar-based implementation of reverse mode automatic differentiation.
//
// Main reference:
// Antoine Savine, "Modern Computational Finance - AAD and Parallel Simuations", 2018
//

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

typedef struct Node {
    size_t     n;
    double    *derivs;     // derivatives of node wrt. parent nodes
    double   **adj_ptrs;   // pointers to adjoints of parent nodes

    double     adjoint;    // derivative of node wrt to variable that "started" the reverse pass
} Node;


// ---------------------- //
// Memory Blocks for Tape //
// ---------------------- //

// Similar to the blocklist in
// "Modern Computational Finance - AAD and Parallel Simuations", p. 363

typedef struct Memory_Blocks {
    struct {
        unsigned char **data;
        size_t    len;
        size_t    cap;
    } raw;

    size_t elementSize;
    size_t numElements;
    size_t arrayCap;

    size_t curBlock;
    size_t curBlockPos;
} Memory_Blocks;


void mb_print( Memory_Blocks mb, const char *name )
{
    printf("(%s) = {\n", name);
    printf("\traw.len   = %zu\n", mb.raw.len);
    printf("\traw.cap   = %zu\n\n", mb.raw.cap);
    printf("\telementSize = %zu\n", mb.elementSize);
    printf("\tnumElements = %zu\n", mb.numElements);
    printf("\tarrayCap    = %zu\n\n", mb.arrayCap);
    printf("\tcurBlock    = %zu\n", mb.curBlock);
    printf("\tcurBlockPos = %zu\n\n", mb.curBlockPos);
    printf("}\n");
}


void *mb_alloc(
    Memory_Blocks *mblocks,
    size_t count,
    size_t size
) {
    assert( count <= mblocks->arrayCap );
    assert( count % mblocks->elementSize == 0 );

    if ( mblocks->curBlockPos + count > mblocks->arrayCap ) {
        // alloc additional block
        if ( mblocks->curBlock >= mblocks->raw.len - 1 ) {
            unsigned char *new_block = (unsigned char *) malloc( mblocks->arrayCap );

            if ( mblocks->raw.cap == 0 || (! mblocks->raw.data) ) {
                mblocks->raw.data = (unsigned char **) malloc( sizeof(*mblocks->raw.data) );
            }
            else if ( mblocks->raw.len >= mblocks->raw.cap - 1 ) {
                size_t new_cap;
                if ( mblocks->raw.cap == 1 ) {
                    new_cap = 2;
                }
                else {
                    new_cap = (size_t) (1.618 * mblocks->raw.cap);
                }

                mblocks->raw.data = (unsigned char **) realloc(mblocks->raw.data, new_cap * sizeof(*mblocks->raw.data));
                mblocks->raw.cap  = new_cap;
            }

            mblocks->raw.data[mblocks->raw.len] = new_block;
            mblocks->raw.len++;
        }

        mblocks->curBlock    += 1;
        mblocks->curBlockPos  = 0;
    }

    unsigned char *ptr = (unsigned char *) &(mblocks->raw.data[mblocks->curBlock][mblocks->curBlockPos]);
    mblocks->curBlockPos += count;

    return ptr;
}


void mb_reset( Memory_Blocks *mblocks )
{
    mblocks->curBlock    = 0;
    mblocks->curBlockPos = 0;
}


Memory_Blocks mb_make( size_t elementSize, size_t numElements )
{
    Memory_Blocks mblocks;

    mblocks.elementSize = elementSize;
    mblocks.numElements = numElements;
    mblocks.arrayCap    = (elementSize * numElements);

    mblocks.raw.data    = (unsigned char **) malloc( sizeof(*mblocks.raw.data) );
    unsigned char *ptr  = (unsigned char *)  malloc( mblocks.arrayCap * sizeof(*ptr) );
    mblocks.raw.data[0] = ptr;
    mblocks.raw.len     = 1;
    mblocks.raw.cap     = 1;

    mblocks.curBlock    = 0;
    mblocks.curBlockPos = 0;

    return mblocks;
}


void mb_free( Memory_Blocks *mblocks )
{
    for ( size_t i=0; i < mblocks->raw.len; ++i ) {
        free( mblocks->raw.data[i] );
    }

    free( mblocks->raw.data );

    memset( mblocks, 0, sizeof(*mblocks) );
}


void mb_clear( Memory_Blocks *mblocks )
{
    mblocks->curBlock    = 0;
    mblocks->curBlockPos = 0;
}


typedef struct Block_Position {
    size_t block;
    size_t element;
} Block_Position;


// NOTE: this position refers to to the next open spot in the memory block.
Block_Position mb_mark( Memory_Blocks *mblocks )
{
    Block_Position pos;
    pos.block   = mblocks->curBlock;
    pos.element = mblocks->curBlockPos;
    return pos;
}


void mb_rewind( Memory_Blocks *mblocks, Block_Position p )
{
    mblocks->curBlock    = p.block;
    mblocks->curBlockPos = p.element;
}


// ---- //
// Node //
// ---- //

Node nd_make( size_t n )
{
    Node node;

    node.n = n;

    node.derivs   = NULL;
    node.adj_ptrs = NULL;
    node.adjoint  = 0.0;

    return node;
}


void nd_propagate_one( Node node )
{
    if ( ! node.n || (node.adjoint == 0.0) )
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

typedef struct Tape {
    Memory_Blocks   nodes;
    Memory_Blocks   derivs;
    Memory_Blocks   arg_ptrs;

    char            _pad[64];
} Tape;


void tp_print( Tape t, const char * name )
{
    printf("(%s) = {\n", name);
    mb_print(t.nodes, "nodes");
    mb_print(t.derivs, "derivs");
    mb_print(t.arg_ptrs, "argPtrs");
    printf("}\n");
}


Tape tp_make( size_t num_nodes, size_t num_derivs )
{
    Tape t;

    t.nodes    = mb_make( sizeof(Node),    num_nodes );
    t.derivs   = mb_make( sizeof(double),  num_derivs );
    t.arg_ptrs = mb_make( sizeof(double*), num_derivs );

    return t;
}


void tp_free( Tape *t )
{
    mb_free( &(t->nodes)   );
    mb_free( &(t->derivs)  );
    mb_free( &(t->arg_ptrs) );
}


Node * tp_alloc_node( Tape *t, size_t n )
{
    Node *node = (Node*) mb_alloc( &(t->nodes), sizeof(*node), 1 );

    if ( n > 0 ) {
        *node = nd_make( n );

        node->derivs   = (double*)  mb_alloc( &(t->derivs),   sizeof(*(node->derivs)), n );
        node->adj_ptrs = (double**) mb_alloc( &(t->arg_ptrs), sizeof(*(node->adj_ptrs)), n );
    }
    else {
        memset( node, 0, sizeof(*node) );
    }

    return node;
}


void tp_reset_adjoints( Tape *t )
{
    size_t nBlocks = t->nodes.curBlock;
    size_t cap     = t->nodes.arrayCap;

    if ( cap == 0 )
        return;

    size_t nodeSize = t->nodes.elementSize;
    size_t lastPos  = t->nodes.curBlockPos;

    size_t blockIdx = 0;
    size_t nodeIdx  = 0;

    Node  *node;
    unsigned char **data = t->nodes.raw.data;

    assert( cap == nodeSize * t->nodes.numElements );

    for ( ; blockIdx < nBlocks; ++blockIdx ) {
        for ( nodeIdx = 0; nodeIdx < cap; nodeIdx += nodeSize ) {
            node = (Node *) (data[blockIdx] + nodeIdx);
            node->adjoint = 0;
        }
    }

    for ( nodeIdx = 0; nodeIdx < lastPos; nodeIdx += nodeSize ) {
        node = (Node *) (data[blockIdx] + nodeIdx);
        node->adjoint = 0;
    }
}


void tp_clear( Tape *t )
{
    mb_clear( &(t->nodes) );
    mb_clear( &(t->derivs) );
    mb_clear( &(t->arg_ptrs) );
}


typedef struct Tape_Position {
    Block_Position nodes;
    Block_Position derivs;
    Block_Position arg_ptrs;
} Tape_Position;


Tape_Position tp_mark( Tape *t )
{
    Tape_Position tp;
    tp.nodes    = mb_mark( &(t->nodes) );
    tp.derivs   = mb_mark( &(t->derivs) );
    tp.arg_ptrs = mb_mark( &(t->arg_ptrs) );
    return tp;
}


void tp_rewind( Tape *t, Tape_Position tp )
{
    mb_rewind( &(t->nodes),    tp.nodes );
    mb_rewind( &(t->derivs),   tp.derivs );
    mb_rewind( &(t->arg_ptrs), tp.arg_ptrs );
}



// TODO: use TLS
Tape   _TAPE      = { 0 };
const
size_t _NUM_NODES = 1e4;



void _TAPE_init( void )
{
    _TAPE = tp_make( _NUM_NODES, 3*_NUM_NODES );
}

void _TAPE_clear( void )
{
    tp_clear( &_TAPE );
}

void acquire_tape( const Tape *newT )
{
    _TAPE = *newT;
}

void release_tape( Tape *old )
{
    *old = _TAPE;
#if DEBUG
    memset(&_TAPE, 0, sizeof(Tape));
#endif
}

Tape_Position _TAPE_mark( void )
{
    return tp_mark( &_TAPE );
}

void _TAPE_rewind( Tape_Position tp )
{
    tp_rewind( &_TAPE, tp );
}

void _TAPE_deinit( void )
{
    tp_free( &_TAPE );
}

void _TAPE_reset_adjoints( void )
{
    tp_reset_adjoints( &_TAPE );
}


// --------------------- //
// Reverse-mode Variable //
// --------------------- //

// Equivalent to the Number type in
// "Modern Computational Finance - AAD and Parallel Simuations", p. 379


struct Rev {
    double  val;
    Node   *node;


    Rev( void ) {};


    // produces a leaf node on tape
    Rev( double val ) : val( val )
    {
        node = tp_alloc_node( &_TAPE, 0 );
    }


    // unary operation
    Rev( Node *arg, double val ) : val( val )
    {
        node = tp_alloc_node( &_TAPE, 1 );
        node->adj_ptrs[0] = &(arg->adjoint);
    }


    // binary operation
    Rev( Node *lhs, Node *rhs, double val ) : val( val )
    {
        node = tp_alloc_node( &_TAPE, 2 );
        node->adj_ptrs[0] = &(lhs->adjoint);
        node->adj_ptrs[1] = &(rhs->adjoint);
    }


    // ternary operation
    Rev( Node *n1, Node *n2, Node *n3, double val ) : val( val )
    {
        node = tp_alloc_node( &_TAPE, 3 );
        node->adj_ptrs[0] = &(n1->adjoint);
        node->adj_ptrs[1] = &(n2->adjoint);
        node->adj_ptrs[2] = &(n3->adjoint);
    }


    // explicit to avoid silent casting / conversion
    explicit operator double& ()       { return val; }
    explicit operator double  () const { return val; }


    // unary operator overloads

    Rev& operator=( double newVal );

    Rev operator-() const
    {
        return 0.0 - *this;
    }

    Rev operator+() const
    {
        return *this;
    }


    // binary operator overloads

    Rev operator+( Rev rhs );
    Rev operator+( double rhs );
    friend
    Rev operator+( double lhs, Rev rhs );


    Rev operator-( Rev rhs );
    Rev operator-( double rhs );
    friend
    Rev operator-( double lhs, Rev rhs );


    Rev operator*( Rev rhs );
    Rev operator*( double rhs );
    friend
    Rev operator*( double lhs, Rev rhs );


    Rev operator/( Rev rhs );
    Rev operator/( double rhs );
    friend
    Rev operator/( double lhs, Rev rhs );


    Rev& operator+=( Rev rhs );
    Rev& operator+=( double rhs );

    Rev& operator-=( Rev rhs );
    Rev& operator-=( double rhs );

    Rev& operator*=( Rev rhs );
    Rev& operator*=( double rhs );

    Rev& operator/=( Rev rhs );
    Rev& operator/=( double rhs );


    friend bool operator==( Rev lhs,    Rev rhs );
    friend bool operator==( Rev lhs,    double rhs );
    friend bool operator==( double lhs, Rev rhs );

    friend bool operator!=( Rev lhs,    Rev rhs );
    friend bool operator!=( Rev lhs,    double rhs );
    friend bool operator!=( double lhs, Rev rhs );

    friend bool operator<( Rev lhs,    Rev rhs );
    friend bool operator<( Rev lhs,    double rhs );
    friend bool operator<( double lhs, Rev rhs );

    friend bool operator>( Rev lhs,    Rev rhs );
    friend bool operator>( Rev lhs,    double rhs );
    friend bool operator>( double lhs, Rev rhs );

    friend bool operator<=( Rev lhs,    Rev rhs );
    friend bool operator<=( Rev lhs,    double rhs );
    friend bool operator<=( double lhs, Rev rhs );

    friend bool operator>=( Rev lhs,    Rev rhs );
    friend bool operator>=( Rev lhs,    double rhs );
    friend bool operator>=( double lhs, Rev rhs );
};


void rev_print( Rev v, const char *name = NULL) {
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


Rev& Rev::operator=( double newVal )
{
    val  = newVal;
    node = tp_alloc_node( &_TAPE, 0 );
    return *this;
}


inline
Rev rev_make( double val, size_t n )
{
    Rev v;

    v.val  = val;
    v.node = tp_alloc_node( &_TAPE, n );

    return v;
}


inline
Rev rev_make_unary( Node *arg, double val )
{
    Rev v;

    v.val  = val;
    v.node = tp_alloc_node( &_TAPE, 1 );
    v.node->adj_ptrs[0] = &(arg->adjoint);

    return v;
}


inline
Rev rev_make_binary( Node *lhs, Node *rhs, double val )
{
    Rev v;

    v.val  = val;
    v.node = tp_alloc_node( &_TAPE, 2 );
    v.node->adj_ptrs[0] = &(lhs->adjoint);
    v.node->adj_ptrs[1] = &(rhs->adjoint);

    return v;
}


// to produce a leaf node on tape
inline
void rev_put_on_tape( Rev *v )
{
    v->node = tp_alloc_node( &_TAPE, 0 );
}


inline
void rev_put_on_tape_n( Rev **a, size_t n_a )
{
    for ( size_t i = 0; i < n_a; ++i ) {
        rev_put_on_tape( a[i] );
    }
}


// access to adjoint
inline
double rev_adjoint( Rev v )
{
    return v.node->adjoint;
}


// accessing local derivatives with bounds check

inline
double& rev_deriv( Rev v )
{
    assert( v.node->n == 1 );
    return v.node->derivs[0];
}


inline
double& rev_deriv_l( Rev v )
{
    assert( v.node->n == 2 );
    return v.node->derivs[0];
}


inline
double& rev_deriv_r( Rev v )
{
    assert( v.node->n == 2 );
    return v.node->derivs[1];
}


inline
double& rev_deriv_n( Rev v, size_t n )
{
    assert( n < v.node->n );
    return v.node->derivs[n];
}


// accessing child node adjoints with bounds check

inline
double * rev_adj( Rev v )
{
    assert( v.node->n == 1 );
    return v.node->adj_ptrs[0];
}


inline
double * rev_adj_l( Rev v )
{
    assert( v.node->n == 2 );
    return v.node->adj_ptrs[0];
}


inline
double * rev_adj_r( Rev v )
{
    assert( v.node->n == 2 );
    return v.node->adj_ptrs[1];
}


inline
double * rev_adj_n( Rev v, size_t n )
{
    assert( n < v.node->n );
    return v.node->adj_ptrs[0];
}


// core function where back propagation of takes place
// "Modern Computational Finance - AAD and Parallel Simuations", p. 389

void backprop_between( Tape_Position start, Tape_Position end )
{

#define NODE_PROPAGATION \
    node = (Node *) (data[blockIdx] + nodeIdx); \
    nd_propagate_one( *node );

#define BACKPROP_LOOP \
    for ( int64_t nodeIdx = cap - nodeSize; nodeIdx >= 0; nodeIdx -= nodeSize ) { \
        NODE_PROPAGATION \
    }


    int64_t block    = (int64_t)end.nodes.block;
    int64_t position = (int64_t)end.nodes.element;

    int64_t blockIdx = (int64_t)start.nodes.block;
    int64_t lastPos  = (int64_t)start.nodes.element;

    int64_t cap      = (int64_t)_TAPE.nodes.arrayCap;

    if ( cap == 0 )
        return;

    int64_t nodeSize = sizeof(Node);
    int64_t nodeIdx  = 0;

    assert( cap == nodeSize * (int64_t)_TAPE.nodes.numElements );

    Node   *node = NULL;
    unsigned char  **data = _TAPE.nodes.raw.data;

    int64_t firstBlockPos = 0;

    if ( blockIdx == block ) {
        firstBlockPos = position;
    }
    for ( nodeIdx = lastPos - nodeSize; nodeIdx >= firstBlockPos; nodeIdx -= nodeSize ) {
        NODE_PROPAGATION
    }

    blockIdx -= 1;

    for ( ; blockIdx > block; --blockIdx ) {
        BACKPROP_LOOP
    }

    if (blockIdx == block) {
        BACKPROP_LOOP
    }

#undef NODE_PROPAGATION
#undef BACKPROP_LOOP
}


inline
void backprop_until( Rev v, Tape_Position tp )
{
    v.node->adjoint = 1.0;
    backprop_between( _TAPE_mark(), tp );
}

inline
void backprop_to_mark( Rev v, Tape_Position tp ) {
    backprop_until( v, tp );
}

inline
void backprop( Rev v ) {
    backprop_until( v, { 0 } );
}



// Binary operations
// Equivalent to
// "Modern Computational Finance - AAD and Parallel Simuations", p. 393

inline
Rev Rev::operator+( Rev rhs )
{
    double resVal = val + rhs.val;

    Rev res( node, rhs.node, resVal );

    rev_deriv_l( res ) = 1.0;
    rev_deriv_r( res ) = 1.0;

    return res;
}

inline
Rev Rev::operator+( double rhs )
{
    double resVal = val + rhs;

    Rev res( node, resVal );

    rev_deriv( res ) = 1.0;

    return res;
}

inline
Rev operator+( double lhs, Rev rhs )
{
    return rhs + lhs;
}


inline
Rev Rev::operator-( Rev rhs )
{
    double resVal = val - rhs.val;

    Rev res( node, rhs.node, resVal );

    rev_deriv_l( res ) = 1.0;
    rev_deriv_r( res ) = -1.0;

    return res;
}

inline
Rev Rev::operator-( double rhs )
{
    double resVal = val - rhs;

    Rev res( node, resVal );

    rev_deriv( res ) = 1.0;

    return res;
}

inline
Rev operator-( double lhs, Rev rhs )
{
    double resVal = lhs - rhs.val;

    Rev res( rhs.node, resVal );

    rev_deriv( res ) = -1.0;

    return res;
}


inline
Rev Rev::operator*( Rev rhs )
{
    double resVal = val * rhs.val;

    Rev res( node, rhs.node, resVal );

    rev_deriv_l( res ) = rhs.val;
    rev_deriv_r( res ) = val;

    assert( isfinite(rev_deriv_l(res)) );
    assert( isfinite(rev_deriv_r(res)) );

    return res;
}

inline
Rev Rev::operator*( double rhs )
{
    double resVal = val * rhs;

    Rev res( node, resVal );

    rev_deriv( res ) = rhs;

    assert( isfinite(rev_deriv(res)) );

    return res;
}

inline
Rev operator*( double lhs, Rev rhs )
{
    return rhs * lhs;
}


inline
Rev Rev::operator/( Rev rhs )
{
    double resVal = val / rhs.val;

    Rev res( node, rhs.node, resVal );

    double irhs = 1.0 / rhs.val;

    rev_deriv_l( res ) = irhs;
    rev_deriv_r( res ) = -val * irhs * irhs;

    assert( isfinite(rev_deriv_l(res)) );
    assert( isfinite(rev_deriv_r(res)) );

    return res;
}

inline
Rev Rev::operator/( double rhs )
{
    double resVal = val / rhs;

    Rev res( node, resVal );

    rev_deriv( res ) = 1.0 / rhs;

    assert( isfinite(rev_deriv(res)) );

    return res;
}

inline
Rev operator/( double lhs, Rev rhs )
{
    double resVal = lhs / rhs.val;

    Rev res( rhs.node, resVal );

    rev_deriv( res ) = -lhs / (rhs.val * rhs.val);

    assert( isfinite(rev_deriv(res)) );

    return res;
}


Rev& Rev::operator+=( Rev rhs )
{
    *this = *this + rhs;
    return *this;
}

Rev& Rev::operator+=( double rhs )
{
    *this = *this + rhs;
    return *this;
}

Rev& Rev::operator-=( Rev rhs )
{
    *this = *this - rhs;
    return *this;
}

Rev& Rev::operator-=( double rhs )
{
    *this = *this - rhs;
    return *this;
}

Rev& Rev::operator*=( Rev rhs )
{
    *this = *this * rhs;
    return *this;
}

Rev& Rev::operator*=( double rhs )
{
    *this = *this * rhs;
    return *this;
}

Rev& Rev::operator/=( Rev rhs )
{
    *this = *this / rhs;
    return *this;
}

Rev& Rev::operator/=( double rhs )
{
    *this = *this / rhs;
    return *this;
}


bool operator==( Rev lhs, Rev rhs )
{
    return lhs.val == rhs.val;
}

bool operator==( Rev lhs, double rhs )
{
    return lhs.val == rhs;
}

bool operator==( double lhs, Rev rhs )
{
    return lhs == rhs.val;
}


bool operator!=( Rev lhs, Rev rhs )
{
    return lhs.val != rhs.val;
}

bool operator!=( Rev lhs, double rhs )
{
    return lhs.val != rhs;
}

bool operator!=( double lhs, Rev rhs )
{
    return lhs != rhs.val;
}


bool operator<( Rev lhs, Rev rhs )
{
    return lhs.val < rhs.val;
}

bool operator<( Rev lhs, double rhs )
{
    return lhs.val < rhs;
}

bool operator<( double lhs, Rev rhs )
{
    return lhs < rhs.val;
}


bool operator>( Rev lhs, Rev rhs )
{
    return lhs.val > rhs.val;
}

bool operator>( Rev lhs, double rhs )
{
    return lhs.val > rhs;
}

bool operator>( double lhs, Rev rhs )
{
    return lhs > rhs.val;
}


bool operator<=( Rev lhs, Rev rhs )
{
    return lhs.val <= rhs.val;
}

bool operator<=( Rev lhs, double rhs )
{
    return lhs.val <= rhs;
}

bool operator<=( double lhs, Rev rhs )
{
    return lhs <= rhs.val;
}


bool operator>=( Rev lhs, Rev rhs )
{
    return lhs.val >= rhs.val;
}

bool operator>=( Rev lhs, double rhs )
{
    return lhs.val >= rhs;
}

bool operator>=( double lhs, Rev rhs )
{
    return lhs >= rhs.val;
}

