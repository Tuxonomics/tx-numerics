// rng.hpp
//
// Initially from https://github.com/Tuxonomics/tx-prng/tree/master
//

#ifndef TX_RNG_HPP
#define TX_RNG_HPP

#ifdef __ARM_NEON
    #include <arm_neon.h>
#endif


#include <stdarg.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>


namespace tx_rng {

using u8  = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using i8  = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

using f32 = float;
using f64 = double;

using b8  = i8;
using b32 = i32;



inline
u32 rotl(const u32 x, i32 k) {
    return (x << k) | (x >> (32 - k));
}

inline
u64 rotl(const u64 x, i32 k) {
    return (x << k) | (x >> (64 - k));
}



/* Convert u64 to uniform f64. See http://xoshiro.di.unimi.it.
 This will cut the number of possible values in half as the lowest significant
 bit will be set to 0 for all returned values. */
inline
f64 to_f64( u64 x )
{
    const union { u64 i; f64 d; }
    u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
}


inline
f64 to_f64( u32 x )
{
    const union { u32 i; f32 d; }
    u = { .i = UINT32_C(0x7F) << 23 | x >> 9 };
    return (f64) (u.d - 1.0f);
}


#define NEXT_FUNC_U32(name) inline u32 name(void *state)
#define NEXT_FUNC_U64(name) inline u64 name(void *state)
template <typename U>
using next_func = U (void *state);

#define NEXT_N_FUNC_U32(name) inline void name(void *state, u32 *output, u64 n)
#define NEXT_N_FUNC_U64(name) inline void name(void *state, u64 *output, u64 n)
template <typename U>
using next_n_func = void (void *state, U *output, u64 n);

#define JUMP_FUNC(name) void name(void *state)
using jump_func = void (void *state);

#define SEED_FUNC(name) void name(void *state, u64 seed)
using seed_func = void (void *state, u64 seed);


template <typename U>
struct Generator {
    void *state;

    next_func<U>    *next;
    next_n_func<U>  *next_n;
    jump_func       *jump;
    seed_func       *seed;
};


template <typename U>
inline
U next( Generator<U> g )
{
    return g.next( g.state );
}


template <typename U>
inline
void next_n( Generator<U> g, U *output, u64 n )
{
    g.next_n( g.state, output, n );
}


template <typename U>
inline
void jump( Generator<U> g )
{
    g.jump( g.state );
}


template <typename U>
inline
void seed( Generator<U> g, u64 seed )
{
    g.seed( g.state, seed );
}


u64 seed_value( void )
{
    return (u64) time( NULL ); // each second a new seed
}


// The SplitMix64 PRNG.
// https://github.com/svaarala/duktape/blob/master/misc/splitmix64.c
// Seed with any value and use to seed other PRNGs.
struct sm64 {
    u64 s;
};


u64 sm64_next( sm64 *state )
{
    u64
    z = (state->s += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}


/* xorshift1024*, see http://vigna.di.unimi.it/ftp/papers/xorshift.pdf */

struct Xorshift1024s {
    u64 s[16];
    i32 p;
};



NEXT_FUNC_U64(xorshift1024s_next)
{
    Xorshift1024s *x = (Xorshift1024s *) state;

    const u64 s0 = x->s[x->p];
    u64       s1 = x->s[ x->p = (x->p + 1) & 15 ];

    s1 ^= s1 << 31;

    x->s[x->p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);

    return x->s[x->p] * UINT64_C(1181783497276652981);
}


NEXT_N_FUNC_U64(xorshift1024s_nextn)
{
    while ( n-- ) {
        output[n] = xorshift1024s_next(state);
    }
}


// NOTE(jonas): the jump is equivalent to 2^512 calls to next
JUMP_FUNC(xorshift1024s_jump)
{
    Xorshift1024s *x = (Xorshift1024s *) state;

    static const u64 JUMP[] = {
        0x84242f96eca9c41d, 0xa3c65b8776f96855, 0x5b34a39f070b5837,
        0x4489affce4f31a1e, 0x2ffeeb0a48316f40, 0xdc2d9891fe68c022,
        0x3659132bb12fea70, 0xaac17d8efa43cab8, 0xc4cb815590989b13,
        0x5ee975283d71c93b, 0x691548c86c1bd540, 0x7910c41d10a1e6a5,
        0x0b5fc64563b3e2a8, 0x047f7684e9fc949d, 0xb99181f2d8f685ca,
        0x284600e3f30e38c3
    };

    u64 t[16] = { 0 };

    for ( i32 i = 0; i < ( sizeof(JUMP) / sizeof(*JUMP) ); ++i ) {
        for (i32 b = 0; b < 64; ++b) {
            if ( JUMP[i] & 1ULL << b ) {
                for (i32 j = 0; j < 16; ++j ) {
                    t[j] ^= x->s[(j + x->p) & 15];
                }
            }
            xorshift1024s_next( state );
        }
    }

    for ( i32 j = 0; j < 16; j++ ) {
        x->s[(j + x->p) & 15] = t[j];
    }
}


SEED_FUNC(xorshift1024s_seed)
{
    Xorshift1024s *x = (Xorshift1024s *) state;

    sm64 sp = (sm64) { .s = seed };

    for ( u32 i=0; i<16; ++i ) {
        x->s[i] = sm64_next( &sp );
    }

    x->p = 0;
}


Generator<u64> xorshift1024s_init( Xorshift1024s *state, u64 seed )
{
    xorshift1024s_seed( state, seed );

    Generator<u64> g = {
        .state  = state,
        .next   = xorshift1024s_next,
        .next_n = xorshift1024s_nextn,
        .jump   = xorshift1024s_jump,
        .seed   = xorshift1024s_seed,
    };

    return g;
}


// xoshiro256** http://xoshiro.di.unimi.it/xoshiro256starstar.c

struct Xoshiro256ss {
    u64 s[4];
};


NEXT_FUNC_U64(xoshiro256ss_next)
{
    Xoshiro256ss *x = (Xoshiro256ss *) state;

    const u64 res = rotl(x->s[1] * 5, 7) * 9;

    const u64 t = x->s[1] << 17;

    x->s[2] ^= x->s[0];
    x->s[3] ^= x->s[1];
    x->s[1] ^= x->s[2];
    x->s[0] ^= x->s[3];

    x->s[2] ^= t;

    x->s[3]  = rotl(x->s[3], 45);

    return res;
}


NEXT_N_FUNC_U64(xoshiro256ss_nextn)
{
    while ( n-- ) {
        output[n] = xoshiro256ss_next(state);
    }
}


// NOTE(jonas): the jump is equivalent to 2^128 calls to next
JUMP_FUNC(xoshiro256ss_jump)
{
    Xoshiro256ss *x = (Xoshiro256ss *) state;

    static const u64 JUMP[] = {
        0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
        0xa9582618e03fc9aa, 0x39abdc4529b1661c
    };

    u64 s0 = 0;
    u64 s1 = 0;
    u64 s2 = 0;
    u64 s3 = 0;
    for ( i32 i = 0; i < sizeof( JUMP ) / sizeof( *JUMP ); ++i ) {
        for ( i32 b = 0; b < 64; ++b ) {
            if ( JUMP[i] & UINT64_C(1) << b ) {
                s0 ^= x->s[0];
                s1 ^= x->s[1];
                s2 ^= x->s[2];
                s3 ^= x->s[3];
            }
            xoshiro256ss_next( state );
        }
    }

    x->s[0] = s0;
    x->s[1] = s1;
    x->s[2] = s2;
    x->s[3] = s3;
}


SEED_FUNC(xoshiro256ss_seed)
{
    Xoshiro256ss *x = (Xoshiro256ss *) state;

    sm64 sp = (sm64) { .s = seed };

    for ( u32 i=0; i<4; ++i ) {
        x->s[i] = sm64_next( &sp );
    }
}


Generator<u64> xoshiro256ss_init( Xoshiro256ss *state, u64 seed )
{
    xoshiro256ss_seed( state, seed );

    Generator<u64> g = {
        .state  = state,
        .next   = xoshiro256ss_next,
        .next_n = xoshiro256ss_nextn,
        .jump   = xoshiro256ss_jump,
        .seed   = xoshiro256ss_seed,
    };

    return g;
}


// xoshiro128+ https://prng.di.unimi.it/xoshiro128plus.c
struct Xoshiro128p {
    u32 s[4];
};


NEXT_FUNC_U32(xoshiro128p_next)
{
    Xoshiro128p *x = (Xoshiro128p *) state;

    const u32 result = x->s[0] + x->s[3];

    const u32 t = x->s[1] << 9;

// #ifdef __ARM_NEON
//     uint32x2_t s01 = vld1_u32(&x->s[0]);
//     uint32x2_t s23 = vld1_u32(&x->s[2]);

//     s23 = veor_u32(s23, s01);

//     vst1_u32(&x->s[2], s23);

//     x->s[1] ^= x->s[2];
//     x->s[0] ^= x->s[3];
// #else
    x->s[2] ^= x->s[0];
    x->s[3] ^= x->s[1];
    x->s[1] ^= x->s[2];
    x->s[0] ^= x->s[3];
// #endif

    x->s[2] ^= t;

    x->s[3] = rotl(x->s[3], 11);

    return result;
}


NEXT_N_FUNC_U32(xoshiro128p_nextn)
{
    while ( n-- ) {
        output[n] = xoshiro128p_next(state);
    }
}


JUMP_FUNC(xoshiro128p_jump)
{
    Xoshiro128p *x = (Xoshiro128p *) state;

    static const u32 JUMP[] = { 0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b };

    u32 s0 = 0;
    u32 s1 = 0;
    u32 s2 = 0;
    u32 s3 = 0;
    for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
        for(int b = 0; b < 32; b++) {
            if (JUMP[i] & UINT32_C(1) << b) {
                s0 ^= x->s[0];
                s1 ^= x->s[1];
                s2 ^= x->s[2];
                s3 ^= x->s[3];
            }
            xoshiro128p_next( state );
        }

    x->s[0] = s0;
    x->s[1] = s1;
    x->s[2] = s2;
    x->s[3] = s3;
}


SEED_FUNC(xoshiro128p_seed)
{
    Xoshiro128p *x = (Xoshiro128p *) state;

    sm64 sp = (sm64) { .s = seed };

    for ( u32 i=0; i<4; ++i ) {
        x->s[i] = (u32)sm64_next( &sp );
    }
}


Generator<u32> xoshiro128p_init( Xoshiro128p *state, u64 seed )
{
    xoshiro128p_seed( state, seed );

    Generator<u32> g = {
        .state  = state,
        .next   = xoshiro128p_next,
        .next_n = xoshiro128p_nextn,
        .jump   = xoshiro128p_jump,
        .seed   = xoshiro128p_seed,
    };

    return g;
}


#define XOSHIRO128_UNROLL (4)

struct Xoshiro128pv {
    u32 s[XOSHIRO128_UNROLL][4];
};


// scalar output only on first state
NEXT_FUNC_U32(xoshiro128pv_next)
{
    Xoshiro128pv *x = (Xoshiro128pv *) state;

    const u32 result = x->s[0][0] + x->s[0][3];

    const u32 t = x->s[0][1] << 9;

    x->s[0][2] ^= x->s[0][0];
    x->s[0][3] ^= x->s[0][1];
    x->s[0][1] ^= x->s[0][2];
    x->s[0][0] ^= x->s[0][3];

    x->s[0][2] ^= t;

    x->s[0][3] = rotl(x->s[0][3], 11);

    return result;
}


NEXT_N_FUNC_U32(xoshiro128pv_nextn)
{
    Xoshiro128pv *x = (Xoshiro128pv *) state;

    u64 batch_count = n / XOSHIRO128_UNROLL;
    u64 left_count = n - batch_count * XOSHIRO128_UNROLL;

    while ( batch_count-- ) {
#if __ARM_NEON
        output[0] = x->s[0][0] + x->s[0][3];
        output[1] = x->s[1][0] + x->s[1][3];
        output[2] = x->s[2][0] + x->s[2][3];
        output[3] = x->s[3][0] + x->s[3][3];

        const u32 t0 = x->s[0][1] << 9;
        const u32 t1 = x->s[1][1] << 9;
        const u32 t2 = x->s[2][1] << 9;
        const u32 t3 = x->s[3][1] << 9;


        x->s[0][2] ^= x->s[0][0];
        x->s[0][3] ^= x->s[0][1];
        x->s[0][1] ^= x->s[0][2];
        x->s[0][0] ^= x->s[0][3];

        x->s[1][2] ^= x->s[1][0];
        x->s[1][3] ^= x->s[1][1];
        x->s[1][1] ^= x->s[1][2];
        x->s[1][0] ^= x->s[1][3];

        x->s[2][2] ^= x->s[2][0];
        x->s[2][3] ^= x->s[2][1];
        x->s[2][1] ^= x->s[2][2];
        x->s[2][0] ^= x->s[2][3];

        x->s[3][2] ^= x->s[3][0];
        x->s[3][3] ^= x->s[3][1];
        x->s[3][1] ^= x->s[3][2];
        x->s[3][0] ^= x->s[3][3];


        x->s[0][2] ^= t0;
        x->s[1][2] ^= t1;
        x->s[2][2] ^= t2;
        x->s[3][2] ^= t3;

        x->s[0][3] = rotl(x->s[0][3], 11);
        x->s[1][3] = rotl(x->s[1][3], 11);
        x->s[2][3] = rotl(x->s[2][3], 11);
        x->s[3][3] = rotl(x->s[3][3], 11);

#elif __AVX__
        output[0] = x->s[0][0] + x->s[0][3];
        output[1] = x->s[1][0] + x->s[1][3];
        output[2] = x->s[2][0] + x->s[2][3];
        output[3] = x->s[3][0] + x->s[3][3];

        const u32 t0 = x->s[0][1] << 9;
        const u32 t1 = x->s[1][1] << 9;
        const u32 t2 = x->s[2][1] << 9;
        const u32 t3 = x->s[3][1] << 9;


        x->s[0][2] ^= x->s[0][0];
        x->s[0][3] ^= x->s[0][1];
        x->s[0][1] ^= x->s[0][2];
        x->s[0][0] ^= x->s[0][3];

        x->s[1][2] ^= x->s[1][0];
        x->s[1][3] ^= x->s[1][1];
        x->s[1][1] ^= x->s[1][2];
        x->s[1][0] ^= x->s[1][3];

        x->s[2][2] ^= x->s[2][0];
        x->s[2][3] ^= x->s[2][1];
        x->s[2][1] ^= x->s[2][2];
        x->s[2][0] ^= x->s[2][3];

        x->s[3][2] ^= x->s[3][0];
        x->s[3][3] ^= x->s[3][1];
        x->s[3][1] ^= x->s[3][2];
        x->s[3][0] ^= x->s[3][3];


        x->s[0][2] ^= t0;
        x->s[1][2] ^= t1;
        x->s[2][2] ^= t2;
        x->s[3][2] ^= t3;

        x->s[0][3] = rotl(x->s[0][3], 11);
        x->s[1][3] = rotl(x->s[1][3], 11);
        x->s[2][3] = rotl(x->s[2][3], 11);
        x->s[3][3] = rotl(x->s[3][3], 11);

#else
        output[0] = xoshiro128p_next(x->s[0]);
        output[1] = xoshiro128p_next(x->s[1]);
        output[2] = xoshiro128p_next(x->s[2]);
        output[3] = xoshiro128p_next(x->s[3]);
#endif

        output += XOSHIRO128_UNROLL;
    }

    while ( left_count-- ) {
        (output++)[0] = xoshiro128p_next(x->s[0]);
    }
}


JUMP_FUNC(xoshiro128pv_jump)
{
    Xoshiro128pv *x = (Xoshiro128pv *) state;

    for ( size_t i = 0; i < XOSHIRO128_UNROLL; ++i ) {
        xoshiro128p_jump(x->s[i]);
    }
}


SEED_FUNC(xoshiro128pv_seed)
{
    Xoshiro128pv *x = (Xoshiro128pv *) state;

    sm64 sp = (sm64) { .s = seed };

    for ( u32 j = 0; j < 4; ++j ) {
        x->s[0][j] = (u32)sm64_next( &sp );
    }

    for ( u32 i = 1; i < XOSHIRO128_UNROLL; ++i ) {
        for ( u32 j = 0; j < 4; ++j ) {
            x->s[i][j] = x->s[i-1][j];
        }
        xoshiro128p_jump(x->s[i]);
    }
}


Generator<u32> xoshiro128pv_init( Xoshiro128pv *state, u64 seed )
{
    xoshiro128pv_seed( state, seed );

    Generator<u32> g = {
        .state  = state,
        .next   = xoshiro128pv_next,
        .next_n = xoshiro128pv_nextn,
        .jump   = xoshiro128pv_jump,
        .seed   = xoshiro128pv_seed,
    };

    return g;
}


// ------------- //
// Distributions //
// ------------- //


// Uniform U[0,1] random number
template <typename U>
inline
f64 uniform( Generator<U> g )
{
    return to_f64( next( g ) );
}



template <typename U>
inline
void uniform_n( Generator<U> g, f64 *output, u64 n )
{
#define BATCH_N 4

    U buffer[BATCH_N];

    U batch_count = U(n) / BATCH_N;
    U left_count = U(n) - batch_count * BATCH_N;

    while ( batch_count-- ) {
        next_n( g, buffer, BATCH_N );

        output[0] = to_f64(buffer[0]);
        output[1] = to_f64(buffer[1]);
        output[2] = to_f64(buffer[2]);
        output[3] = to_f64(buffer[3]);

        output += BATCH_N;
    }

    while ( left_count-- ) {
        U tmp = next( g );
        (output++)[0] = to_f64( tmp );
    }

#undef BATCH_N
}


template <typename U>
inline
f64 uniform_positive( Generator<U> g )
{
    f64 res;
    do {
        res = uniform( g );
    } while( res == 0 );
    return res;
}


f64 uniform_pdf( f64 a, f64 b )
{
    return 1.0 / ( b - a );
}


f64 uniform_lpdf( f64 a, f64 b )
{
    return - std::log( b - a );
}


f64 uniform_cdf( f64 x, f64 a, f64 b )
{
    if ( x >= a || x <= b ) {
        return ( x - a ) / ( b - a );
    }
    else if ( x < a ) {
        return 0.0;
    }
    else {
        return 1.0;
    }
}


// Standard normal random number with Box-Muller transform.
// from Numerical Recipes in C
// prepared for threadsafety without static variables

template <typename U>
f64 normal_bm( Generator<U> g )
{
    f64 u, v, rsq;

    do {
        u = 2.0 * uniform( g ) - 1.0;
        v = 2.0 * uniform( g ) - 1.0;

        rsq = u*u + v*v;
    } while ( rsq >= 1.0 || rsq == 0.0 );


    return u * std::sqrt( -2.0 * std::log(rsq)/rsq );
}


template <typename U>
inline
f64 normal( Generator<U> g )
{
    return normal_bm( g );
}


f64 normal_pdf( f64 x, f64 m, f64 s )
{
    f64 tmp1 = 2.0 * s * s;
    f64 tmp2 = x - m;

    f64 p1 = 1.0 / std::sqrt( M_PI * tmp1 );
    f64 p2 = std::exp( - tmp2 * tmp2 / tmp1 );

    return p1 * p2;
}


f64 normal_lpdf( f64 x, f64 m, f64 s )
{
    f64 tmp1 = 2.0 * s * s;
    f64 tmp2 = x - m;

    return -0.5 * std::log( M_PI * tmp1 ) - tmp2 * tmp2 / tmp1;
}


f64 normal_cdf( f64 x, f64 m, f64 s )
{
    return 0.5 * std::erfc( -M_SQRT1_2 * (x-m)/s );
}


} // end of namespace tx_rng





#endif
