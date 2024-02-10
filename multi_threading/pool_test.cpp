// pool_test.cpp

#include "../utilities/tx_tests.hpp"
#include "pool.hpp"

using namespace std;

#define ASSERT TX_ASSERT
#define ASSERT_MSG TX_ASSERT_MSG
#define ASSERT_MSG_VA TX_ASSERT_MSG

#define APPROX TX_APPROX
#define APPROX_ARRAY TX_APPROX_ARRAY


Thread_Pool tp;

std::atomic_int _COUNTER;

void test_thread_pool( void )
{
    int workloads = 1e6;

    tp_init( &tp, 220 );

    for ( int i = 0; i < workloads; ++i ) {
        if ( !tp_push( &tp, [] {
            ++_COUNTER;
        } ) ) {
            printf("couldn't queue workload #%d\n", i);
        }
    }

    tp_main_work( &tp );
    tp_deinit( &tp );

    ASSERT( _COUNTER.load() == workloads );

    printf("counter = %d\n", _COUNTER.load());
}


TX_TEST_LIST = {
    TX_ADD_TEST(test_thread_pool),
};


TX_TESTS_MAIN
