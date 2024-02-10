// pool.hpp
//
// A thread pool that is lock-free dependent on the architecture.
//
// Reference: A. Williams - "C++ Concurrency in Action", 2012
//

#ifndef TX_THREAD_POOL_HPP
#define TX_THREAD_POOL_HPP

#include <stdlib.h>
#include <thread>
#include <functional>
#include <atomic>

#include "third_party/concurrentqueue-1.0.4/concurrentqueue.h"


using TP_Item = std::function<void(void)>;


struct Worker {
    std::thread::id id;
    std::thread     t;
};


struct Thread_Pool {
    moodycamel::ConcurrentQueue<TP_Item> queue;
    std::atomic_bool continue_run;

    Worker          *workers;
    size_t           n_workers;
    std::atomic_int  n_active;

    bool             is_initialized;
};


bool tp_push( Thread_Pool *tp, TP_Item item )
{
    if ( !tp->queue.enqueue(item) ) {
        return false;
    }

    tp->continue_run.store(true);
    tp->continue_run.notify_one();

    return true;
}


void tp_wkr_work( Thread_Pool *tp )
{
    while (true) {
        TP_Item f;
        ++tp->n_active;
        while (tp->queue.try_dequeue(f)) {
            f();
        }

        if ( tp->queue.size_approx() == 0 ) {
            tp->continue_run.store(false);
        }

        --tp->n_active;
        tp->continue_run.wait(false);
    }
}


void tp_main_work( Thread_Pool *tp )
{
    while ( tp->n_active || tp->queue.size_approx() ) {
        TP_Item f;
        while (tp->queue.try_dequeue(f)) {
            f();
        }
    }
}


void tp_deinit( Thread_Pool *tp )
{
    // NOTE: let the OS clean up the threads
    free(tp->workers);
    tp->n_workers = 0;

    tp->is_initialized = false;
}


bool tp_init( Thread_Pool *tp, size_t n_workers = 0 )
{
    if ( tp->is_initialized ) {
        return false;
    }

    tp->n_active.store(0);

    if ( n_workers ) {
        tp->n_workers = n_workers;
    }
    else {
        size_t available_threads = std::thread::hardware_concurrency();
        tp->n_workers = (available_threads == 0) ? 4 : available_threads-1;
    }

    tp->workers = (Worker *) calloc( tp->n_workers, sizeof(*tp->workers) );

    for ( size_t i = 0; i < tp->n_workers; ++i ) {
        Worker *worker = &tp->workers[i];

        tp->workers[i].t = std::thread(
            [tp] {
                tp_wkr_work( tp );
            }
        );
        tp->workers[i].id = tp->workers[i].t.get_id();
    }

    tp->is_initialized = true;

    return true;
}


#endif
