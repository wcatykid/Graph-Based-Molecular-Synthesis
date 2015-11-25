/*
 *  This file is part of esynth.
 *
 *  esynth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  esynth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with esynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _THREAD_POOL_GUARD
#define _THREAD_POOL_GUARD 1


#include <iostream>
#include <queue>
#include <unistd.h>
#include <pthread.h>


template <class In_Type, class Out_Type> class Thread_Pool;
template <class In_Type, class Out_Type> void *process_data_func(void * This);
template <class In_Type, class Out_Type> void *worker_func(void * This);


template <class In_Type, class Out_Type>
class Thread_Pool
{
  public:
    Thread_Pool(Out_Type (*p)(In_Type)); // function pointer only
    Thread_Pool(int num_threads, Out_Type (*p)(In_Type)); // num threads and function pointer
    ~Thread_Pool(); // destructor
    void print(); // display some debuggin info

    void push(In_Type data); // push
    int in_q_size(); // size of un-processed item q
    int out_q_size(); // size of processed item q
    Out_Type front(); // front
    void pop(); // pop

  private:
    pthread_t _manager; // for manager thread
    pthread_t *threads; // for worker threads
    int num_threads; // total threads
    std::queue<In_Type> in_q; // items to be processed
    std::queue<Out_Type> out_q; // finished results
    pthread_mutex_t lock_in_q, lock_out_q; // mutex locks for q's

    Out_Type (*process)(In_Type); // misc processing function entered by user
    bool processing_data; // flag to stop processing data - threads will all die, currently no way to start again excpet making a new thread pool
    void process_data(); // asynchronously process all incoming data
    void stop_processing(); // finish with in-progress data then stop processing

    //template <class In_Type, class Out_Type> friend void *process_data_func(void * This);
    friend void *process_data_func<In_Type, Out_Type>(void * This); // manager thread
    friend void *worker_func<In_Type, Out_Type>(void * This); // worker thread
};

template <class In_Type, class Out_Type>
Thread_Pool<In_Type, Out_Type>::Thread_Pool(Out_Type (*p)(In_Type))
{
    num_threads = 10;
    threads=new pthread_t[num_threads];

    process=p;

    pthread_mutex_init(&lock_in_q, NULL);
    pthread_mutex_init(&lock_out_q, NULL);

    process_data();
}

template <class In_Type, class Out_Type>
Thread_Pool<In_Type, Out_Type>::Thread_Pool(int num_threads, Out_Type (*p)(In_Type))
{
    this->num_threads = num_threads;
    threads = new pthread_t[num_threads];

    process = p;

    pthread_mutex_init(&lock_in_q, NULL);
    pthread_mutex_init(&lock_out_q, NULL);

    process_data();
}

template <class In_Type, class Out_Type>
Thread_Pool<In_Type, Out_Type>::~Thread_Pool()
{
    stop_processing();
}

//
// display some debugging info
//
template <class In_Type, class Out_Type>
void Thread_Pool<In_Type, Out_Type>::print()
{
    std::cout << "num_threads = " << num_threads << std::endl;
    std::cout << "in_q length = " << in_q_size() << std::endl;
    std::cout << "out_q length = " << out_q_size() << std::endl;
}

template <class In_Type, class Out_Type>
void Thread_Pool<In_Type, Out_Type>::push(In_Type data)
{
    pthread_mutex_lock(&lock_in_q);
    in_q.push(data);
    pthread_mutex_unlock(&lock_in_q);
}

//
// size of unprocessed item queue
//
template <class In_Type, class Out_Type>
int Thread_Pool<In_Type, Out_Type>::in_q_size()
{
    int size;

    pthread_mutex_lock(&lock_in_q);
    size = in_q.size();
    pthread_mutex_unlock(&lock_in_q);

    return size;
}

//
// size of processed item queue
//
template <class In_Type, class Out_Type>
int Thread_Pool<In_Type, Out_Type>::out_q_size()
{
    int size;

    pthread_mutex_lock(&lock_out_q);
    size = out_q.size();
    pthread_mutex_unlock(&lock_out_q);

    return size;
}

template <class In_Type, class Out_Type>
Out_Type Thread_Pool<In_Type, Out_Type>::front()
{
    Out_Type result;

    pthread_mutex_lock(&lock_out_q);
    result = out_q.front();
    pthread_mutex_unlock(&lock_out_q);

    return result;
}

template <class In_Type, class Out_Type>
void Thread_Pool<In_Type, Out_Type>::pop()
{
    pthread_mutex_lock(&lock_out_q);
    out_q.pop();
    pthread_mutex_unlock(&lock_out_q);
}

//
// asynchronously process all incoming data
//
template <class In_Type, class Out_Type>
void Thread_Pool<In_Type, Out_Type>::process_data()
{
    processing_data = true;

    if (~pthread_create(&_manager, NULL, process_data_func<In_Type, Out_Type>, this))
    {
        std::cout << "Thread pool manager created." << std::endl;
    }
    else
    {
        std::cout << "Thread pool manager creation failed." << std::endl;
        processing_data = false;
    }
}

//
// Finish with in-progress data then stop processing
//
template <class In_Type, class Out_Type>
void Thread_Pool<In_Type, Out_Type>::stop_processing()
{
// std::cerr << "Stop Processing signal given." << std::endl;

    processing_data = false;

//std::cerr << "Before Manager join." << std::endl;

    (void) pthread_join(_manager, NULL);

//std::cerr << "After manager join." << std::endl;
}

template <class In_Type, class Out_Type>
void *process_data_func(void *This_void)
{
    // manager thread
    // init
    Thread_Pool<In_Type, Out_Type>* This = (Thread_Pool<In_Type, Out_Type> *)This_void;

    // spawn workers
    for (int x = 0; x < This->num_threads; x++)
    {
        if(~pthread_create(&This->threads[x], NULL, worker_func<In_Type, Out_Type>, This_void))
            std::cout << "worker " << x << " created" << std::endl;
        else
            std::cout << "worker " << x << " creation failed" << std::endl;
    }

    // manage workers
    while(This->processing_data)
    {
        sleep(1); // sleep 1 second
    } // not much to manage atm, sleep() maybe?

    // de-spawn workers
    for (int x = 0; x < This->num_threads; x++)
    {
        (void) pthread_join(This->threads[x], NULL);
        std::cerr << "worker " << x << " removed" << std::endl;
    }

    // exit
    pthread_exit(NULL);
}

template <class In_Type, class Out_Type>
void *worker_func(void *This_void)
{ // worker thread
    // init
    Thread_Pool<In_Type, Out_Type> * This=(Thread_Pool<In_Type, Out_Type> *)This_void;
    In_Type in;
    Out_Type out;

    while(This->processing_data)
    {
        // pop
        pthread_mutex_lock(&This->lock_in_q);
        if (This->in_q.size()>0)
        {
            in=This->in_q.front();
            This->in_q.pop();
            pthread_mutex_unlock(&This->lock_in_q);
        }
        else
        {
            pthread_mutex_unlock(&This->lock_in_q);
            continue;
        }

        // process
        out=This->process(in);

        // push
        pthread_mutex_lock(&This->lock_out_q);
        This->out_q.push(out);
        pthread_mutex_unlock(&This->lock_out_q);
    }

    //exit
    pthread_exit(NULL);
}

// list all possile templates or move implementation into .h file...
// template class Thread_Pool<int, int>;


#endif
