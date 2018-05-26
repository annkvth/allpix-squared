/**
 * @file
 * @brief Definition of thread pool used for module multithreading
 * @copyright MIT License
 */

#include "ThreadPool.hpp"

#include "Module.hpp"

using namespace allpix;

/**
 * The threads are created in an exception-safe way and all of them will be destroyed when creation of one fails
 */
ThreadPool::ThreadPool(unsigned int num_threads, std::function<void()> worker_init_function) {
    // Create threads
    try {
        for(unsigned int i = 0u; i < num_threads; ++i) {
            threads_.emplace_back(&ThreadPool::worker, this, worker_init_function);
        }
    } catch(...) {
        destroy();
        throw;
    }

    run_cnt_ = 0;
}

ThreadPool::~ThreadPool() {
    destroy();
}

void ThreadPool::submit_event_function(std::function<void()> event_function) {
    event_queue_.push(std::make_unique<std::packaged_task<void()>>(std::move(event_function)));
}

/**
 * Run by the \ref ModuleManager to ensure all tasks and modules are completed before moving to the next instantiations.
 * Besides waiting for the queue to empty this will also wait for all the tasks to be completed. If an exception is
 * thrown by another thread, the exception will be propagated to the main thread by this function.
 */
bool ThreadPool::execute_all() {
    while (true) {
        Task task{nullptr};

        if (event_queue_.pop(task, false)) {
            try {
                // Execute task
                (*task)();
                // Fetch the future to propagate exceptions
                task->get_future().get();
            } catch (...) {
                // Check if the first exception thrown
                if(has_exception_.test_and_set()) {
                    // Check if the first exception thrown
                    exception_ptr_ = std::current_exception();
                    // Invalidate the queue to terminate other threads
                    event_queue_.invalidate();
                }
            }
        }

        // Wait for the threads to complete their task, continue helping if a new task was pushed
        std::unique_lock<std::mutex> lock{run_mutex_};
        run_condition_.wait(lock, [this]() { return !event_queue_.empty() || run_cnt_ == 0; });

        // Only stop when both the queue is empty and the run count is zero
        if (event_queue_.empty() && run_cnt_ == 0) {
            break;
        }
    }

    // If exception has been thrown, destroy pool and propagate it
    if (exception_ptr_) {
        destroy();
        Log::setSection("");
        std::rethrow_exception(exception_ptr_);
    }

    return event_queue_.isValid();
}

/**
 * If an exception is thrown by a module, the first exception is saved to propagate in the main thread
 */
void ThreadPool::worker(const std::function<void()>& init_function) {
    // Initialize the worker
    init_function();

    // Safe lambda to increase the atomic run count
    auto increase_run_cnt_func = [this]() { ++run_cnt_; };

    // Continue running until the thread pool is finished
    while (!done_) {
        Task task{nullptr};

        if (event_queue_.pop(task, true, increase_run_cnt_func)) {
            // Try to run the task
            try {
                // Execute task
                (*task)();
                // Fetch the future to propagate exceptions
                task->get_future().get();
            } catch (...) {
                // Check if the first exception thrown
                if (has_exception_.test_and_set()) {
                    // Save first exception
                    exception_ptr_ = std::current_exception();
                    // Invalidate the queue to terminate other threads
                    event_queue_.invalidate();
                }
            }
        }

        // Propagate that the task has been finished
        std::lock_guard<std::mutex> lock{run_mutex_};
        --run_cnt_;
        run_condition_.notify_all();
    }
}

void ThreadPool::destroy() {
    done_ = true;
    event_queue_.invalidate();

    for(auto& thread : threads_) {
        if(thread.joinable()) {
            thread.join();
        }
    }
}
