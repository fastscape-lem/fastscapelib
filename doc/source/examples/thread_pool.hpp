#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <atomic>
#include <atomic_queue/atomic_queue.h>
#include <chrono>
#include <functional>
#include <future>
#include <iostream>
#include <thread>

struct Job
{
    std::function<void()> job = nullptr;

    void run()
    {
        if (job != nullptr)
            job();
    }
};

class ThreadPool
{
public:
    explicit ThreadPool(size_t numThreads)
        : stop(false)
    {
        assert(numThreads <= todos.size());

        workers.reserve(numThreads);

        tasks_count = numThreads;
    }

    void set_tasks(std::vector<Job>& jobs_)
    {
        jobs = &jobs_;
        assert(tasks_count == jobs_.size());
    }

    void run_tasks()
    {
        if (!m_started)
            start();

        for (auto i = 0; i < tasks_count; ++i)
            if ((*jobs)[i].job != nullptr)
                todos[i].push(1);
    }

    bool was_empty() const
    {
        for (auto i = 0; i < tasks_count; ++i)
        {
            if (!todos[i].was_empty())
                return false;
        }
        return true;
    }

    void wait() const
    {
        while (!was_empty())
        {
        }
    }

    void close()
    {
        stop = true;

        for (std::thread& worker : workers)
        {
            worker.join();
        }
    }

    std::size_t size() const
    {
        return tasks_count;
    }

    bool stopped() const
    {
        return stop;
    }

    void start()
    {
        m_started = true;

        for (size_t i = 0; i < tasks_count; ++i)
        {
            workers.emplace_back(
                [this, i]
                {
                    while (!stop.load(std::memory_order_relaxed))
                    {
                        if (!todos[i].was_empty())
                        {
                            (*jobs)[i].run();
                            todos[i].pop();
                        }
                    }
                });
        }
    }

    bool started() const
    {
        return m_started;
    }

private:
    std::vector<std::thread> workers;
    std::vector<Job>* jobs;  // Adjust the size as needed
    std::atomic_bool stop;
    std::array<atomic_queue::AtomicQueue<int, 1>, 24> todos;
    std::size_t tasks_count;
    bool m_started = false;
};

#endif  // THREAD_POOL_HPP
