#ifndef FASTSCAPELIB_UTILS_THREAD_POOL_HPP
#define FASTSCAPELIB_UTILS_THREAD_POOL_HPP

#include <atomic_queue/atomic_queue.h>

#include <atomic>
#include <chrono>
#include <functional>
#include <future>
#include <iostream>
#include <thread>
#include <condition_variable>


namespace fastscapelib
{
    struct Job
    {
        Job() = default;

        Job(const std::function<void()>& j)
            : job(j)
        {
        }

        std::function<void()> job = nullptr;

        void run()
        {
            if (job != nullptr)
                job();
        }
    };

    class thread_pool
    {
    public:
        explicit thread_pool(size_t size)
            : m_stopped(false)
            , m_size(size)
        {
            assert(size <= todos.size());
            workers.reserve(size);
            init_pause_jobs();
        }

        ~thread_pool()
        {
            stop();
        }

        void init_pause_jobs()
        {
            pause_jobs.resize(m_size);
            for (std::size_t i = 0; i < m_size; ++i)
                pause_jobs[i] = Job(
                    [this, i]()
                    {
                        std::unique_lock<std::mutex> lk(m_cv_m);
                        ++m_paused_count;
                        m_cv.wait(lk);
                        --m_paused_count;
                    });
        }
        void set_tasks(std::vector<Job>& jobs_)
        {
            jobs = &jobs_;
            assert(m_size == jobs_.size());
        }

        void run_tasks()
        {
            if (!m_started)
                start();

            for (std::size_t i = 0; i < m_size; ++i)
                if ((*jobs)[i].job != nullptr)
                    todos[i].push(1);
        }

        void pause()
        {
            if (!m_paused)
            {
                wait();
                set_tasks(pause_jobs);
                run_tasks();
                m_paused = true;

                while (m_paused_count != m_size)
                {
                }
            }
        }

        void resume()
        {
            if (m_paused)
            {
                m_cv.notify_all();
                m_paused = false;
                wait();
            }
        }

        bool paused() const
        {
            return m_paused;
        }

        bool was_empty() const
        {
            for (std::size_t i = 0; i < m_size; ++i)
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

        void stop()
        {
            if (!m_stopped)
            {
                m_stopped = true;

                if (m_paused)
                    resume();

                for (std::thread& worker : workers)
                    worker.join();
            }
        }

        std::size_t size() const
        {
            return m_size;
        }

        bool stopped() const
        {
            return m_stopped;
        }

        void start()
        {
            if (!m_started)
            {
                m_started = true;

                for (size_t i = 0; i < m_size; ++i)
                {
                    workers.emplace_back(
                        [this, i]
                        {
                            while (!m_stopped.load(std::memory_order_relaxed))
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
        }

        bool started() const
        {
            return m_started;
        }

        void resize(std::size_t size)
        {
            if (size != m_size)
            {
                // std::cout << "resizing" << std::endl;
                m_size = size;
                stop();
                m_stopped = false;
                workers.clear();
                workers.reserve(size);
                init_pause_jobs();
                m_started = false;
            }
        }

        template <typename T>
        class [[nodiscard]] blocks
        {
        public:
            /**
             * @brief Construct a `blocks` object with the given specifications.
             *
             * @param first_index_ The first index in the range.
             * @param index_after_last_ The index after the last index in the range.
             * @param num_blocks_ The desired number of blocks to divide the range into.
             */
            blocks(const T& first_index_, const T& index_after_last_, const std::size_t num_blocks_)
                : m_first_index(first_index_)
                , m_index_after_last(index_after_last_)
                , m_num_blocks(num_blocks_)
            {
                if (m_index_after_last > m_first_index)
                {
                    const std::size_t total_size
                        = static_cast<size_t>(m_index_after_last - m_first_index);
                    if (m_num_blocks > total_size)
                        m_num_blocks = total_size;
                    m_block_size = total_size / m_num_blocks;
                    m_remainder = total_size % m_num_blocks;
                    if (m_block_size == 0)
                    {
                        m_block_size = 1;
                        m_num_blocks = (total_size > 1) ? total_size : 1;
                    }
                }
                else
                {
                    m_num_blocks = 0;
                }
            }

            /**
             * @brief Get the first index of a block.
             *
             * @param block The block number.
             * @return The first index.
             */
            [[nodiscard]] T start(const std::size_t block) const
            {
                return m_first_index + static_cast<T>(block * m_block_size)
                       + static_cast<T>(block < m_remainder ? block : m_remainder);
            }

            /**
             * @brief Get the index after the last index of a block.
             *
             * @param block The block number.
             * @return The index after the last index.
             */
            [[nodiscard]] T end(const std::size_t block) const
            {
                return (block == m_num_blocks - 1) ? m_index_after_last : start(block + 1);
            }

            /**
             * @brief Get the number of blocks. Note that this may be different than the desired
             * number of blocks that was passed to the constructor.
             *
             * @return The number of blocks.
             */
            [[nodiscard]] std::size_t num_blocks() const
            {
                return m_num_blocks;
            }

        private:
            /**
             * @brief The size of each block (except possibly the last block).
             */
            std::size_t m_block_size = 0;

            /**
             * @brief The first index in the range.
             */
            T m_first_index = 0;

            /**
             * @brief The index after the last index in the range.
             */
            T m_index_after_last = 0;

            /**
             * @brief The number of blocks.
             */
            std::size_t m_num_blocks = 0;

            /**
             * @brief The remainder obtained after dividing the total size by the number of blocks.
             */
            std::size_t m_remainder = 0;
        };  // class blocks
    private:
        std::vector<std::thread> workers;
        std::vector<Job>* jobs;
        std::vector<Job> pause_jobs;
        std::atomic_bool m_stopped;
        std::array<atomic_queue::AtomicQueue<int, 1>, 24> todos;
        std::size_t m_size;
        bool m_started = false, m_paused = false;
        std::atomic<std::size_t> m_paused_count = 0;

        std::condition_variable m_cv;
        std::mutex m_cv_m;
    };
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_UTILS_THREAD_POOL_HPP
