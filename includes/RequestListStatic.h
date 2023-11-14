#ifndef _REQUEST_LIST_STATIC_H
#define _REQUEST_LIST_STATIC_H

#include <Request.h>
#include <Partition.h>
#include <cassert>
#include <algorithm>

class RequestListStatic
{
public:

    RequestListStatic(const PartitionStatic& p, const NODE_ID num_entries) noexcept
        : partition(p),
          num_threads(p.get_num_threads()),
          buffer(num_threads, std::vector<std::vector<UpdateRequest>>(num_threads, std::vector<UpdateRequest>(std::max(static_cast<unsigned int>(1), (8 * (static_cast<unsigned int>(sqrt(num_entries) + 1))) / num_threads)))),
          buffer_pos(num_threads, std::vector<unsigned int>(num_threads, 0))
    {}

    void reset() noexcept
    {
        for(auto& buffer_pos_i : buffer_pos)
            std::fill(buffer_pos_i.begin(), buffer_pos_i.end(), 0);
    }

    void insert(const UpdateRequest& curr) noexcept
    {
        const unsigned int thread_v = partition.get_thread_id(curr.parent_node);
        const unsigned int thread_w = partition.get_thread_id(curr.node);
        assert(static_cast<unsigned int>(omp_get_thread_num()) == thread_v);

        if(buffer_pos[thread_v][thread_w] == buffer[thread_v][thread_w].size())
            extend_buffer_size(thread_v, thread_w);

        buffer[thread_v][thread_w][buffer_pos[thread_v][thread_w]] = curr;
        buffer_pos[thread_v][thread_w]++;
    }

    unsigned int get_num_threads() const noexcept
    {
        return num_threads;
    }

    const UpdateRequest& get_request(const unsigned int i, const unsigned int thread_no, const unsigned int j) const noexcept
    {
        return buffer[i][thread_no][j];
    }

    /**
     * Returns the size of buffer[i][thread_no] and sets it to zero after that. The requests are not gone as long as the
     * returned value is not gone.
     * @param i
     * @param thread_no
     * @return the size of buffer[i][thread_no] and sets it to zero after that.
     */
    unsigned int pop_buffer_size(const unsigned int i, const unsigned int thread_no) noexcept
    {
        const auto tmp = buffer_pos[i][thread_no];
        buffer_pos[i][thread_no] = 0;
        return tmp;
    }

    /**
     * This will not reset anything. The requests should all be relaxed by now.
     */
    void reset_buffer() noexcept
    {
        const auto b = static_cast<unsigned int>(omp_get_thread_num());
        std::fill(buffer_pos[b].begin(), buffer_pos[b].end(), 0);
    }

private:

    inline void extend_buffer_size(const unsigned int i, const unsigned int j) noexcept
    {
        buffer[i][j].resize(buffer[i][j].size() << 1);
    }

    const PartitionStatic& partition;
    const unsigned int num_threads;
    std::vector<std::vector<std::vector<UpdateRequest>>> buffer;
    std::vector<std::vector<unsigned int>> buffer_pos;
};

#endif  // _REQUEST_LIST_STATIC_H
