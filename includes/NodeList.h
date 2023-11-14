#ifndef _NODE_LIST_H
#define _NODE_LIST_H

class NodeList
{
public:
    NodeList(const unsigned int num_threads, const unsigned int size) noexcept
        : num_threads(num_threads),
          r(num_threads, std::vector<NODE_ID>(size, NULL_NODE)),
          r_end(num_threads, 0)
    {}

    // the node list should be empty at this point
    void clear_buffer() const noexcept
    {
        for(unsigned int i = 0; i < num_threads; i++)
            assert(r_end[i] == 0);
    }

    /**
     * Adds node to the end of R[b].
     * @param node
     */
    void add_node(const NODE_ID node) noexcept
    {
        const unsigned int t = static_cast<unsigned int>(omp_get_thread_num());
        if(r_end[t] == r[t].size())
            r[t].resize(r[t].size() << 1);

        r[t][r_end[t]] = node;
        r_end[t]++;
    }

    unsigned int get_r_end(const unsigned int t) const noexcept
    {
        return r_end[t];
    }

    void reset_r_end(const unsigned int t) noexcept
    {
        r_end[t] = 0;
    }

    NODE_ID get_node(const unsigned int t, const unsigned int i) const noexcept
    {
        assert(i < r_end[t]);
        assert(r[t][i] != NULL_NODE);
        return r[t][i];
    }

private:
    const unsigned int num_threads;
    std::vector<std::vector<NODE_ID>> r;
    std::vector<unsigned int> r_end;
};

#endif  // _NODE_LIST_H