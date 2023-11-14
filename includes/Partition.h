#ifndef _PARTITION_H
#define _PARTITION_H

#include <vector>

class Partition
{
protected:
    unsigned int num_threads;

public:
    explicit Partition(const unsigned int num_threads)
        : num_threads(num_threads)
    {}

    inline unsigned int get_num_threads() const
    {
        return num_threads;
    }

    /**
     * Sets the number of threads to no_threads and updates num_threads_changed.
     * @param new_num_threads
     */
    virtual bool reset_num_threads(const unsigned int new_num_threads) noexcept        // todo can this really change?
    {
        if(new_num_threads != num_threads)
        {
            num_threads = new_num_threads;
            return true;
        }

        return false;
    }
};

class PartitionStatic final : public Partition
{
private:
    const unsigned int size;
    std::vector<unsigned int> ind;
    std::vector<unsigned int> rand_nodes;
    /**
     * Assigns nodes to threads in chunks of size floor(#node / #threads).
     * The remaining nodes are assigned one by one
     */
    void assign_nodes_to_threads() noexcept
    {
        const unsigned int chunk = size / num_threads;
        unsigned int start = 0;
        unsigned int end = start + chunk;

        for(unsigned int i = 0; i < num_threads; i++)
        {
            for(unsigned int j = start; j < end; j++)
            {
                ind[rand_nodes[j]] = i;
            }

            start = end;
            end = start + chunk;
        }

        // assign the nodes that have not been assigned yet in a round robin fassion
        for(unsigned int i = start, t_id = 0; i < size; i++, t_id++)
        {
            ind[rand_nodes[i]] = t_id;
        }
    }

public:
    PartitionStatic(const unsigned int num_threads, const unsigned int num_nodes) noexcept
        : Partition(num_threads),
          size(num_nodes),
          ind(num_nodes),        // gets initialized in the constructor with assign_nodes_to_threads()
          rand_nodes(num_nodes)  // gets also initialized in the constructor
    {
        // todo this can be changed to a random value? But every integer in [0,size) needs to bis in rand_nodes, so maybe this array can only be shuffled
        for(unsigned int j = 0; j < num_nodes; j++)
        {
            rand_nodes[j] = j;
        }

        assign_nodes_to_threads();
    }

    unsigned int get_thread_id(const unsigned int node) const noexcept
    {
        return ind[node];
    }

    /**
     * checks if the node is assigned to the thread with the ID tid
     * @param node
     * @param tid
     * @return
     */
    bool assigned_to_thread(unsigned int node, unsigned int tid) const noexcept
    {
        return ind[node] == tid;
    }

    bool reset_num_threads(const unsigned int num_threads) noexcept override
    {
        if(Partition::reset_num_threads(num_threads))
        {
            assign_nodes_to_threads();
            return true;
        }

        return false;
    }
};

class PartitionDynamic final : public Partition
{
public:
    unsigned int get_thread_id(const unsigned int /*node*/) const noexcept     // todo Does this even work? this should be used hidden in a template
    {
        return 0;
    }

    explicit PartitionDynamic(unsigned int num_threads) noexcept
        : Partition(num_threads)
    {}
};

#endif  // _PARTITION_H
