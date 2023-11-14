#ifndef _BUCKET_LIST_STATIC_H
#define _BUCKET_LIST_STATIC_H

#include <Partition.h>
#include <KSPGraph.h>
#include <cassert>

template<typename Partition>
class BucketListStatic
{
    unsigned int step_counter = 0;
public:
    static constexpr unsigned int NULL_BUCKET = std::numeric_limits<unsigned int>::max();

    BucketListStatic(const Partition& p, const unsigned int num_buckets, const NODE_ID num_entries) noexcept
        : partition(p),
          num_threads(p.get_num_threads()),
          num_entries(num_entries),
          num_buckets(num_buckets),
          buckets(num_buckets, std::vector<BucketListStatic::bucket_list>(num_threads)),
          bucket_assignment(num_entries),
          old_bucket(num_entries)
    {
        reset();
    }

    inline unsigned int get_bucket(const NODE_ID node) const noexcept
    {
        return old_bucket[node];
    }

    unsigned int get_step_counter() const noexcept
    {
        return step_counter;
    }

    /**
     * Checks if the bucket is empty for all threads
     * @param bucket_id
     * @return
     */
    bool is_empty(const unsigned int bucket_id) const noexcept
    {
        for(unsigned int t = 0; t < num_threads; t++)
        {
            if(!buckets[bucket_id][t].is_empty())
                return false;
        }

        return true;
    }

    void reset() noexcept
    {
        step_counter = 0;
        std::fill(old_bucket.begin(), old_bucket.end(), NULL_BUCKET);

#pragma omp parallel for num_threads(num_threads)   // reset()
        for(NODE_ID j = 0; j < num_entries; j++)
        {
            bucket_assignment[j].clear();
        }

#pragma omp parallel for num_threads(num_threads)   // reset()
        for(unsigned int b = 0; b < num_buckets; b++)
        {
            for(unsigned int t = 0; t < num_threads; t++)
                buckets[b][t].clear();
        }
    }

    /**
     * Moves a node from its old bucked to the new one.
     * @param node
     * @param new_bucket
     */
    void move(const NODE_ID node, const unsigned int new_bucket) noexcept
    {
        const auto bucket_id = get_bucket(node);
        if(bucket_id != new_bucket)
        {
            if(bucket_id != NULL_BUCKET)
            {
                remove_from_bucket_list(partition.get_thread_id(node), bucket_id, node);
            }

            insert_to_bucket_list(new_bucket, node);
        }
    }

    unsigned int find_smallest_bucket(const unsigned int last_bucket = NULL_BUCKET) noexcept
    {
        step_counter++;
        // start = last_bucket + 1 if this exists else 0;
        const unsigned int start = (last_bucket + 1) * static_cast<unsigned int>(last_bucket < num_buckets);

        for(unsigned int b = start; b < num_buckets; b++, step_counter++)
        {
            for(unsigned int t = 0; t < num_threads; t++)
            {
                if(!buckets[b][t].is_empty())
                    return b;
            }
        }

        for(unsigned int b = 0; b < start; b++, step_counter++)
        {
            for(unsigned int t = 0; t < num_threads; t++)
            {
                if(!buckets[b][t].is_empty())
                    return b;
            }
        }

        return NULL_BUCKET;
    }

    void insert_to_bucket_list(const unsigned int bucket_id, const NODE_ID node) noexcept
    {
        const unsigned int thread = partition.get_thread_id(node);

        if(buckets[bucket_id][thread].is_empty())
        {
            buckets[bucket_id][thread].start_node = node;
        }
        else
        {
            bucket_assignment[buckets[bucket_id][thread].start_node].prev_node = node;
            bucket_assignment[node].next_node = buckets[bucket_id][thread].start_node;
            buckets[bucket_id][thread].start_node = node;
        }

        set_bucket(node, bucket_id);
    }

    NODE_ID pop(const unsigned int bucket_id) noexcept
    {
        const auto thread_id = static_cast<unsigned int>(omp_get_thread_num());
        const NODE_ID node = buckets[bucket_id][thread_id].start_node;

        if(node != NULL_NODE)
        {
            buckets[bucket_id][thread_id].start_node = bucket_assignment[node].next_node;
            bucket_assignment[node].clear();
            old_bucket[node] = NULL_BUCKET;
        }

        return node;
    }

    /**
     * Counts (!) the number of nodes in the current bucket over all threads.
     * @param bucket_id
     * @return
     */
    unsigned int length(const unsigned int bucket_id) const noexcept
    {
        unsigned int len = 0;

        for(unsigned int t = 0; t < num_threads; t++)
        {
            NODE_ID node = buckets[bucket_id][t].start_node;

            while(node != NULL_NODE)
            {
                len++;
                node = bucket_assignment[node].next_node;
            }
        }

        return len;
    }

    bool contains_entry(const NODE_ID node) const noexcept
    {
        return node != NULL_NODE && old_bucket[node] != NULL_NODE;
    }

    bool local_is_empty(const unsigned int bucket_id) const noexcept
    {
        return buckets[bucket_id][omp_get_thread_num()].is_empty();
    }

    template<typename P>
    friend std::ostream &operator<<(std::ostream& o, const BucketListStatic<P>& b);

private:
    // structures required by BucketListStatic
    struct bucket_node
    {
        NODE_ID prev_node = NULL_NODE;
        NODE_ID next_node = NULL_NODE;

        bucket_node() = default;
        bucket_node(NODE_ID prev_node, NODE_ID next_node) noexcept  : prev_node(prev_node), next_node(next_node) {}

        void clear() noexcept
        {
            prev_node = NULL_NODE;
            next_node = NULL_NODE;
        }

        bool is_first() const noexcept
        {
            return prev_node == NULL_NODE;
        }

        bool is_last() const noexcept
        {
            return next_node == NULL_NODE;
        }

        bool assigned() const noexcept
        {
            return prev_node != NULL_NODE || next_node != NULL_NODE;
        }
    };

    struct bucket_list
    {
        NODE_ID start_node = NULL_NODE;

        void clear() noexcept
        {
            start_node = NULL_NODE;
        }

        bool is_empty() const noexcept
        {
            return start_node == NULL_NODE;
        }
    };

    inline void set_bucket(const NODE_ID node, const unsigned int b) noexcept
    {
        old_bucket[node] = b;
    }

    // Private Interface Functions
    void remove_from_bucket_list(const unsigned int thread, const unsigned int bucket, const NODE_ID node) noexcept
    {
        const bucket_node& n = bucket_assignment[node];
        if(n.next_node != NULL_NODE)
        {
            assert(partition.get_thread_id(n.next_node) == thread);
        }

        if(n.prev_node != NULL_NODE)
        {
            assert(partition.get_thread_id(n.prev_node) == thread);
        }

        if(n.is_first() && n.is_last())
        {
            buckets[bucket][thread].start_node = NULL_NODE;
        }
        else if(n.is_first())
        {
            buckets[bucket][thread].start_node = n.next_node;
            bucket_assignment[n.next_node].prev_node = NULL_NODE;
        }
        else if(n.is_last())
        {
            bucket_assignment[n.prev_node].next_node = NULL_NODE;
        }
        else
        {
            bucket_assignment[n.prev_node].next_node = n.next_node;
            bucket_assignment[n.next_node].prev_node = n.prev_node;
        }

        bucket_assignment[node].clear();
    }

    // Private Fields
    const Partition& partition;
    const unsigned int num_threads;
    const NODE_ID num_entries;
    const unsigned int num_buckets;

    std::vector<std::vector<bucket_list>> buckets;
    std::vector<bucket_node> bucket_assignment;
    std::vector<unsigned int> old_bucket;
};

template<typename Partition>
constexpr unsigned int BucketListStatic<Partition>::NULL_BUCKET;

template<typename Partition>
std::ostream &operator<<(std::ostream& o, const BucketListStatic<Partition>& b)
{
    for(unsigned int j = 0; j < b.num_buckets; j++)
    {
        o << "(";
        for(unsigned int i = 0; i < b.num_threads; i++)
        {
            NODE_ID node = b.buckets[j][i].start_node;
            while(node != NULL_NODE)
            {
                o << node;
                if(b.bucket_assignment[node].next_node != NULL_NODE)
                    o << ",";
                node = b.bucket_assignment[node].next_node;
            }
        }
        o << ")";
    }
    return o;
}

#endif  // _BUCKET_LIST_STATIC_H
