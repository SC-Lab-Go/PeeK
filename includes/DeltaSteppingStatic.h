#ifndef _DELTA_STEPPING_STATIC_H
#define _DELTA_STEPPING_STATIC_H

#include <DeltaStepping.h>
#include <BucketListStatic.h>
#include <RequestListStatic.h>
#include <Profiling.h>
#include <NodeList.h>

template<typename GraphType, template<typename> class KSPGraphType = KSPGraph>
class DeltaSteppingStatic final : public DeltaStepping<PartitionStatic, BucketListStatic, RequestListStatic, GraphType, KSPGraphType>
{
private:
    using parent_type = DeltaStepping<PartitionStatic, BucketListStatic, RequestListStatic, GraphType, KSPGraphType>;
    using parent_type::g;
    using parent_type::bucket_list;
    using parent_type::request_list;
    using parent_type::node_list;
    using parent_type::node_is_removed;
    using parent_type::partition;
    using parent_type::sssp_tree;
#ifdef PROF
    using parent_type::touched;
#endif

public:
    DeltaSteppingStatic(const KSPGraphType<GraphType>& g, const unsigned int num_threads) noexcept
            : parent_type(g, {num_threads, g.get_num_nodes()})
    {}

    /**
     * todo update comment
     * @param t
     * @param source
     * @param dest
     * @param pt
     * @param start
     * @param end
     */
    template <bool stop_early = false>
    void compute(
#ifdef PROF
            ProfilerStatistics& t,
#endif
            const NODE_ID source, const NODE_ID dest, SsspTree::path_direction pt,
            std::vector<NODE_ID>::iterator start, std::vector<NODE_ID>::iterator end, const w_type
#ifndef PROF
            max_distance
#endif
            = INFINITE_DISTANCE
                ) noexcept
    {
        sssp_tree.define_parent_path(start, end, pt);

        {
#ifdef PROF
            scoped_timer scopedTimer(t.time_sssp_compute);
            compute<stop_early>(t, source, dest);
#else
            compute<stop_early>(source, dest, max_distance);
#endif
        }

        {
#ifdef PROF
            t.num_branch_updates++;
            scoped_timer scopedTimer(t.time_branch_updates);
#endif
            sssp_tree.update_branch();
        }
    }

#ifdef PROF

    template<bool stop_early = false>
    void compute(ProfilerStatistics& t, const NODE_ID source, const NODE_ID destination = NULL_NODE, const w_type max_distance = INFINITE_DISTANCE) noexcept
    {
        touched = 1;
        compute<stop_early>(source, destination, max_distance);
        t.nodes_touched_in_ds += touched;
    }

#endif

    template <bool stop_early = false>
    void compute(const NODE_ID source, const NODE_ID destination = NULL_NODE, const w_type max_distance = INFINITE_DISTANCE) noexcept
    {
        assert(!stop_early || destination != NULL_NODE);

        reset();

        if(stop_early)
            sssp_tree.set_destination(destination);

        const unsigned int max_steps = !stop_early || max_distance == INFINITE_DISTANCE
            ? std::numeric_limits<unsigned int>::max()
            : max_distance / g.get_delta() + 1;

        // add source node to bucket[0]
        bucket_list.insert_to_bucket_list(0, source);
        sssp_tree[source].tent = 0;
        sssp_tree[source].parent = source;

        const unsigned int threads = partition.get_num_threads();   // used in omp macro
        unsigned int smallest_bucket = bucket_list.find_smallest_bucket();

        while(smallest_bucket != BucketListStatic<PartitionStatic>::NULL_BUCKET)
        {
            if(stop_early && bucket_list.get_step_counter() > max_steps)
                break;

            bool destination_in_bucket = false;
            node_list.clear_buffer();
            while(!bucket_list.is_empty(smallest_bucket))
            {
#pragma omp parallel num_threads(threads) reduction(||:destination_in_bucket)   // find_path()
                {
#pragma omp barrier
                    // process the smallest bucket of this thread
                    NODE_ID node = bucket_list.pop(smallest_bucket);
                    while(node != NULL_NODE)
                    {
                        if(stop_early)
                            destination_in_bucket = destination_in_bucket || node == destination;

#pragma omp flush(node_is_removed)
                        if(!node_is_removed[node])
                        {
                            node_is_removed[node] = true;
                            node_list.add_node(node);
                        }

#ifdef TECHNIQUE_BOUND
                        parent_type::assign_requests_light_with_k_bound(node);
#else
                        parent_type::assign_requests_light(node);
#endif
                        node = bucket_list.pop(smallest_bucket);
                    }
#pragma omp barrier
                    relax_requests();
#pragma omp barrier
                }
            }
#pragma omp parallel num_threads(threads)   // find_path()
            {
                request_list.reset_buffer();
#pragma omp barrier

#ifdef TECHNIQUE_BOUND
                parent_type::assign_request_list_with_k_bound();
#else
                parent_type::assign_request_list();
#endif

#pragma omp barrier
                relax_requests();
            }

            if(stop_early && destination_in_bucket)
                break;

            smallest_bucket = bucket_list.find_smallest_bucket(smallest_bucket);
        }
    }

    void compute_sssp(const NODE_ID source) noexcept
    {
        reset();

        // add source node to bucket[0]
        bucket_list.insert_to_bucket_list(0, source);
        sssp_tree[source].tent = 0;
        sssp_tree[source].parent = source;

        const unsigned int threads = partition.get_num_threads();   // used in omp macro
        unsigned int smallest_bucket = bucket_list.find_smallest_bucket();

        while(smallest_bucket != BucketListStatic<PartitionStatic>::NULL_BUCKET)
        {
            node_list.clear_buffer();
            while(!bucket_list.is_empty(smallest_bucket))
            {
#pragma omp parallel num_threads(threads)   // find_path()
                {
#pragma omp barrier
                    // process the smallest bucket of this thread
                    NODE_ID node = bucket_list.pop(smallest_bucket);
                    while(node != NULL_NODE)
                    {
#pragma omp flush(node_is_removed)
                        if(!node_is_removed[node])
                        {
                            node_is_removed[node] = true;
                            node_list.add_node(node);
                        }

                        parent_type::assign_requests_light(node);
                        node = bucket_list.pop(smallest_bucket);
                    }
#pragma omp barrier
                    relax_requests();
#pragma omp barrier
                }
            }
#pragma omp parallel num_threads(threads)   // find_path()
            {
                request_list.reset_buffer();
#pragma omp barrier
                parent_type::assign_request_list();
#pragma omp barrier
                relax_requests();
            }

            smallest_bucket = bucket_list.find_smallest_bucket(smallest_bucket);
        }
    }

private:
    /**
     * This function can be called in parallel. Each thread relaxes its own requests.
     */
    void relax_requests() noexcept
    {
        const auto thread_no = static_cast<unsigned int>(omp_get_thread_num());
        const unsigned int num_threads = request_list.get_num_threads();

        for(unsigned int i = 0; i < num_threads; i++)
        {
            unsigned int buffer_size = request_list.pop_buffer_size(i, thread_no);
            for(unsigned int j = 0; j < buffer_size; j++)
            {
                parent_type::relax(request_list.get_request(i, thread_no, j));
            }
        }
    }

    void reset() noexcept
    {
        parent_type::reset();
        bucket_list.reset();
        request_list.reset();

        std::fill(node_is_removed.begin(), node_is_removed.end(), false);
    }
};

#endif // _DELTA_STEPPING_STATIC_H
