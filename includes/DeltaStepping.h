#ifndef _DELTA_STEPPING_H
#define _DELTA_STEPPING_H

#include <SsspTree.h>
#include <Request.h>
#include <NodeList.h>
#include <Profiling.h>
#include <pthread.h>

template<typename Partition, template<typename> class BucketList, typename RequestList, typename GraphType, template<typename> class KSPGraphType = KSPGraph>
class DeltaStepping
{
protected:
#ifdef PROF
    EDGE_ID touched = 0;
#endif
    const KSPGraphType<GraphType>& g;
    Partition partition;
    BucketList<Partition> bucket_list;
    RequestList request_list;
    NodeList node_list; // todo rename this. It only stores heavy nodes/edges to add/relax them later. so it should be a heavy_node_buffer or something like this
    SsspTree sssp_tree;
#ifdef TECHNIQUE_BOUND
    w_type sub_k_bound = INFINITE_DISTANCE;
#endif

public:
#ifdef TECHNIQUE_BOUND
    void setKBound(w_type sub_k_bound=INFINITE_DISTANCE)
    {
        this->sub_k_bound = sub_k_bound;
    }
#endif

    std::vector<std::pair<NODE_ID, unsigned int>> get_reachable_nodes(unsigned int path_node_num_min = 3, unsigned int path_node_num_max = 8) {
        return sssp_tree.get_reachable_nodes(path_node_num_min,path_node_num_max);
    }

    void reset_num_partition_threads(const unsigned int new_num_threads) noexcept
    {
        partition.reset_num_threads(new_num_threads);
    }

    const SsspTree& get_sssp_tree() const noexcept
    {
        return sssp_tree;
    }

    std::vector<NODE_ID> get_shortest_path_to(const NODE_ID node, const bool get_reverse_path = false) const
    {
        return sssp_tree.get_shortest_path(node, get_reverse_path);
    }

    bool is_connected_to_source(const NODE_ID node) const noexcept
    {
        return sssp_tree.get_distance(node) < INFINITE_DISTANCE;
    }

    w_type get_distance(const NODE_ID node) const noexcept
    {
        return sssp_tree.get_distance(node);
    }

    std::vector<NODE_ID> get_unreachable_nodes(){
        return sssp_tree.get_unreachable_nodes();
    }

    void get_all_distance(w_type *distance_from_src_2_nodes){
        sssp_tree.get_all_distance(distance_from_src_2_nodes);
    }

protected:
#ifdef TECHNIQUE_BOUND
    void assign_requests_light_with_k_bound(const NODE_ID node) noexcept
    {
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
        if (g.get_removed_node()[node]) {
            return;
        }
#endif
        const auto node_end = g.end_light(node);
        for(EDGE_ID neighor = g.begin_light(node); neighor < node_end; ++neighor)
        {
            if(g.get_removed_edge()[neighor])
            {
                continue;
            }
#ifdef PROF
#pragma omp atomic
            touched++;
#endif
            const auto csr = g.get_original_graph().get_csr_graph();
            w_type  weight  = csr->value[neighor];
            NODE_ID node_to = csr->adj[neighor];

#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
            if (g.get_removed_node()[node_to]) {
                continue;
            }
#endif
            
            w_type new_weight = sssp_tree.at(node).tent + weight;

            if(sub_k_bound < new_weight)
            {
                continue;
            }

            if(new_weight < sssp_tree.at(node_to).tent)
            {
                const UpdateRequest curr(node, node_to, new_weight);
                request_list.insert(curr);
            }
        }
    }
    
    void assign_requests_heavy_with_k_bound(const NODE_ID node) noexcept
    {
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
        if (g.get_removed_node()[node]) {
            return;
        }
#endif
        const auto node_end = g.end_heavy(node);
        for(EDGE_ID neighor = g.begin_heavy(node); neighor < node_end; ++neighor)
        {
            if(g.get_removed_edge()[neighor])
            {
                continue;
            }
#ifdef PROF
#pragma omp atomic
            touched++;
#endif
            const auto csr = g.get_original_graph().get_csr_graph();
            w_type  weight  = csr->value[neighor];
            NODE_ID node_to = csr->adj[neighor];

#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
            if (g.get_removed_node()[node_to]) {
                continue;
            }
#endif
            
            w_type new_weight = sssp_tree.at(node).tent + weight;

            if(sub_k_bound < new_weight)
            {
                continue;
            }

            if(new_weight < sssp_tree.at(node_to).tent)
            {
                const UpdateRequest curr(node, node_to, new_weight);
                request_list.insert(curr);
            }
        }
    }
    
    void assign_request_list_with_k_bound() noexcept
    {
        const auto t = static_cast<unsigned int>(omp_get_thread_num());
        const unsigned int R_end = node_list.get_r_end(t);

        for(unsigned int i = 0; i < R_end; i++)
        {
            assign_requests_heavy_with_k_bound(node_list.get_node(t, i));
        }

        node_list.reset_r_end(t);
    }
#endif

    void assign_request_list() noexcept
    {
        const auto t = static_cast<unsigned int>(omp_get_thread_num());
        const unsigned int R_end = node_list.get_r_end(t);

        for(unsigned int i = 0; i < R_end; i++)
        {
            assign_requests_heavy(node_list.get_node(t, i));
        }

        node_list.reset_r_end(t);
    }

    /**
     * Assigns all outgoing edges of node of type edge_type to the request list, if the can possible improve the distance
     * to the target.
     * @param node
     * @param edge_type
     */
    void assign_requests_light(const NODE_ID node) noexcept
    {
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
        if (g.get_removed_node()[node]) {
            return;
        }
#endif
        const auto node_end = g.end_light(node);
        for(EDGE_ID neighor = g.begin_light(node); neighor < node_end; ++neighor)
        {
            if(g.get_removed_edge()[neighor])
            {
                continue;
            }
#ifdef PROF
#pragma omp atomic
            touched++;
#endif
            const auto csr = g.get_original_graph().get_csr_graph();
            w_type  weight  = csr->value[neighor];
            NODE_ID node_to = csr->adj[neighor];

#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
            if (g.get_removed_node()[node_to]) {
                continue;
            }
#endif
            
            w_type new_weight = sssp_tree.at(node).tent + weight;

            if(new_weight < sssp_tree.at(node_to).tent)
            {
                const UpdateRequest curr(node, node_to, new_weight);
                request_list.insert(curr);
            }
        }
    }

    void assign_requests_heavy(const NODE_ID node) noexcept
    {
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
        if (g.get_removed_node()[node]) {
            return;
        }
#endif
        const auto node_end = g.end_heavy(node);
        for(EDGE_ID neighor = g.begin_heavy(node); neighor < node_end; ++neighor)
        {
            if(g.get_removed_edge()[neighor])
            {
                continue;
            }
#ifdef PROF
#pragma omp atomic
            touched++;
#endif
            const auto csr = g.get_original_graph().get_csr_graph();
            w_type  weight  = csr->value[neighor];
            NODE_ID node_to = csr->adj[neighor];
            
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
            if (g.get_removed_node()[node_to]) {
                continue;
            }
#endif
            
            w_type new_weight = sssp_tree.at(node).tent + weight;

            if(new_weight < sssp_tree.at(node_to).tent)
            {
                const UpdateRequest curr(node, node_to, new_weight);
                request_list.insert(curr);
            }
        }
    }

    /**
     * Checks if r improves the distance to r.node and if so, the update is executed.
     */
    void relax(const UpdateRequest& r) noexcept
    {
        if(r.proposed_tent < sssp_tree.at(r.node).tent)
        {
            bucket_list.move(r.node, static_cast<unsigned int>(r.proposed_tent / g.get_delta()) % b_size);
            sssp_tree[r.node].tent = r.proposed_tent;
            sssp_tree[r.node].parent = r.parent_node;
        }
    }

    DeltaStepping(const KSPGraphType<GraphType>& g, Partition p) noexcept
        : g(g),
          partition(std::move(p)),
          bucket_list({partition, g.compute_bucket_size(), g.get_num_nodes()}),
          request_list({partition, g.get_num_nodes()}),
          node_list(partition.get_num_threads(), g.get_num_nodes()),
          sssp_tree(g.get_num_nodes()),
          b_size(g.compute_bucket_size()),
          node_is_removed(g.get_num_nodes(), false)
    {
        omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(partition.get_num_threads());
    }

    virtual ~DeltaStepping() = default;

    virtual void reset() noexcept
    {
        std::fill(node_is_removed.begin(), node_is_removed.end(), false);
        sssp_tree.reset();
    }

    const unsigned int b_size;
    std::vector<bool> node_is_removed;   // each node has a flag indicating
                                         // whether this node has been removed or not
};

#endif  // _DELTA_STEPPING_H
