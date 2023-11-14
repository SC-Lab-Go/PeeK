#ifndef _KSP_GRAPH_H
#define _KSP_GRAPH_H

#include <BasicGraph.h>
#include <vector>
#include <queue>

template<typename GraphType>
class KSPGraph
{
public:
    const GraphType& g;
    std::queue<EDGE_ID> removed;
    std::vector<bool> edge_is_removed;
    
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
    std::vector<bool> node_is_removed;
#endif

public:
    explicit KSPGraph(const GraphType& graph) noexcept
        : g(graph),
          edge_is_removed(g.get_num_edges(), false)
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
        , node_is_removed(g.get_num_nodes(), false)
#endif
    {}

    NODE_ID get_num_nodes() const noexcept { return g.get_num_nodes(); }

    EDGE_ID get_num_edges() const noexcept { return g.get_num_edges(); }

    EDGE_ID begin(const NODE_ID node) const noexcept
    {
        return g.get_edge_list_begin(node);
    }

    EDGE_ID begin_light(const NODE_ID node) const noexcept
    {
        return g.get_edge_list_begin(node);
    }

    EDGE_ID begin_heavy(const NODE_ID node) const noexcept
    {
        return g.get_heavy_edge_list_begin(node);
    }

    EDGE_ID end(const NODE_ID node) const noexcept
    {
        return g.get_edge_list_end(node);
    }

    EDGE_ID end_light(const NODE_ID node) const noexcept
    {
        return g.get_light_edge_list_end(node);
    }

    EDGE_ID end_heavy(const NODE_ID node) const noexcept
    {
        return end(node);
    }

    EDGE_ID get_edge_id(const NODE_ID src, const NODE_ID target) const noexcept
    {
        EDGE_ID end_edge = end(src);
        for(auto neighor = begin(src); neighor < end_edge; ++neighor)
        {
            const auto csr = g.get_csr_graph();
            if(csr->adj[neighor] == target)
                return neighor;
        }

        return end_edge;
    }

    void remove_node(const NODE_ID node) noexcept
    {
        const EDGE_ID edge_list_end = g.get_edge_list_end(node);
        for(EDGE_ID e_id = g.get_edge_list_begin(node); e_id < edge_list_end; e_id++)
        {
            if(edge_is_removed[e_id])
                continue;

            removed.push(e_id);

            edge_is_removed[e_id] = true;
        }
    }

    void remove_edge(const EDGE_ID e_id) noexcept
    {
        if(edge_is_removed[e_id])
            return;

        removed.push(e_id);        // todo can we avoid this?

        edge_is_removed[e_id] = true;

    }

    /**
     * Restores every removed edge.
     */
    void restore_graph() noexcept
    {
        while(!removed.empty())
        {
            edge_is_removed[removed.front()] = false;
            removed.pop();
        }
    }

    const GraphType& get_original_graph() const noexcept
    {
        return g;
    }

    w_type get_delta() const noexcept { return g.get_delta(); }

    double get_avg_degree() const noexcept { return g.get_avg_degree(); }

    NODE_ID get_num_neighbors(const NODE_ID node) const noexcept { return g.get_num_neighbors(node); }

    unsigned int compute_bucket_size() const noexcept { return g.compute_bucket_size(); }

    void print_graph() const noexcept
    {
        const auto csr = g.get_csr_graph();
        for(NODE_ID i = 0; i < get_num_nodes(); i++)
        {
            std::cout << "node:" << i << "[";
            for(EDGE_ID neighor = begin(i); neighor < end(i); ++neighor)
            {
                if(!edge_is_removed[neighor])
                {
                    w_type  weight  = csr->value[neighor];
                    NODE_ID node_to = csr->adj[neighor];
                    std::cout << "("<<node_to<<","<<weight << ") ";
                }
            }
            std::cout << "]" << std::endl;
        }
    }
    
    const std::vector<bool>& get_removed_edge() const noexcept
    {
        return edge_is_removed;
    }

#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
    const std::vector<bool>& get_removed_node() const noexcept { return node_is_removed; }

    void remove_node_permanently(const NODE_ID node) noexcept { node_is_removed[node] = true; }
#endif

#ifdef ADAPTIVE_GRAPH_COMPACT_EDGE
    void remove_edge_permanently(const NODE_ID node) noexcept {
        const EDGE_ID edge_list_end = g.get_edge_list_end(node);
        for (EDGE_ID e_id = g.get_edge_list_begin(node); e_id < edge_list_end; e_id++) {
            edge_is_removed[e_id] = true;
        }
    }
#endif

    double get_remain_edge_ratio() {
        EDGE_ID remain_edge_num = 0;
        EDGE_ID num_edges = g.get_num_edges();
        for (auto elem : edge_is_removed) {
            if (!elem) {
                remain_edge_num++;
            }
        }

        return static_cast<double>(remain_edge_num) / num_edges;
    }
};

template<typename GraphType>
class KSPGraphEdgeSwap
{
private:
    const GraphType& g;
    std::queue<EDGE_ID> removed;
    std::vector<bool> edge_is_removed;

public:
    explicit KSPGraphEdgeSwap(const GraphType& graph) noexcept
        : g(graph),
          edge_is_removed(g.get_num_edges(), false)
    {}

    NODE_ID get_num_nodes() const noexcept { return g.get_num_nodes(); }

    EDGE_ID get_num_edges() const noexcept { return g.get_num_edges(); }

    EDGE_ID begin(const NODE_ID node) const noexcept
    {
        return g.get_edge_list_begin(node);
    }

    EDGE_ID begin_light(const NODE_ID node) const noexcept
    {
        return g.get_edge_list_begin(node);
    }

    EDGE_ID begin_heavy(const NODE_ID node) const noexcept
    {
        return g.get_heavy_edge_list_begin(node);
    }

    EDGE_ID end(const NODE_ID node) const noexcept
    {
        return g.get_edge_list_end(node);
    }

    EDGE_ID end_light(const NODE_ID node) const noexcept
    {
        return g.get_light_edge_list_end_edge_swap(node);
    }

    EDGE_ID end_heavy(const NODE_ID node) const noexcept
    {
        return g.get_heavy_edge_list_end_edge_swap(node);
    }

    void remove_node(const NODE_ID node) noexcept
    {
        auto node_end = end_light(node);
        for(EDGE_ID neighor = begin_light(node); neighor < node_end; neighor++)
        {
            remove_edge(neighor);
        }

        node_end = end_heavy(node);
        for(EDGE_ID neighor = begin_heavy(node); neighor < node_end; neighor++)
        {
            remove_edge(neighor);
        }
    }

    void remove_edge(const EDGE_ID e_id) noexcept
    {
        if(edge_is_removed[e_id])
            return;

        removed.push(e_id);        // todo can we avoid this?

        edge_is_removed[e_id] = true;

    }

    void remove_edge(const NODE_ID src, const NODE_ID dest) noexcept
    {
        const auto csr = g.get_csr_graph();
        auto node_end = end_light(src);
        for(EDGE_ID neighor = begin_light(src); neighor < node_end; neighor++)
        {
            if(csr->adj[neighor] == dest)
            {
                remove_edge(neighor);
                return;
            }
        }

        node_end = end_heavy(src);
        for(EDGE_ID neighor = begin_heavy(src); neighor < node_end; neighor++)
        {
            if(csr->adj[neighor] == dest)
            {
                remove_edge(neighor);
                return;
            }
        }
    }

    /**
     * Restores every removed edge.
     */
    void restore_graph() noexcept
    {
        while(!removed.empty())
        {
            edge_is_removed[removed.front()] = false;
            removed.pop();
        }
    }

    const GraphType& get_original_graph() const noexcept
    {
        return g;
    }

    w_type get_delta() const noexcept { return g.get_delta(); }

    double get_avg_degree() const noexcept { return g.get_avg_degree(); }

    NODE_ID get_num_neighbors(const NODE_ID node) const noexcept { return g.get_num_neighbors(node); }

    unsigned int compute_bucket_size() const noexcept { return g.compute_bucket_size(); }

    const std::vector<bool>& get_removed_edge() const noexcept
    {
        return edge_is_removed;
    }

    //backward
    void init_edge_swap(std::unique_ptr<KSPGraph<GraphType>>& reverse_graph, bool* is_k_bound_nodes, const unsigned int num_threads) noexcept {
        const auto csr = g.get_csr_graph();
        const auto csr_reverse = reverse_graph->get_original_graph().get_csr_graph();
        bool* affected_nodes = new bool[csr->num_nodes]();

        // set memory for edge swap in csr
        csr->initEdgeSwap();

        // set initial value
        const NODE_ID num_nodes = get_num_nodes();
#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID i = 0; i < num_nodes; i++) {
            g.restore_node(i);
        }

#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID node = 0; node < csr->num_nodes; node++) {
            if (!is_k_bound_nodes[node]) {
                // delete nodes
                g.delete_node_by_edge_swap(node);

                // delete edges
                for (EDGE_ID neighor_reverse = csr_reverse->begin[node]; neighor_reverse < csr_reverse->begin[node + 1]; neighor_reverse++) {
                    NODE_ID delete_node_from = csr_reverse->adj[neighor_reverse];

                    if (!affected_nodes[delete_node_from]) {
                        affected_nodes[delete_node_from] = true;
                    }
                }
            }
        }

#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID node = 0; node < csr->num_nodes; node++) {
            if (!affected_nodes[node] || !is_k_bound_nodes[node]) {
                continue;
            }

            for (EDGE_ID neighor = begin_light(node); neighor < end_light(node); ) {
                NODE_ID node_to = csr->adj[neighor];
                if (!is_k_bound_nodes[node_to]) {
                    g.delete_edge_by_edge_swap(node, neighor, true);
                }else{
                    neighor++;
                }
            }

            for (EDGE_ID neighor = begin_heavy(node); neighor < end_heavy(node); ) {
                NODE_ID node_to = csr->adj[neighor];
                if (!is_k_bound_nodes[node_to]) {
                    g.delete_edge_by_edge_swap(node, neighor, false);
                }else{
                    neighor++;
                }
            }
        }

        delete[] affected_nodes;
    }
    //forward
    void init_edge_swap(bool* is_k_bound_nodes, const unsigned int num_threads) noexcept {
        const auto csr = g.get_csr_graph();

        // set memory for edge swap in csr
        csr->initEdgeSwap();

#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID node = 0; node < csr->num_nodes; node++) {
            if (is_k_bound_nodes[node]) {
                g.restore_node(node);
                // delete edges
                for (EDGE_ID neighor = begin_light(node); neighor < end_light(node); ) {
                    NODE_ID node_to = csr->adj[neighor];
                    if (!is_k_bound_nodes[node_to]) {
                        g.delete_edge_by_edge_swap(node, neighor, true);
                    }else{
                        neighor++;
                    }
                }

                for (EDGE_ID neighor = begin_heavy(node); neighor < end_heavy(node); ) {
                    NODE_ID node_to = csr->adj[neighor];
                    if (!is_k_bound_nodes[node_to]) {
                        g.delete_edge_by_edge_swap(node, neighor, false);
                    }else{
                        neighor++;
                    }
                }
            } else {
                // delete nodes
                g.delete_node_by_edge_swap(node);
            }
        }
    }

    //forward
    void init_edge_swap_with_2_pointers(bool* is_k_bound_nodes, const unsigned int num_threads) noexcept {
        const auto csr = g.get_csr_graph();

        // set memory for edge swap in csr
        csr->initEdgeSwap();

#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID node = 0; node < csr->num_nodes; node++) {
            if (is_k_bound_nodes[node]) {
                // delete light edges
                int i = csr->begin[node];
                int j = csr->begin[node] + csr->num_light_edges[node] - 1;
                int vertex_deleted = 0;
                while (i <= j) {
                    NODE_ID i_node = csr->adj[i];
                    if (!is_k_bound_nodes[i_node]) {
                        while (!is_k_bound_nodes[csr->adj[j]] && i <= j) {
                            vertex_deleted++;
                            j--;
                        }
                        if (i < j) {
                            // edge swap
                            w_type old_weight = csr->value[i];
                            csr->adj[i] = csr->adj[j];
                            csr->value[i] = csr->value[j];

                            csr->adj[j] = i_node;
                            csr->value[j] = old_weight;

                            vertex_deleted++;
                            i++;
                            j--;
                        }
                    } else {
                        i++;
                    }
                }
                csr->valid_light_num[node] = csr->num_light_edges[node] - vertex_deleted;

                // delete heavy edges
                i = csr->begin[node] + csr->num_light_edges[node];
                j = csr->begin[node + 1] - 1;
                vertex_deleted = 0;
                while (i <= j) {
                    NODE_ID i_node = csr->adj[i];
                    if (!is_k_bound_nodes[i_node]) {
                        while (!is_k_bound_nodes[csr->adj[j]] && i <= j) {
                            vertex_deleted++;
                            j--;
                        }
                        if (i < j) {
                            // edge swap
                            w_type old_weight = csr->value[i];
                            csr->adj[i] = csr->adj[j];
                            csr->value[i] = csr->value[j];

                            csr->adj[j] = i_node;
                            csr->value[j] = old_weight;

                            vertex_deleted++;
                            i++;
                            j--;
                        }
                    } else {
                        i++;
                    }
                }
                csr->valid_heavy_num[node] = csr->begin[node + 1] - csr->begin[node] - csr->num_light_edges[node] - vertex_deleted;

            } else {
                // delete nodes
                g.delete_node_by_edge_swap(node);
            }
        }
    }
    
};


#endif  // _KSP_GRAPH_H
