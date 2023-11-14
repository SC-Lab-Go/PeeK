#ifndef _BASIC_GRAPH_H
#define _BASIC_GRAPH_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <Node.h>
#include <omp.h>
#include <cassert>
#include <memory>

/**
 * Stores nodes and edges of a graph
 */
template<bool directed, bool weighted>
class BasicGraph
{
    friend class GraphRW;

private:
    const NODE_ID num_nodes;        // holds the number of actual nodes, without the sentinel
    const EDGE_ID num_edges;        // holds the number of edges (edges.size())

    const w_type heaviest_weight;
    const w_type delta;

    std::unique_ptr<const BasicGraph<directed, weighted>> reverse_graph;
    
    CSRGraph* csr;

public:
    static constexpr bool is_directed = directed;
    static constexpr bool is_weighted = weighted;

    BasicGraph(CSRGraph* csr, w_type heaviest_weight, w_type delta,
               std::unique_ptr<const BasicGraph<directed, weighted>>&& reverse_graph = std::unique_ptr<const BasicGraph<directed, weighted>>(nullptr))
        : csr(csr),
          num_nodes(csr->num_nodes),
          num_edges(csr->num_edges),
          heaviest_weight(heaviest_weight),
          delta(delta),
          reverse_graph(std::move(reverse_graph))
    {}
    
    CSRGraph* get_csr_graph() const noexcept
    {
        return csr;
    }

    /**
     * Returns the reverse graph. This is for undirected graphs the original graph.
     * For the directed case the reverse graph needs to be computed explicitly and provided to the constructor.
     * @return
     */
    const BasicGraph<directed, weighted>& get_reverse_graph() const noexcept
    {
        if(!directed)
            return *this;

        assert(reverse_graph != nullptr);

        return *reverse_graph;
    }

    /**
     * Returns the number of actual nodes without the sentinel.
     * @return The number of actual nodes without the sentinel.
     */
    NODE_ID get_num_nodes() const noexcept
    {
        return csr->num_nodes;
    }

    /**
     * Returns the number of edges.
     * @return The number of edges.
     */
    EDGE_ID get_num_edges() const noexcept
    {
        return csr->num_edges;
    }

    /**
     * Returns the number of outgoing edges of node.
     * @param node_id
     * @return The number of outgoing edges of node.
     */
    NODE_ID get_num_neighbors(const NODE_ID node_id) const noexcept
    {
        return static_cast<NODE_ID>(csr->begin[node_id + 1] - csr->begin[node_id]);   // works thanks to the sentinel
    }

    /**
     * Returns the weight of the edge (u,v).
     * @param u
     * @param v
     * @return The weight of the edge (u,v). If the edge does not exists, infinity is returned.
     */
    w_type get_edge_weight(const NODE_ID u, const NODE_ID v) const noexcept
    {
        for ( EDGE_ID neighor = csr->begin[u]; neighor < csr->begin[u + 1]; neighor++ )
        {
            w_type  weight  = csr->value[neighor];
            NODE_ID node_to = csr->adj[neighor];
            if(node_to == v)
            {
                return weight;
            }
        }

        return INFINITE_DISTANCE;
    }


    /**
     * Returns the begin of the edge list of the node.
     * @param node
     * @return begin, such that all edges starting in node are in [begin, end)
     */
    EDGE_ID get_edge_list_begin(const NODE_ID node) const noexcept
    {
        return csr->begin[node];
    }

    /**
     * Returns the begin of the heavy edge list of the node. (alias for get_light_edge_list_end())
     * @param node
     * @return begin, such that all heavy edges starting in node are in [begin, end)
     */
    EDGE_ID get_heavy_edge_list_begin(const NODE_ID node) const noexcept
    {
        return get_light_edge_list_end(node);
    }

    /**
     * Returns the end of the edge list of the node.
     * @param node
     * @return end, such that all edges starting in node are in [begin, end)
     */
    EDGE_ID get_edge_list_end(const NODE_ID node) const noexcept
    {
        return csr->begin[node + 1];      // works thanks to the sentinel
    }

    /**
     * Returns the end of the light edge list of the node.
     * @param node
     * @return end, such that all light edges starting in node are in [begin, end)
     */
    EDGE_ID get_light_edge_list_end(const NODE_ID node) const noexcept
    {
        return csr->begin[node] + csr->num_light_edges[node];
    }

    EDGE_ID get_light_edge_list_end_edge_swap(const NODE_ID node) const noexcept { return csr->begin[node] + csr->valid_light_num[node]; }

    EDGE_ID get_heavy_edge_list_end_edge_swap(const NODE_ID node) const noexcept { return csr->begin[node] + csr->num_light_edges[node] + csr->valid_heavy_num[node]; }

    void delete_node_by_edge_swap(const NODE_ID node) const noexcept {
        csr->valid_light_num[node] = 0;
        csr->valid_heavy_num[node] = 0;
    }

    void delete_edge_by_edge_swap(const NODE_ID node, const EDGE_ID e_id, bool is_light) const noexcept {
        if (e_id >= get_edge_list_begin(node) && e_id < get_edge_list_end(node)) {
            NODE_ID old_node_to = csr->adj[e_id];
            w_type old_weight = csr->value[e_id];
            EDGE_ID last_valid_index;

            if (is_light) {
                assert(csr->valid_light_num[node] > 0);
                last_valid_index = get_light_edge_list_end_edge_swap(node) - 1;
                csr->valid_light_num[node]--;
            } else {
                assert(csr->valid_heavy_num[node] > 0);
                last_valid_index = get_heavy_edge_list_end_edge_swap(node) - 1;
                csr->valid_heavy_num[node]--;
            }

            csr->adj[e_id] = csr->adj[last_valid_index];
            csr->value[e_id] = csr->value[last_valid_index];

            csr->adj[last_valid_index] = old_node_to;
            csr->value[last_valid_index] = old_weight;
        }
    }

    void restore_node(const NODE_ID node) const noexcept {
        EDGE_ID total_neighors = get_edge_list_end(node) - get_edge_list_begin(node);
        NODE_ID light_nodes = csr->num_light_edges[node];
        NODE_ID heavy_nodes = total_neighors - light_nodes;
        csr->valid_light_num[node] = light_nodes;
        csr->valid_heavy_num[node] = heavy_nodes;
    }

    /**
     * Returns the delta value the edges were classified with.
     * @return The delta value the edges were classified with.
     */
    w_type get_delta() const noexcept
    {
        return delta;
    }

    /**
     * Returns the bucket size used in DeltaStepping.
     * @return
     */
    unsigned int compute_bucket_size() const noexcept
    {
        return static_cast<unsigned int>(ceil(heaviest_weight / delta)) + 1;
    }

    /**
     * Retruns the average degree calculated as num_edges / num_nodes.
     * @return The average degree.
     */
    double get_avg_degree() const noexcept
    {
        return static_cast<double>(csr->num_edges) / static_cast<double>(csr->num_nodes);
    }

    w_type get_path_length(const std::vector<NODE_ID>& node_list) const noexcept
    {
        return get_path_length(node_list, node_list.size() - 1);
    }

    w_type get_path_length(const std::vector<NODE_ID>& node_list, const size_t end) const noexcept
    {
        if(node_list.empty())
            return 0;

        assert(end < node_list.size());

        w_type length = 0;

        for(size_t i = 0; i < end; i++)
        {
            length += get_edge_weight(node_list[i], node_list[i + 1]);
        }

        return length;
    }

};

#endif  // _BASIC_GRAPH_H
