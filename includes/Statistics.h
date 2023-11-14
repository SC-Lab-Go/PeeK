//
// Created by Alex Schickedanz <alex@ae.cs.uni-frankfurt.de> on 14.08.18.
//

#ifndef K_SHORTEST_PATHS_CODE_STATISTICS_H
#define K_SHORTEST_PATHS_CODE_STATISTICS_H

#include <vector>
#include "Node.h"

struct Statistics
{
    // helper to find the right place for statistics
    unsigned int current_k = 0;

    // Source and destination to be able to reproduce the results
    NODE_ID source, destination;
    // The k shortest paths them self
    std::vector<path_type> k_shortest_paths;
    // Number of yellow nodes for each of the k shortest paths and each deviation node
    std::vector<std::vector<size_t>> num_yellow_nodes;
    // Number of edges in the yellow graph for each of the k shortest paths and each deviation node
    std::vector<std::vector<size_t>> num_yellow_edges;
    // Number of red nodes for each of the k shortest paths and each deviation node
    std::vector<std::vector<size_t>> num_red_nodes;
    // Number of in edges to all red notes (most of them are not part of the yellow graph)
    std::vector<std::vector<size_t>> num_red_edges;
    // The length of each of each of the shortest paths
    std::vector<w_type> lengths;
    // The number of hops of each of the shortest paths
    std::vector<unsigned int> hops;
    // Number of nodes a path and its parend share
    std::vector<unsigned int> shared_nodes;
    // Number of nodes a path and its parend share
    std::vector<unsigned int> deviation_node_id;
    // Number of nodes a path and its parend share
    std::vector<unsigned int> parent_id;
    // depth of shortest path tree for the ith node of the kth shortest path
    std::vector<std::vector<unsigned int>> shortest_path_tree_depth;
    // the number of nodes hanging on the ith node of the kth shortest path in the shortest path tree
    std::vector<std::vector<unsigned int>> shortest_path_tree_hanging_nodes;

    explicit Statistics(const NODE_ID source = 0, const NODE_ID destination = 0, const unsigned int k = 0)
        : source(source),
          destination(destination),
          k_shortest_paths(k),
          num_yellow_nodes(k),
          num_yellow_edges(k),
          num_red_nodes(k),
          num_red_edges(k),
          lengths(k, 0),
          hops(k, 0),
          shared_nodes(k, 0),
          deviation_node_id(k, 0),
          parent_id(k, 0),
          shortest_path_tree_depth(k),
          shortest_path_tree_hanging_nodes(k)
    {}

    void print_as_JSON(const std::string& algorithm = "", const std::string& graph = "")
    {
        printf(R"(,{"algorithm":"%s","graph":"%s","source":%u,"destination":%u,"data":[)", algorithm.c_str(), graph.c_str(), source, destination);
        for(unsigned int i = 0; i <= current_k; ++i)
        {
            if(i > 0)
                printf(",");

            printf(R"({"length":%.5f,"hops":%u,"shared_nodes":%u,"deviation_node_id":%u,"parent_path_id":%u)",
                   lengths[i], hops[i], shared_nodes[i], deviation_node_id[i], parent_id[i]);
            printf(",\"num_yellow_nodes\":");
            vector_to_JSON(num_yellow_nodes[i]);
            printf(",\"num_yellow_edges\":");
            vector_to_JSON(num_yellow_edges[i]);
            printf(",\"num_red_nodes\":");
            vector_to_JSON(num_red_nodes[i]);
            printf(",\"num_red_edges\":");
            vector_to_JSON(num_red_edges[i]);
            printf(",\"depth_in_reverse_sssp_tree\":");
            vector_to_JSON(shortest_path_tree_depth[i]);
            printf(",\"num_hanging_nodes\":");
            vector_to_JSON(shortest_path_tree_hanging_nodes[i]);
            printf(",\"path\":");
            vector_to_JSON(k_shortest_paths[i]);
            printf("}");
        }

        printf("]}");
    }

private:
    template <typename T>
    void vector_to_JSON(const std::vector<T>& vec)
    {
        if(vec.empty())
        {
            printf("[]");
            return;
        }

        printf("[%lu", static_cast<unsigned long>(vec[0]));
        for(unsigned int j = 1; j < vec.size(); j++)
            printf(",%lu", static_cast<unsigned long>(vec[j]));
        printf("]");
    }
};


#endif //K_SHORTEST_PATHS_CODE_STATISTICS_H
