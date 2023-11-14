#ifndef _PROFILING_H
#define _PROFILING_H

#include "omp.h"
#include "Node.h"

class scoped_timer
{
    double& timer;
    double start;

public:
    scoped_timer(double& timer)
        : timer(timer),
          start(omp_get_wtime())
    {}

    ~scoped_timer()
    {
        timer += omp_get_wtime() - start;
    }
};

struct ProfilerStatistics
{
    unsigned int num_graph_copy;
    unsigned int num_sssp;
    unsigned int num_FSP;
    unsigned int num_branch_updates;
    unsigned int num_tree_copy;
    unsigned int num_graph_restore;
    unsigned int num_removed_edges;
    unsigned int num_simple_paths;
    double avg_path_size;
    EDGE_ID nodes_touched_in_ds; // TODO does not work in pd mode

    double time_graph_copy;
    double time_sssp_constructor;
    double time_sssp_compute;
    double time_sssp_total;
    double time_FSP;
    double time_branch_updates;
    double time_tree_copy;
    double time_graph_restore;
    double time_remove_edges;

    ProfilerStatistics()
        : num_graph_copy(0),
          num_sssp(0),
          num_FSP(0),
          num_branch_updates(0),
          num_tree_copy(0),
          num_graph_restore(0),
          num_removed_edges(0),
          num_simple_paths(0),
          avg_path_size(0),
          nodes_touched_in_ds(0),
          time_graph_copy(0),
          time_sssp_constructor(0),
          time_sssp_compute(0),
          time_sssp_total(0),
          time_FSP(0),
          time_branch_updates(0),
          time_tree_copy(0),
          time_graph_restore(0),
          time_remove_edges(0)
    {}

    std::string return_stats()
    {
        return "\"avg_path_size\":" + std::to_string(avg_path_size)
               + ",\"num_sssp\":" + std::to_string(num_sssp)
               + ",\"nodes_touched_in_ds\":" + std::to_string(nodes_touched_in_ds);
    }
/*
    void print_times()
    {
        std::cout << "graph copy:\nno of calls = " << num_graph_copy << "\ntotal time = " << time_graph_copy << std::endl;
        std::cout << "sssp:\nno calls = " << num_sssp << "\ntotal time = " << time_sssp_total << std::endl;
        std::cout << "sssp compute time = " << time_sssp_compute << std::endl;
        std::cout << "sssp constructor time = " << time_sssp_constructor << std::endl;
        std::cout << "sssp branch no calls = " << num_branch_updates << ", time = " << time_branch_updates << std::endl;
        std::cout << "simple paths: \nno calls = " << num_simple_paths << std::endl;
        std::cout << "FSP: \nno calls = " << num_FSP << "\ntotal time = " << time_FSP << std::endl;
        std::cout << "Graph restore: no calls = " << num_graph_restore << "\ntotal time = " << time_graph_restore << std::endl;
        std::cout << "remove edges: \nno calls = " << num_removed_edges << "\n total time = " << time_remove_edges << std::endl;
    }*/
};

#endif // _PROFILING_H
