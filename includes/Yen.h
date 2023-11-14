#ifndef _YEN_H
#define _YEN_H

#include <KSP.h>
#include <omp.h>
#include "Statistics.h"

template<template<typename, template<typename> class> class DeltaSteppingType, typename GraphType, unsigned int num_threads, bool pps = false>
class Yen : public KSP<GraphType, num_threads, pps>
{
    using ParentType = KSP<GraphType, num_threads, pps>;
    using DS = DeltaSteppingType<GraphType, KSPGraph>;
    using ParentType::parallel_ps;
    using ParentType::orig;
    using ParentType::paths;
    using ParentType::candidates;
    using ParentType::k;
#ifdef PROF
    using ParentType::profiler_statistics;
#endif

#ifdef DO_STATISTICS
    Statistics stats;
#endif

    std::vector<KSPGraph<GraphType>> graphs;

public:
    explicit Yen(const GraphType& g, const unsigned int k) noexcept
        : ParentType(g, k)
    {
        if(parallel_ps)
        {
#ifdef DO_STATISTICS
            std::cerr << "Statistics are not supported during parallel computation." << std::endl;
#endif

            graphs.reserve(num_threads);
            for(unsigned int i = 1; i < num_threads; i++)
                graphs.push_back(KSPGraph<GraphType>(g));

            graphs.push_back(orig);
        }
    }

#ifdef DO_STATISTICS
    Statistics get_stats() const
    {
        return stats;
    }
#endif

    std::vector<Path> compute(const NODE_ID source, const NODE_ID destination)
    {
        assert(k > 1);

        paths.reserve(k);

#ifdef TECHNIQUE_BRIDGE_WCC
        //check whether source and destination are in the same wcc
        CSRGraph* csr=orig.get_original_graph().get_csr_graph();
        if(csr->wcc && csr->wcc[source]!=csr->wcc[destination])
        {
            std::cout<<"wcc is different between deviation_node:"<<source<<" and destination:"<<destination<<std::endl;
            return paths;
        }
#endif

        if(parallel_ps)
            return find_path_par(source, destination);

        return find_path_seq(source, destination);
    }

private:
    std::vector<Path> find_path_seq(const NODE_ID source, const NODE_ID destination)
    {
#ifdef PROF
        profiler_statistics.num_sssp++;
        double start = omp_get_wtime();
#endif

        DS d(orig, num_threads);

#ifdef PROF
        const double c_start = omp_get_wtime();
        profiler_statistics.time_sssp_constructor += c_start - start;
        d.template compute<SSSP_EARLY_STOPPING>(profiler_statistics, source, destination);
        profiler_statistics.time_sssp_compute += omp_get_wtime() - c_start;
#else
        d.template compute<SSSP_EARLY_STOPPING>(source, destination);
#endif

        candidates.insert({d.get_distance(destination), 0, 0, d.get_shortest_path_to(destination)});

#ifdef PROF
        profiler_statistics.time_sssp_total += omp_get_wtime() - start;
#endif

#ifdef DO_STATISTICS
        stats = Statistics(source, destination, k);
#endif
        for(unsigned int kp = 0; kp < k && !candidates.empty(); kp++)
        {
            paths.push_back(candidates.get_best());

#ifdef DO_STATISTICS
            if(kp == 0 && paths[0].p.size() == 2)
                throw NoPathException();    // not thread safe

            stats.current_k = kp;
            stats.k_shortest_paths[stats.current_k] = paths[stats.current_k].p;
            stats.hops[stats.current_k] = static_cast<unsigned int>(paths[stats.current_k].p.size() - 1);
            stats.lengths[stats.current_k] = paths[stats.current_k].length;
            if(stats.current_k > 0)
            {
                for(const auto u : paths[stats.current_k].p)
                {
                    if(std::find(paths[paths[stats.current_k].parent_index].p.begin(),
                                 paths[paths[stats.current_k].parent_index].p.end(), u) !=
                       paths[paths[stats.current_k].parent_index].p.end())
                    {
                        stats.shared_nodes[stats.current_k]++;
                    }
                }
                stats.parent_id[stats.current_k] = paths[stats.current_k].parent_index;
                stats.deviation_node_id[stats.current_k] = paths[stats.current_k].alpha;
            }
#endif

            if(kp == k - 1)
                break;

            for(unsigned int i = 0; i < paths[kp].alpha; i++)
            {
#ifdef PROF
                profiler_statistics.num_removed_edges += orig.get_num_neighbors(paths[kp].p[i]) * 2;
                scoped_timer t(profiler_statistics.time_remove_edges);
#endif

                orig.remove_node(paths[kp].p[i]);
            }

            for(unsigned int i = 0; i < kp; i++)
            {
                if(i == paths[kp].parent_index)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    scoped_timer t(profiler_statistics.time_remove_edges);
#endif

                    remove_edge(orig, paths[i].p[paths[kp].alpha], paths[i].p[paths[kp].alpha + 1]);
                }

                if(paths[i].parent_index == paths[kp].parent_index && paths[i].alpha == paths[kp].alpha)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    scoped_timer t(profiler_statistics.time_remove_edges);
#endif

                    remove_edge(orig, paths[i].p[paths[i].alpha], paths[i].p[paths[i].alpha + 1]);
                }
            }

            // now compute intermediate paths from alpha to end
            for(unsigned int i = paths[kp].alpha; i < paths[kp].p.size() - 1; i++)
            {
                const NODE_ID node = paths[kp].p[i];

#ifdef TECHNIQUE_BOUND
                std::vector<NODE_ID> prefix_path = paths[kp].getSubPath(i+1);
                w_type prefix_path_length=orig.get_original_graph().get_path_length(prefix_path);
                w_type k_bound=candidates.getKBound();
                w_type sub_k_bound=k_bound-prefix_path_length;

                d.setKBound(sub_k_bound);
#endif

#ifdef PROF
                profiler_statistics.num_removed_edges += 2;
                start = omp_get_wtime();
#endif

#ifdef TECHNIQUE_BRIDGE_WCC
                //check whether the edge of deviation_node and next one is a bridge
                CSRGraph* csr=orig.get_original_graph().get_csr_graph();        
                if(csr->bridge && csr->bridge->count(std::make_pair(node,paths[kp].p[i + 1])))
                {
                    std::cout<<"find_path_seq the deleting edge ("<<node<<","<<paths[kp].p[i + 1]<<") is a bridge"<<std::endl;
                    remove_edge(orig, node, paths[kp].p[i + 1]);
                    orig.remove_node(node);
                    continue;
                }
#endif

                remove_edge(orig, node, paths[kp].p[i + 1]);

#ifdef PROF
                profiler_statistics.time_remove_edges += omp_get_wtime() - start;
                const double t1 = omp_get_wtime();
                d.template compute<SSSP_EARLY_STOPPING>(profiler_statistics, node, destination, candidates.get_length_threshold() - orig.get_original_graph().get_path_length(paths[kp].p, i));
                profiler_statistics.time_sssp_compute += omp_get_wtime() - t1;
                profiler_statistics.time_sssp_total += omp_get_wtime() - t1;
                profiler_statistics.num_sssp++;
#else
                d.template compute<SSSP_EARLY_STOPPING>(node, destination, candidates.get_length_threshold() - orig.get_original_graph().get_path_length(paths[kp].p, i));
#endif

#ifdef TECHNIQUE_BOUND
                if (sub_k_bound <= d.get_distance(destination))
                {
                    // it could get a path which is not the shortest if early terminate with TECHNIQUE_BOUND
                    orig.remove_node(node);
                    continue;
                }
#endif
                // concatenate paths
                if(d.get_distance(destination) < candidates.get_length_threshold() - orig.get_original_graph().get_path_length(paths[kp].p, i))
                {
                    // get the end of the path
                    std::vector<NODE_ID> new_path_end = d.get_shortest_path_to(destination);
                    // prepare a temporary path (we be removed by the optimizer)
                    std::vector<NODE_ID> new_path;
                    new_path.reserve(i + new_path_end.size());
                    // put the paths together
                    new_path.insert(new_path.end(), paths[kp].p.begin(), paths[kp].p.begin() + i);
                    new_path.insert(new_path.end(), new_path_end.begin(), new_path_end.end());

                    candidates.insert({orig.get_original_graph().get_path_length(new_path),
                                            i == paths[kp].alpha ? paths[kp].parent_index : kp, i, new_path});
                }

                // remove this node from graph
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += orig.get_num_neighbors(node) * 2;
                    scoped_timer t(profiler_statistics.time_remove_edges);
#endif
                    orig.remove_node(node);
                }
            }

            {
#ifdef PROF
                profiler_statistics.num_graph_restore++;
                scoped_timer t(profiler_statistics.time_graph_restore);
#endif
                orig.restore_graph();
            }
        }

#ifdef PROF
        profiler_statistics.avg_path_size = this->get_avg_path_size();
#endif

        return paths;
    }

    std::vector<Path> find_path_par(const NODE_ID source, const NODE_ID destination)
    {
#if !KSP_PD_L1
        omp_set_nested(1);
#endif
#ifdef PROF
        profiler_statistics.num_sssp++;
        const double start = omp_get_wtime();
        const double construct_time = start;
#endif

        std::vector<std::unique_ptr<DS>> d;
#if KSP_PD_L1
        d.push_back(std::make_unique<DS>(graphs[0], num_threads));
        for(unsigned int i = 1; i < num_threads; i++)
            d.push_back(std::make_unique<DS>(graphs[i], 1));
#else
        for(unsigned int i = 0; i < num_threads; i++)
            d.push_back(std::make_unique<DS>(graphs[i], num_threads));
#endif

#ifdef PROF
        const double c_start = omp_get_wtime();
        profiler_statistics.time_sssp_constructor += c_start - construct_time;
#endif

#ifdef PROF
        d[0]->template compute<SSSP_EARLY_STOPPING>(profiler_statistics, source, destination);
#else
        d[0]->template compute<SSSP_EARLY_STOPPING>(source, destination);
#endif

#ifdef PROF
        profiler_statistics.time_sssp_compute += omp_get_wtime() - c_start;
#endif

        candidates.insert({d[0]->get_distance(destination), 0, 0, d[0]->get_shortest_path_to(destination)});

#ifdef PROF
        profiler_statistics.time_sssp_total += omp_get_wtime() - start;
#endif
#if KSP_PD_L1
        d[0]->reset_num_partition_threads(1);
#else
        std::vector<unsigned int> threads(num_threads);
#endif

        for(unsigned int kp = 0; kp < k && !candidates.empty(); kp++)
        {
#if !KSP_PD_L1
            std::fill(threads.begin(), threads.end(), 1);
#endif

            paths.push_back(candidates.get_best());

            if(kp == k - 1)
                break;

            const auto int_paths = static_cast<unsigned int>(paths[kp].p.size() - paths[kp].alpha - 1);

            if(int_paths == 0)
                continue;

#if !KSP_PD_L1
            if(int_paths < num_threads)
            {
                const unsigned int num_t = num_threads / int_paths;
                for(size_t i = 0; i < int_paths; i++)
                    threads[i] = num_t;
                for(size_t i = num_threads - (num_t * int_paths); i > 0; i--)
                    threads[i]++;
            }
#endif

            unsigned int paral_paths = int_paths > num_threads ? num_threads : int_paths;
            NODE_ID end_paral_paths = paral_paths + paths[kp].alpha;
            calculate_candidates_parallel(kp, d
#if !KSP_PD_L1
                , threads
#endif
                , destination, paths[kp].alpha, end_paral_paths);

            if(paral_paths == int_paths)        // by definition paral_paths <= int_paths holds
                continue;

            unsigned int num_rest_paths = int_paths - paral_paths;
            while(num_rest_paths > 0)
            {
#if !KSP_PD_L1
                if(num_rest_paths < num_threads)
                {
                    const unsigned int num_t = num_threads / num_rest_paths;

                    for(size_t i = 0; i < num_rest_paths; i++)
                        threads[i] = num_t;

                    for(size_t i = num_threads - (num_t * num_rest_paths); i > 0; i--)
                        threads[i]++;
                }
#endif
                paral_paths = num_rest_paths > num_threads ? num_threads : num_rest_paths;
                calculate_candidates_parallel(kp, d
#if !KSP_PD_L1
                    , threads
#endif
                    , destination, end_paral_paths, end_paral_paths + paral_paths);

                end_paral_paths += paral_paths;
                num_rest_paths -= paral_paths;
            }
        }

#ifdef PROF
        profiler_statistics.avg_path_size = this->get_avg_path_size();
#endif

        return paths;
    }

    void calculate_candidates_parallel(const unsigned int kp, const std::vector<std::unique_ptr<DS>>& d
#if !KSP_PD_L1
        , const std::vector<unsigned int>& threads
#endif
        , const NODE_ID destination, const unsigned int deviation_start, const size_t deviation_end) noexcept
    {
        const size_t num_deviations = deviation_end - deviation_start;

        std::vector<Path> candidates_buffer(num_deviations);

#pragma omp parallel for num_threads(num_deviations)   // find_path_par()
        for(unsigned int inter_p = deviation_start; inter_p < deviation_end; inter_p++)
        {
            const auto t_id = static_cast<unsigned int>(omp_get_thread_num());
            const NODE_ID node = paths[kp].p[inter_p];

#ifdef TECHNIQUE_BRIDGE_WCC
            //check whether the edge of deviation_node and next one is a bridge
            CSRGraph* csr=orig.get_original_graph().get_csr_graph();        
            if(csr->bridge && csr->bridge->count(std::make_pair(node,paths[kp].p[inter_p + 1])))
            {
                std::cout<<"find_path_par the deleting edge ("<<node<<","<<paths[kp].p[inter_p + 1]<<") is a bridge"<<std::endl;
                continue;
            }
#endif

#ifdef TECHNIQUE_BOUND
            std::vector<NODE_ID> prefix_path = paths[kp].getSubPath(inter_p+1);
            w_type prefix_path_length=orig.get_original_graph().get_path_length(prefix_path);
            w_type k_bound=candidates.getKBound();
            w_type sub_k_bound=k_bound-prefix_path_length;

            d[t_id]->setKBound(sub_k_bound);
#endif


#if !KSP_PD_L1
            d[t_id]->reset_num_partition_threads(threads[t_id]);
#endif
            remove_nodes(kp, t_id, inter_p);

            // remove previous node from graph
            for(unsigned int i = paths[kp].alpha; i < inter_p; i++)
            {
                graphs[t_id].remove_node(paths[kp].p[i]);
            }

            // now compute intermediate paths
            
            remove_edge(graphs[t_id], node, paths[kp].p[inter_p + 1]);

#ifdef PROF
            d[t_id]->template compute<SSSP_EARLY_STOPPING>(profiler_statistics, node, destination, candidates.get_length_threshold() - orig.get_original_graph().get_path_length(paths[kp].p, inter_p));
#else
            d[t_id]->template compute<SSSP_EARLY_STOPPING>(node, destination, candidates.get_length_threshold() - orig.get_original_graph().get_path_length(paths[kp].p, inter_p));
#endif

#ifdef TECHNIQUE_BOUND
            if (sub_k_bound <= d[t_id]->get_distance(destination))
            {
                // it could get a path which is not the shortest if early terminate with TECHNIQUE_BOUND
                graphs[t_id].restore_graph();
                continue;
            }
#endif

            // concatenate paths
            if(d[t_id]->get_distance(destination) < candidates.get_length_threshold())
            {
                const std::vector<NODE_ID> new_path_end = d[t_id]->get_shortest_path_to(destination);
                // prepare temporary new path
                std::vector<NODE_ID> new_path;
                new_path.reserve(inter_p + new_path_end.size());
                // concatenate the paths
                new_path.insert(new_path.end(), paths[kp].p.begin(), paths[kp].p.begin() + inter_p);
                new_path.insert(new_path.end(), new_path_end.begin(), new_path_end.end());

                candidates_buffer[inter_p - deviation_start] = {orig.get_original_graph().get_path_length(new_path),
                                                                  inter_p == paths[kp].alpha ? paths[kp].parent_index : kp, inter_p, new_path};
            }

            graphs[t_id].restore_graph();
        }

        candidates.insert(candidates_buffer);
    }

    void remove_edge(KSPGraph<GraphType>& this_g, const NODE_ID src, const NODE_ID dest) const noexcept
    {
        EDGE_ID edge_id = this_g.get_edge_id(src, dest);

        if(edge_id != this_g.end(src))
            this_g.remove_edge(edge_id);
    }

    void remove_nodes(const unsigned int kp, const unsigned int t_id, const unsigned int deviating_at) noexcept
    {
        if(paths[kp].alpha > 0)
        {
            for(unsigned int i = 0; i < paths[kp].alpha; i++)
            {
                graphs[t_id].remove_node(paths[kp].p[i]);
            }
        }

        if(paths[kp].parent_index < kp)
        {
            remove_edge(graphs[t_id], paths[paths[kp].parent_index].p[paths[kp].alpha],
                        paths[paths[kp].parent_index].p[paths[kp].alpha + 1]);
        }

        for(unsigned int i = 0; i < kp; i++)
        {
            if(paths[i].parent_index == paths[kp].parent_index
               && paths[i].alpha == paths[kp].alpha
               && paths[kp].alpha >= deviating_at)
            {
                remove_edge(graphs[t_id], paths[i].p[paths[i].alpha], paths[i].p[paths[i].alpha + 1]);
            }
        }
    }
};

#endif  // _YEN_H
