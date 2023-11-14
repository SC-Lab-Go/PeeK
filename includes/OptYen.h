#ifndef _OPT_YEN_H
#define _OPT_YEN_H

#include <KSP.h>
#include <algorithm>
#include <omp.h>

template<template<typename, template<typename> class> class DeltaSteppingType, typename GraphType, unsigned int num_threads, bool pps = false>
class OptYen : public KSP<GraphType, num_threads, pps>
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

    std::vector<KSPGraph<GraphType>> graphs;
    std::unique_ptr<KSPGraph<GraphType>> reverse_graph;
    SsspTree sssp_t;

public:
    explicit OptYen(const GraphType& g, const unsigned int k) noexcept
        : ParentType(g, k),
          reverse_graph(GraphType::is_directed ? std::make_unique<KSPGraph<GraphType>>(g.get_reverse_graph()) : nullptr),
          sssp_t(g.get_num_nodes())
    {
        if(parallel_ps)
        {
            graphs.reserve(num_threads);
            for(unsigned int i = 1; i < num_threads; i++)
                graphs.push_back(KSPGraph<GraphType>(g));

            graphs.push_back(orig);
        }
    }

    explicit OptYen(const GraphType& g, const unsigned int k, const SsspTree& sssp_t, std::vector<Path> path) noexcept
        : ParentType(g, k),
          sssp_t(sssp_t)
    {
        candidates.insert(path);
        if(parallel_ps)
        {
            graphs.reserve(num_threads);
            for(unsigned int i = 1; i < num_threads; i++)
                graphs.push_back(KSPGraph<GraphType>(g));

            graphs.push_back(orig);
        }
    }

    explicit OptYen(const GraphType& g, const unsigned int k, const SsspTree& sssp_t, std::vector<Path> path, bool* is_k_bound_nodes) noexcept : ParentType(g, k), sssp_t(sssp_t) {
        candidates.insert(path);

#if defined(ADAPTIVE_GRAPH_COMPACT_VERTEX) || defined(ADAPTIVE_GRAPH_COMPACT_EDGE)
        const NODE_ID num_nodes = orig.get_num_nodes();
        const EDGE_ID num_edges = orig.get_num_edges();
        auto csr = orig.get_original_graph().get_csr_graph();
#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID node = 0; node < num_nodes; node++) {
            if (is_k_bound_nodes[node]) {
#ifdef ADAPTIVE_GRAPH_COMPACT_EDGE
                const EDGE_ID edge_list_end = orig.end(node);
                for (EDGE_ID neighor = orig.begin(node); neighor < edge_list_end; neighor++) {
                    NODE_ID node_to = csr->adj[neighor];
                    if (!is_k_bound_nodes[node_to]) {
                        orig.edge_is_removed[neighor] = true;
                    }
                }
#endif
            } else {
#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
                orig.node_is_removed[node] = true;
#endif
#ifdef ADAPTIVE_GRAPH_COMPACT_EDGE
                const EDGE_ID edge_list_end = orig.end(node);
                for (EDGE_ID neighor = orig.begin(node); neighor < edge_list_end; neighor++) {
                    orig.edge_is_removed[neighor] = true;
                }
#endif
            }
        }
#endif

        if (parallel_ps) {
            graphs.reserve(num_threads);
            for (unsigned int i = 1; i < num_threads; i++) {
                KSPGraph<GraphType> new_ksp_graph = KSPGraph<GraphType>(g);

#ifdef ADAPTIVE_GRAPH_COMPACT_VERTEX
#pragma omp parallel for num_threads(num_threads)
                for (NODE_ID node = 0; node < num_nodes; node++) {
                    if (orig.node_is_removed[node]) {
                        new_ksp_graph.node_is_removed[node] = true;
                    }
                }
#endif
#ifdef ADAPTIVE_GRAPH_COMPACT_EDGE
#pragma omp parallel for num_threads(num_threads)
                for (EDGE_ID edge = 0; edge < num_edges; edge++) {
                    if (orig.edge_is_removed[edge]) {
                        new_ksp_graph.edge_is_removed[edge] = true;
                    }
                }
#endif

                graphs.push_back(new_ksp_graph);
            }

            graphs.push_back(orig);
        }
    }

    std::vector<Path> compute(const NODE_ID source, const NODE_ID destination, bool run_reverse_graph=true)
    {
        assert(k > 1);

        paths.reserve(k);

        if(parallel_ps)
            return compute_pd(source, destination, run_reverse_graph);

        return compute_psp(source, destination, run_reverse_graph);
    }

private:
    std::vector<Path> compute_psp(const NODE_ID source, const NODE_ID destination, bool run_reverse_graph=true)
    {
#ifdef PROF
        profiler_statistics.num_sssp++;
        double start = omp_get_wtime();
        double c_start;
        const double construct_time = omp_get_wtime();
#endif
        std::unique_ptr<DS> d;

        if(run_reverse_graph){

            std::unique_ptr<DS> d_rev = std::make_unique<DS>(GraphType::is_directed ? *reverse_graph : orig, num_threads);

#ifdef PROF
            c_start = omp_get_wtime();
            profiler_statistics.time_sssp_constructor += c_start - construct_time;
#endif

            d_rev->template compute<false>(
#ifdef PROF
            profiler_statistics,
#endif
            destination, source);

#ifdef PROF
            profiler_statistics.time_sssp_compute += omp_get_wtime() - c_start;
#endif

            sssp_t = d_rev->get_sssp_tree();

#ifdef PROF
            profiler_statistics.time_sssp_total += omp_get_wtime() - start;
#endif
            candidates.insert({d_rev->get_distance(source), 0, 0, d_rev->get_shortest_path_to(source, true)});
        
            if(GraphType::is_directed)
                d = std::make_unique<DS>(orig, num_threads);
            else
                d = std::move(d_rev);
        }else{
            d = std::make_unique<DS>(orig, num_threads);
        }

        for(unsigned int kp = 0; kp < k && !candidates.empty(); kp++)
        {
            paths.push_back(candidates.get_best());

            if(kp == k - 1)
                break;

            if(paths[kp].alpha > 0)
            {
                for(unsigned int i = 0; i < paths[kp].alpha; i++)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += orig.get_num_neighbors(paths[kp].p[i]) * (1 + static_cast<unsigned int>(!GraphType::is_directed));
                    const double rm_node_start = omp_get_wtime();
#endif

                    orig.remove_node(paths[kp].p[i]);

#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - rm_node_start;
#endif
                }
            }

            std::vector<NODE_ID> avoid;

            for(unsigned int i = 0; i < kp; i++)
            {
                if(i == paths[kp].parent_index)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    start = omp_get_wtime();
#endif

                    remove_edge(paths[i].p[paths[kp].alpha], paths[i].p[paths[kp].alpha + 1]);

#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - start;
#endif

                    avoid.push_back(paths[i].p[paths[kp].alpha + 1]);
                }

                if(paths[i].parent_index == paths[kp].parent_index && paths[i].alpha == paths[kp].alpha)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    start = omp_get_wtime();
#endif
                    remove_edge(paths[i].p[paths[i].alpha], paths[i].p[paths[i].alpha + 1]);
#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - start;
#endif
                    avoid.push_back(paths[i].p[paths[i].alpha + 1]);
                }
            }

            for(unsigned int i = paths[kp].alpha; i < paths[kp].p.size() - 1; i++)
            {
                const NODE_ID node = paths[kp].p[i];


                avoid.push_back(paths[kp].p[i + 1]);
                Path p = compute_simple_path(orig, destination, paths[kp].p, i + 1, avoid);

                if(p.length < INFINITE_DISTANCE)
                {
#ifdef PROF
                    profiler_statistics.num_simple_paths++;
#endif
                }
                else
                {
                    p.p.clear();
                    
#ifdef TECHNIQUE_BOUND
                    std::vector<NODE_ID> prefix_path = paths[kp].getSubPath(i+1);
                    w_type prefix_path_length=orig.get_original_graph().get_path_length(prefix_path);
                    w_type k_bound=candidates.getKBound();
                    w_type sub_k_bound=k_bound-prefix_path_length;

                    d->setKBound(sub_k_bound);
#endif

#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    start = omp_get_wtime();
#endif

                    remove_edge(node, paths[kp].p[i + 1]);

#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - start;
                    const double t1 = omp_get_wtime();
#endif

                    d->template compute<SSSP_EARLY_STOPPING>(
#ifdef PROF
                        profiler_statistics,
#endif
                        node, destination, candidates.get_length_threshold() -  - orig.get_original_graph().get_path_length(paths[kp].p, i)
                    );

#ifdef PROF
                    profiler_statistics.time_sssp_compute += omp_get_wtime() - t1;
                    profiler_statistics.time_sssp_total += omp_get_wtime() - t1;
                    profiler_statistics.num_sssp++;
#endif

#ifdef TECHNIQUE_BOUND
                    if (sub_k_bound < d->get_distance(destination))
                    {
                        // it could get a path which is not the shortest if early terminate with TECHNIQUE_BOUND
                        orig.remove_node(node);
                        avoid.clear();
                        continue;
                    }
#endif

                    // concatenate paths
                    if(d->get_distance(destination) < candidates.get_length_threshold())
                    {
                        p.p.insert(p.p.end(), paths[kp].p.begin(), paths[kp].p.begin() + i);

                        const std::vector<NODE_ID> new_path_end = d->get_shortest_path_to(destination);
                        p.p.insert(p.p.end(), new_path_end.begin(), new_path_end.end());
                        p.length = orig.get_original_graph().get_path_length(p.p);
                        p.alpha = i;
                    }
                }

                if(p.length < candidates.get_length_threshold())
                {
                    p.parent_index = p.alpha == paths[kp].alpha ? paths[kp].parent_index : kp;
                    candidates.insert(std::move(p));
                }

#ifdef PROF
                profiler_statistics.num_removed_edges += orig.get_num_neighbors(node) * (1 + static_cast<unsigned int>(!GraphType::is_directed));
                start = omp_get_wtime();
#endif

                orig.remove_node(node);

#ifdef PROF
                profiler_statistics.time_remove_edges += omp_get_wtime() - start;
#endif

                avoid.clear();
            }
            // the end of intermediate paths

#ifdef PROF
            const double t = omp_get_wtime();
#endif

            orig.restore_graph();

#ifdef PROF
            profiler_statistics.time_graph_restore += omp_get_wtime() - t;
#endif
        }

#ifdef PROF
        profiler_statistics.avg_path_size = this->get_avg_path_size();
#endif

        return paths;
    }

    std::vector<Path> compute_pd(const NODE_ID source, const NODE_ID destination, bool run_reverse_graph=true)
    {
#if !KSP_PD_L1
        omp_set_nested(1);
#endif

        std::vector<std::unique_ptr<DS>> d;
        d.reserve(num_threads);
#if KSP_PD_L1
        d.push_back(std::make_unique<DS>(graphs[0], num_threads));
        for(unsigned int i = 1; i < num_threads; i++)
            d.push_back(std::make_unique<DS>(graphs[i], 1));
#else
        for(unsigned int i = 0; i < num_threads; i++)
            d.push_back(std::make_unique<DS>(graphs[i], num_threads));
#endif

        if(run_reverse_graph){
            std::unique_ptr<DS> d_rev(nullptr);
            if(GraphType::is_directed)
                d_rev = std::make_unique<DS>(*reverse_graph, num_threads);
            else
                d_rev.swap(d[0]);

            // compute reverse sssp tree
            d_rev->reset_num_partition_threads(num_threads);
#ifdef PROF
            d_rev->template compute<false>(profiler_statistics, destination, source);
#else
            d_rev->template compute<false>(destination, source);
#endif
            sssp_t = d_rev->get_sssp_tree();

            candidates.insert({d_rev->get_distance(source), 0, 0, d_rev->get_shortest_path_to(source, true)});

            if(!GraphType::is_directed)
                d_rev.swap(d[0]);
        }

#if KSP_PD_L1
        d[0]->reset_num_partition_threads(1);
#else
        std::vector<unsigned int> threads(num_threads);
#endif
        std::vector<bool> simple_paths(num_threads, false);
        std::vector<unsigned int> curr_alpha(num_threads);

        for(unsigned int kp = 0; kp < k && !candidates.empty(); kp++)
        {
            paths.push_back(candidates.get_best());

            if(kp == k - 1)
                break;

            auto int_paths = static_cast<unsigned int>(paths[kp].p.size() - paths[kp].alpha - 1);

            if(int_paths == 0)
                continue;

            std::vector<Path> candidates_buffer(int_paths);
#if !KSP_PD_L1
            std::fill(threads.begin(), threads.end(), 1);
#endif
            const size_t paral_paths_1 = int_paths > num_threads ? num_threads : int_paths;

            int_paths = int_paths - paral_paths_1;
            unsigned int paths_done = 0;

#pragma omp parallel for num_threads(paral_paths_1)   // compute_pd()
            for(unsigned int inter_p = paths[kp].alpha; inter_p < paral_paths_1 + paths[kp].alpha; inter_p++)
            {
                const auto t_id = static_cast<unsigned int>(omp_get_thread_num());
                simple_paths[t_id] = find_simple_paths(kp, inter_p, t_id, candidates_buffer[inter_p - paths[kp].alpha], destination);
                curr_alpha[t_id] = inter_p;
            }

            paths_done += paral_paths_1 + paths[kp].alpha;
            std::vector<unsigned int> full_paths;
            for(unsigned int s = 0; s < paral_paths_1; s++)
            {
                if(!simple_paths[s])
                    full_paths.push_back(s);
            }

            std::fill(simple_paths.begin(), simple_paths.end(), false);

            const size_t full_path_size = full_paths.size();
            if(full_path_size > 0)
            {
#if !KSP_PD_L1
                const auto num_t = static_cast<unsigned int>(num_threads / full_path_size);
                for(unsigned int i = 0; i < full_path_size; i++)
                    threads[full_paths[i]] = num_t;

                const size_t left = num_threads - (num_t * full_path_size);
                for(size_t i = left; i > 0; i--)
                    threads[full_paths[i]]++;
#endif
                // compute full paths
#pragma omp parallel for num_threads(full_path_size)    // compute_pd()
                for(unsigned int i = 0; i < full_path_size; i++)
                {
                    const auto t_id = full_paths[i];
                    find_full_paths(*d[t_id], graphs[t_id], destination
#if !KSP_PD_L1
                                , threads[t_id]
#endif
                                , kp, curr_alpha[t_id], candidates_buffer[curr_alpha[t_id] - paths[kp].alpha]);
                }
#if !KSP_PD_L1
                std::fill(threads.begin(), threads.end(), 1);
#endif
            }

#pragma omp parallel for num_threads(paral_paths_1)   // compute_pd()
            for(size_t i = 0; i < paral_paths_1; i++)
                graphs[i].restore_graph();

            while(int_paths > 0)
            {
                const size_t paral_paths_2 = int_paths > num_threads ? num_threads : int_paths;

                int_paths = int_paths - paral_paths_2;

#pragma omp parallel for num_threads(paral_paths_2)   // compute_pd()
                for(unsigned int inter_p = paths_done; inter_p < paths_done + paral_paths_2; inter_p++)
                {
                    simple_paths[omp_get_thread_num()] = find_simple_paths(kp, inter_p, static_cast<unsigned int>(omp_get_thread_num()),
                                                                           candidates_buffer[inter_p - paths[kp].alpha], destination);
                    curr_alpha[omp_get_thread_num()] = inter_p;
                }

                paths_done += paral_paths_2;
                full_paths.clear();
                for(unsigned int s = 0; s < paral_paths_2; s++)
                {
                    if(!simple_paths[s])
                        full_paths.push_back(s);
                }

                std::fill(simple_paths.begin(), simple_paths.end(), false);

                if(!full_paths.empty())
                {
                    const auto full_paths_size = full_paths.size();
#if !KSP_PD_L1
                    const auto num_t = static_cast<unsigned int>(num_threads / full_paths_size);
                    for(unsigned int full_path : full_paths)
                        threads[full_path] = num_t;

                    const size_t left = num_threads - (num_t * full_paths_size);
                    for(size_t i = left; i > 0; i--)
                        threads[full_paths[i]]++;
#endif
                    // compute full paths
#pragma omp parallel for num_threads(full_paths_size)   // compute_pd()
                    for(unsigned int i = 0; i < full_paths_size; i++)
                    {
                        const auto t_id = full_paths[i];
                        find_full_paths(*d[t_id], graphs[t_id], destination
#if !KSP_PD_L1
                            , threads[t_id]
#endif
                            , kp, curr_alpha[t_id], candidates_buffer[curr_alpha[t_id] - paths[kp].alpha]);
                    }
                }

#pragma omp parallel for num_threads(paral_paths_2)   // compute_pd()
                for(size_t i = 0; i < paral_paths_2; i++)
                    graphs[i].restore_graph();
            }

            candidates.insert(candidates_buffer);
        }
#ifdef PROF
        profiler_statistics.avg_path_size = this->get_avg_path_size();
#endif

        return paths;
    }

    bool has_loops(std::vector<NODE_ID> temp_path) const noexcept // todo can this be done faster?
    {
        std::sort(temp_path.begin(), temp_path.end());
        return std::adjacent_find(temp_path.begin(), temp_path.end()) != temp_path.end();
    }

    Path compute_simple_path(const KSPGraph<GraphType>& this_g, const NODE_ID dest, const std::vector<NODE_ID>& prefix,
                             const unsigned int prefix_size, const std::vector<NODE_ID>& nodes_to_exclude) const noexcept
    {
        const NODE_ID branching_node = prefix[prefix_size - 1];
        w_type len = INFINITE_DISTANCE;
        NODE_ID ngh_node = NULL_NODE;
        
        
        EDGE_ID end_edge = this_g.end(branching_node);
        for(auto neighor = this_g.begin(branching_node); neighor < end_edge; ++neighor)
        {
            if(this_g.get_removed_edge()[neighor])
            {
                continue;
            }
            
            const auto csr = this_g.get_original_graph().get_csr_graph();
            const NODE_ID target = csr->adj[neighor];
            const w_type weight = csr->value[neighor];
            
            if(contains(nodes_to_exclude, target))
                continue;

            const w_type temp_len = sssp_t.at(target).tent + weight;
            if(temp_len < len)
            {
                len = temp_len;
                ngh_node = target;
            }
        }

        Path p;

        if(len == INFINITE_DISTANCE || len > candidates.get_length_threshold())
            return p;

        p.p.insert(p.p.end(), prefix.begin(), prefix.begin() + prefix_size);
        assert(ngh_node != NULL_NODE);
        p.p.push_back(ngh_node);

        while(ngh_node != dest)
        {
            ngh_node = sssp_t.at(ngh_node).parent;
            p.p.push_back(ngh_node);
        }

        if(!has_loops(p.p))
        {
            p.length = orig.get_original_graph().get_path_length(p.p);
            p.alpha = prefix_size - 1;
        }

        return p;
    }

    bool find_simple_paths(const unsigned int kp, const unsigned int alpha, const unsigned int t_id, Path& candidate_buffer, const NODE_ID destination) noexcept
    {
        std::vector<NODE_ID> avoid;

        if(alpha == paths[kp].alpha)
            avoid = remove_nodes(kp, t_id, alpha);
        else
            remove_nodes(kp, t_id, alpha);

        // remove previous node from graph
        for(unsigned int i = paths[kp].alpha; i < alpha; i++)
        {
            graphs[t_id].remove_node(paths[kp].p[i]);
        }

        // now compute intermediate paths
        avoid.push_back(paths[kp].p[alpha + 1]);
        candidate_buffer = compute_simple_path(graphs[t_id], destination, paths[kp].p, alpha + 1, avoid);

        if(candidate_buffer.length == INFINITE_DISTANCE || candidate_buffer.length > candidates.get_length_threshold())
            return false;

#if defined(PROF) && !defined(KSP_PARALLEL_DEVIATIONS_L2)
        profiler_statistics.num_simple_paths++;
#endif

        candidate_buffer.parent_index = candidate_buffer.alpha == paths[kp].alpha
                                        ? paths[kp].parent_index
                                        : kp;

        return true;
    }

    void find_full_paths(DS& d, KSPGraph<GraphType>& g, const NODE_ID dest
#if !KSP_PD_L1
                         , const unsigned int new_num_partition_threads
#endif
                         , const unsigned int parent_path_id, const unsigned int alpha, Path& candidate_buffer) noexcept
    {
        Path& parent_path = paths[parent_path_id];

#if !KSP_PD_L1
        d.reset_num_partition_threads(new_num_partition_threads);
#endif
        const NODE_ID node = parent_path.p[alpha];
        remove_undirected_edge(g, node, parent_path.p[alpha + 1]);

#ifdef TECHNIQUE_BOUND
        std::vector<NODE_ID> prefix_path = parent_path.getSubPath(alpha+1);
        w_type prefix_path_length=g.get_original_graph().get_path_length(prefix_path);
        w_type k_bound=candidates.getKBound();
        w_type sub_k_bound=k_bound-prefix_path_length;

        d.setKBound(sub_k_bound);
#endif


        d.template compute<SSSP_EARLY_STOPPING>(
#ifdef PROF
            profiler_statistics,
#endif
            node, dest, candidates.get_length_threshold() - orig.get_original_graph().get_path_length(parent_path.p, alpha));

#ifdef TECHNIQUE_BOUND
        if (sub_k_bound < d.get_distance(dest))
        {
            // it could get a path which is not the shortest if early terminate with TECHNIQUE_BOUND
            return;
        }
#endif

        if(d.get_distance(dest) < candidates.get_length_threshold())
        {
            std::vector<NODE_ID> new_path;
            new_path.insert(new_path.end(), parent_path.p.begin(), parent_path.p.begin() + alpha);

            const std::vector<NODE_ID> new_path_end = d.get_shortest_path_to(dest);
            new_path.insert(new_path.end(), new_path_end.begin(), new_path_end.end());
            candidate_buffer = {orig.get_original_graph().get_path_length(new_path),
                                                              parent_path.alpha == alpha
                                                                ? parent_path.parent_index
                                                                : parent_path_id,
                                                              alpha, new_path};
        }
    }

    bool contains(const std::vector<NODE_ID>& nodes, const NODE_ID node) const noexcept
    {
        return std::find(nodes.begin(), nodes.end(), node) != nodes.end();
    }

    void remove_edge(const NODE_ID src, const NODE_ID dest) noexcept
    {
        EDGE_ID edge_id = orig.get_edge_id(src, dest);

        if(edge_id != orig.end(src))
            orig.remove_edge(edge_id);
    }

    std::vector<NODE_ID> remove_nodes(const unsigned int kp, const unsigned int t_id, const unsigned int deviating_at) noexcept
    {
        std::vector<NODE_ID> avoid;
        if(paths[kp].alpha > 0)
        {
            for(unsigned int i = 0; i < paths[kp].alpha; i++)
            {
                graphs[t_id].remove_node(paths[kp].p[i]);
            }
        }

        for(unsigned int i = 0; i < kp; i++)
        {
            if(i == paths[kp].parent_index)
            {
                remove_undirected_edge(graphs[t_id], paths[i].p[paths[kp].alpha], paths[i].p[paths[kp].alpha + 1]);
                avoid.push_back(paths[i].p[paths[kp].alpha + 1]);
            }

            if(paths[i].parent_index == paths[kp].parent_index && paths[i].alpha == paths[kp].alpha && paths[kp].alpha >= deviating_at)
            {
                remove_undirected_edge(graphs[t_id], paths[i].p[paths[i].alpha], paths[i].p[paths[i].alpha + 1]);
                avoid.push_back(paths[i].p[paths[i].alpha + 1]);
            }
        }

        return avoid;
    }

    void remove_undirected_edge(KSPGraph<GraphType>& this_g, const NODE_ID u, const NODE_ID v) const noexcept
    {
        EDGE_ID edge_id = this_g.get_edge_id(u, v);

        if(edge_id != this_g.end(u))
            this_g.remove_edge(edge_id);
    }
};



template<template<typename, template<typename> class> class DeltaSteppingType, typename GraphType, unsigned int num_threads, bool pps = false>
class OptYenEdgeSwap : public KSPEdgeSwap<GraphType, num_threads, pps>
{
    using ParentType = KSPEdgeSwap<GraphType, num_threads, pps>;
    using DS = DeltaSteppingType<GraphType, KSPGraphEdgeSwap>;
    using ParentType::parallel_ps;
    using ParentType::orig;
    using ParentType::paths;
    using ParentType::candidates;
    using ParentType::k;
#ifdef PROF
    using ParentType::profiler_statistics;
#endif

    std::vector<KSPGraphEdgeSwap<GraphType>> graphs;
    SsspTree sssp_t;

public:
    explicit OptYenEdgeSwap(const GraphType& g, const unsigned int k, const SsspTree& sssp_t, std::vector<Path> path, std::unique_ptr<KSPGraph<GraphType>>& reverse_graph, bool* is_k_bound_nodes) noexcept : ParentType(g, k), sssp_t(sssp_t)
    {
        orig.init_edge_swap_with_2_pointers(is_k_bound_nodes, num_threads);

        candidates.insert(path);

        if(parallel_ps)
        {
            graphs.reserve(num_threads);
            for(unsigned int i = 1; i < num_threads; i++)
                graphs.push_back(KSPGraphEdgeSwap<GraphType>(g));

            graphs.push_back(orig);
        }
    }

    std::vector<Path> compute(const NODE_ID source, const NODE_ID destination)
    {
        assert(k > 1);

        paths.reserve(k);

        if(parallel_ps)
            return compute_pd(source, destination);

        return compute_psp(source, destination);
    }

private:
    std::vector<Path> compute_psp(const NODE_ID source, const NODE_ID destination)
    {
#ifdef PROF
        profiler_statistics.num_sssp++;
        double start = omp_get_wtime();
        double c_start;
        const double construct_time = omp_get_wtime();
#endif
        std::unique_ptr<DS> d = std::make_unique<DS>(orig, num_threads);

        for(unsigned int kp = 0; kp < k && !candidates.empty(); kp++)
        {
            paths.push_back(candidates.get_best());

            if(kp == k - 1)
                break;

            if(paths[kp].alpha > 0)
            {
                for(unsigned int i = 0; i < paths[kp].alpha; i++)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += orig.get_num_neighbors(paths[kp].p[i]) * (1 + static_cast<unsigned int>(!GraphType::is_directed));
                    const double rm_node_start = omp_get_wtime();
#endif

                    orig.remove_node(paths[kp].p[i]);

#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - rm_node_start;
#endif
                }
            }

            std::vector<NODE_ID> avoid;

            for(unsigned int i = 0; i < kp; i++)
            {
                if(i == paths[kp].parent_index)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    start = omp_get_wtime();
#endif

                    remove_edge(paths[i].p[paths[kp].alpha], paths[i].p[paths[kp].alpha + 1]);

#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - start;
#endif

                    avoid.push_back(paths[i].p[paths[kp].alpha + 1]);
                }

                if(paths[i].parent_index == paths[kp].parent_index && paths[i].alpha == paths[kp].alpha)
                {
#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    start = omp_get_wtime();
#endif
                    remove_edge(paths[i].p[paths[i].alpha], paths[i].p[paths[i].alpha + 1]);
#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - start;
#endif
                    avoid.push_back(paths[i].p[paths[i].alpha + 1]);
                }
            }

            for(unsigned int i = paths[kp].alpha; i < paths[kp].p.size() - 1; i++)
            {
                const NODE_ID node = paths[kp].p[i];


                avoid.push_back(paths[kp].p[i + 1]);
                Path p = compute_simple_path(orig, destination, paths[kp].p, i + 1, avoid);

                if(p.length < INFINITE_DISTANCE)
                {
#ifdef PROF
                    profiler_statistics.num_simple_paths++;
#endif
                }
                else
                {
                    p.p.clear();
                    
#ifdef TECHNIQUE_BOUND
                    std::vector<NODE_ID> prefix_path = paths[kp].getSubPath(i+1);
                    w_type prefix_path_length=orig.get_original_graph().get_path_length(prefix_path);
                    w_type k_bound=candidates.getKBound();
                    w_type sub_k_bound=k_bound-prefix_path_length;

                    d->setKBound(sub_k_bound);
#endif

#ifdef PROF
                    profiler_statistics.num_removed_edges += 2;
                    start = omp_get_wtime();
#endif

                    remove_edge(node, paths[kp].p[i + 1]);

#ifdef PROF
                    profiler_statistics.time_remove_edges += omp_get_wtime() - start;
                    const double t1 = omp_get_wtime();
#endif

                    d->template compute<SSSP_EARLY_STOPPING>(
#ifdef PROF
                        profiler_statistics,
#endif
                        node, destination, candidates.get_length_threshold() -  - orig.get_original_graph().get_path_length(paths[kp].p, i)
                    );

#ifdef PROF
                    profiler_statistics.time_sssp_compute += omp_get_wtime() - t1;
                    profiler_statistics.time_sssp_total += omp_get_wtime() - t1;
                    profiler_statistics.num_sssp++;
#endif

#ifdef TECHNIQUE_BOUND
                    if (sub_k_bound < d->get_distance(destination))
                    {
                        // it could get a path which is not the shortest if early terminate with TECHNIQUE_BOUND
                        orig.remove_node(node);
                        avoid.clear();
                        continue;
                    }
#endif

                    // concatenate paths
                    if(d->get_distance(destination) < candidates.get_length_threshold())
                    {
                        p.p.insert(p.p.end(), paths[kp].p.begin(), paths[kp].p.begin() + i);

                        const std::vector<NODE_ID> new_path_end = d->get_shortest_path_to(destination);
                        p.p.insert(p.p.end(), new_path_end.begin(), new_path_end.end());
                        p.length = orig.get_original_graph().get_path_length(p.p);
                        p.alpha = i;
                    }
                }

                if(p.length < candidates.get_length_threshold())
                {
                    p.parent_index = p.alpha == paths[kp].alpha ? paths[kp].parent_index : kp;
                    candidates.insert(std::move(p));
                }

#ifdef PROF
                profiler_statistics.num_removed_edges += orig.get_num_neighbors(node) * (1 + static_cast<unsigned int>(!GraphType::is_directed));
                start = omp_get_wtime();
#endif

                orig.remove_node(node);

#ifdef PROF
                profiler_statistics.time_remove_edges += omp_get_wtime() - start;
#endif

                avoid.clear();
            }
            // the end of intermediate paths

#ifdef PROF
            const double t = omp_get_wtime();
#endif

            orig.restore_graph();

#ifdef PROF
            profiler_statistics.time_graph_restore += omp_get_wtime() - t;
#endif
        }

#ifdef PROF
        profiler_statistics.avg_path_size = this->get_avg_path_size();
#endif

        return paths;
    }

    std::vector<Path> compute_pd(const NODE_ID source, const NODE_ID destination)
    {
#if !KSP_PD_L1
        omp_set_nested(1);
#endif

        std::vector<std::unique_ptr<DS>> d;
        d.reserve(num_threads);
#if KSP_PD_L1
        d.push_back(std::make_unique<DS>(graphs[0], num_threads));
        for(unsigned int i = 1; i < num_threads; i++)
            d.push_back(std::make_unique<DS>(graphs[i], 1));
#else
        for(unsigned int i = 0; i < num_threads; i++)
            d.push_back(std::make_unique<DS>(graphs[i], num_threads));
#endif

#if KSP_PD_L1
        d[0]->reset_num_partition_threads(1);
#else
        std::vector<unsigned int> threads(num_threads);
#endif
        std::vector<bool> simple_paths(num_threads, false);
        std::vector<unsigned int> curr_alpha(num_threads);

        for(unsigned int kp = 0; kp < k && !candidates.empty(); kp++)
        {
            paths.push_back(candidates.get_best());

            if(kp == k - 1)
                break;

            auto int_paths = static_cast<unsigned int>(paths[kp].p.size() - paths[kp].alpha - 1);

            if(int_paths == 0)
                continue;

            std::vector<Path> candidates_buffer(int_paths);
#if !KSP_PD_L1
            std::fill(threads.begin(), threads.end(), 1);
#endif
            const size_t paral_paths_1 = int_paths > num_threads ? num_threads : int_paths;

            int_paths = int_paths - paral_paths_1;
            unsigned int paths_done = 0;

#pragma omp parallel for num_threads(paral_paths_1)   // compute_pd()
            for(unsigned int inter_p = paths[kp].alpha; inter_p < paral_paths_1 + paths[kp].alpha; inter_p++)
            {
                const auto t_id = static_cast<unsigned int>(omp_get_thread_num());
                simple_paths[t_id] = find_simple_paths(kp, inter_p, t_id, candidates_buffer[inter_p - paths[kp].alpha], destination);
                curr_alpha[t_id] = inter_p;
            }

            paths_done += paral_paths_1 + paths[kp].alpha;
            std::vector<unsigned int> full_paths;
            for(unsigned int s = 0; s < paral_paths_1; s++)
            {
                if(!simple_paths[s])
                    full_paths.push_back(s);
            }

            std::fill(simple_paths.begin(), simple_paths.end(), false);

            const size_t full_path_size = full_paths.size();
            if(full_path_size > 0)
            {
#if !KSP_PD_L1
                const auto num_t = static_cast<unsigned int>(num_threads / full_path_size);
                for(unsigned int i = 0; i < full_path_size; i++)
                    threads[full_paths[i]] = num_t;

                const size_t left = num_threads - (num_t * full_path_size);
                for(size_t i = left; i > 0; i--)
                    threads[full_paths[i]]++;
#endif
                // compute full paths
#pragma omp parallel for num_threads(full_path_size)    // compute_pd()
                for(unsigned int i = 0; i < full_path_size; i++)
                {
                    const auto t_id = full_paths[i];
                    find_full_paths(*d[t_id], graphs[t_id], destination
#if !KSP_PD_L1
                                , threads[t_id]
#endif
                                , kp, curr_alpha[t_id], candidates_buffer[curr_alpha[t_id] - paths[kp].alpha]);
                }
#if !KSP_PD_L1
                std::fill(threads.begin(), threads.end(), 1);
#endif
            }

#pragma omp parallel for num_threads(paral_paths_1)   // compute_pd()
            for(size_t i = 0; i < paral_paths_1; i++)
                graphs[i].restore_graph();

            while(int_paths > 0)
            {
                const size_t paral_paths_2 = int_paths > num_threads ? num_threads : int_paths;

                int_paths = int_paths - paral_paths_2;

#pragma omp parallel for num_threads(paral_paths_2)   // compute_pd()
                for(unsigned int inter_p = paths_done; inter_p < paths_done + paral_paths_2; inter_p++)
                {
                    simple_paths[omp_get_thread_num()] = find_simple_paths(kp, inter_p, static_cast<unsigned int>(omp_get_thread_num()),
                                                                           candidates_buffer[inter_p - paths[kp].alpha], destination);
                    curr_alpha[omp_get_thread_num()] = inter_p;
                }

                paths_done += paral_paths_2;
                full_paths.clear();
                for(unsigned int s = 0; s < paral_paths_2; s++)
                {
                    if(!simple_paths[s])
                        full_paths.push_back(s);
                }

                std::fill(simple_paths.begin(), simple_paths.end(), false);

                if(!full_paths.empty())
                {
                    const auto full_paths_size = full_paths.size();
#if !KSP_PD_L1
                    const auto num_t = static_cast<unsigned int>(num_threads / full_paths_size);
                    for(unsigned int full_path : full_paths)
                        threads[full_path] = num_t;

                    const size_t left = num_threads - (num_t * full_paths_size);
                    for(size_t i = left; i > 0; i--)
                        threads[full_paths[i]]++;
#endif
                    // compute full paths
#pragma omp parallel for num_threads(full_paths_size)   // compute_pd()
                    for(unsigned int i = 0; i < full_paths_size; i++)
                    {
                        const auto t_id = full_paths[i];
                        find_full_paths(*d[t_id], graphs[t_id], destination
#if !KSP_PD_L1
                            , threads[t_id]
#endif
                            , kp, curr_alpha[t_id], candidates_buffer[curr_alpha[t_id] - paths[kp].alpha]);
                    }
                }

#pragma omp parallel for num_threads(paral_paths_2)   // compute_pd()
                for(size_t i = 0; i < paral_paths_2; i++)
                    graphs[i].restore_graph();
            }

            candidates.insert(candidates_buffer);
        }
#ifdef PROF
        profiler_statistics.avg_path_size = this->get_avg_path_size();
#endif

        return paths;
    }

    bool has_loops(std::vector<NODE_ID> temp_path) const noexcept // todo can this be done faster?
    {
        std::sort(temp_path.begin(), temp_path.end());
        return std::adjacent_find(temp_path.begin(), temp_path.end()) != temp_path.end();
    }

    Path compute_simple_path(const KSPGraphEdgeSwap<GraphType>& this_g, const NODE_ID dest, const std::vector<NODE_ID>& prefix,
                             const unsigned int prefix_size, const std::vector<NODE_ID>& nodes_to_exclude) const noexcept
    {
        const NODE_ID branching_node = prefix[prefix_size - 1];
        w_type len = INFINITE_DISTANCE;
        NODE_ID ngh_node = NULL_NODE;
        
        
        const auto csr = this_g.get_original_graph().get_csr_graph();
        auto node_end = this_g.end_light(branching_node);
        for(EDGE_ID neighor = this_g.begin_light(branching_node); neighor < node_end; neighor++)
        {
            if(this_g.get_removed_edge()[neighor])
            {
                continue;
            }

            const NODE_ID target = csr->adj[neighor];
            const w_type weight = csr->value[neighor];
            
            if(contains(nodes_to_exclude, target))
                continue;

            const w_type temp_len = sssp_t.at(target).tent + weight;
            if(temp_len < len)
            {
                len = temp_len;
                ngh_node = target;
            }
        }

        node_end = this_g.end_heavy(branching_node);
        for(EDGE_ID neighor = this_g.begin_heavy(branching_node); neighor < node_end; neighor++)
        {
            if(this_g.get_removed_edge()[neighor])
            {
                continue;
            }

            const NODE_ID target = csr->adj[neighor];
            const w_type weight = csr->value[neighor];
            
            if(contains(nodes_to_exclude, target))
                continue;

            const w_type temp_len = sssp_t.at(target).tent + weight;
            if(temp_len < len)
            {
                len = temp_len;
                ngh_node = target;
            }
        }

        Path p;

        if(len == INFINITE_DISTANCE || len > candidates.get_length_threshold())
            return p;

        p.p.insert(p.p.end(), prefix.begin(), prefix.begin() + prefix_size);
        assert(ngh_node != NULL_NODE);
        p.p.push_back(ngh_node);

        while(ngh_node != dest)
        {
            ngh_node = sssp_t.at(ngh_node).parent;
            p.p.push_back(ngh_node);
        }

        if(!has_loops(p.p))
        {
            p.length = orig.get_original_graph().get_path_length(p.p);
            p.alpha = prefix_size - 1;
        }

        return p;
    }

    bool find_simple_paths(const unsigned int kp, const unsigned int alpha, const unsigned int t_id, Path& candidate_buffer, const NODE_ID destination) noexcept
    {
        std::vector<NODE_ID> avoid;

        if(alpha == paths[kp].alpha)
            avoid = remove_nodes(kp, t_id, alpha);
        else
            remove_nodes(kp, t_id, alpha);

        // remove previous node from graph
        for(unsigned int i = paths[kp].alpha; i < alpha; i++)
        {
            graphs[t_id].remove_node(paths[kp].p[i]);
        }

        // now compute intermediate paths
        avoid.push_back(paths[kp].p[alpha + 1]);
        candidate_buffer = compute_simple_path(graphs[t_id], destination, paths[kp].p, alpha + 1, avoid);

        if(candidate_buffer.length == INFINITE_DISTANCE || candidate_buffer.length > candidates.get_length_threshold())
            return false;

#if defined(PROF) && !defined(KSP_PARALLEL_DEVIATIONS_L2)
        profiler_statistics.num_simple_paths++;
#endif

        candidate_buffer.parent_index = candidate_buffer.alpha == paths[kp].alpha
                                        ? paths[kp].parent_index
                                        : kp;

        return true;
    }

    void find_full_paths(DS& d, KSPGraphEdgeSwap<GraphType>& g, const NODE_ID dest
#if !KSP_PD_L1
                         , const unsigned int new_num_partition_threads
#endif
                         , const unsigned int parent_path_id, const unsigned int alpha, Path& candidate_buffer) noexcept
    {
        Path& parent_path = paths[parent_path_id];

#if !KSP_PD_L1
        d.reset_num_partition_threads(new_num_partition_threads);
#endif
        const NODE_ID node = parent_path.p[alpha];
        remove_undirected_edge(g, node, parent_path.p[alpha + 1]);

#ifdef TECHNIQUE_BOUND
        std::vector<NODE_ID> prefix_path = parent_path.getSubPath(alpha+1);
        w_type prefix_path_length=g.get_original_graph().get_path_length(prefix_path);
        w_type k_bound=candidates.getKBound();
        w_type sub_k_bound=k_bound-prefix_path_length;

        d.setKBound(sub_k_bound);
#endif


        d.template compute<SSSP_EARLY_STOPPING>(
#ifdef PROF
            profiler_statistics,
#endif
            node, dest, candidates.get_length_threshold() - orig.get_original_graph().get_path_length(parent_path.p, alpha));

#ifdef TECHNIQUE_BOUND
        if (sub_k_bound < d.get_distance(dest))
        {
            // it could get a path which is not the shortest if early terminate with TECHNIQUE_BOUND
            return;
        }
#endif

        if(d.get_distance(dest) < candidates.get_length_threshold())
        {
            std::vector<NODE_ID> new_path;
            new_path.insert(new_path.end(), parent_path.p.begin(), parent_path.p.begin() + alpha);

            const std::vector<NODE_ID> new_path_end = d.get_shortest_path_to(dest);
            new_path.insert(new_path.end(), new_path_end.begin(), new_path_end.end());
            candidate_buffer = {orig.get_original_graph().get_path_length(new_path),
                                                              parent_path.alpha == alpha
                                                                ? parent_path.parent_index
                                                                : parent_path_id,
                                                              alpha, new_path};
        }
    }

    bool contains(const std::vector<NODE_ID>& nodes, const NODE_ID node) const noexcept
    {
        return std::find(nodes.begin(), nodes.end(), node) != nodes.end();
    }

    void remove_edge(const NODE_ID src, const NODE_ID dest) noexcept
    {
        orig.remove_edge(src, dest);
    }

    std::vector<NODE_ID> remove_nodes(const unsigned int kp, const unsigned int t_id, const unsigned int deviating_at) noexcept
    {
        std::vector<NODE_ID> avoid;
        if(paths[kp].alpha > 0)
        {
            for(unsigned int i = 0; i < paths[kp].alpha; i++)
            {
                graphs[t_id].remove_node(paths[kp].p[i]);
            }
        }

        for(unsigned int i = 0; i < kp; i++)
        {
            if(i == paths[kp].parent_index)
            {
                remove_undirected_edge(graphs[t_id], paths[i].p[paths[kp].alpha], paths[i].p[paths[kp].alpha + 1]);
                avoid.push_back(paths[i].p[paths[kp].alpha + 1]);
            }

            if(paths[i].parent_index == paths[kp].parent_index && paths[i].alpha == paths[kp].alpha && paths[kp].alpha >= deviating_at)
            {
                remove_undirected_edge(graphs[t_id], paths[i].p[paths[i].alpha], paths[i].p[paths[i].alpha + 1]);
                avoid.push_back(paths[i].p[paths[i].alpha + 1]);
            }
        }

        return avoid;
    }

    void remove_undirected_edge(KSPGraphEdgeSwap<GraphType>& this_g, const NODE_ID u, const NODE_ID v) const noexcept
    {
        this_g.remove_edge(u, v);
    }

};


#endif  // _OPT_YEN_H
