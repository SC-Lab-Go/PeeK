#ifndef _PEEK_ADAPTIVE_WITH_EDGE_SWAP_PARALLEL_H
#define _PEEK_ADAPTIVE_WITH_EDGE_SWAP_PARALLEL_H

#include <DeltaSteppingStatic.h>
#include <GraphRW.h>
#include <KSP.h>
#include <MaxHeapNoDuplicate.h>
#include <OptYen.h>
#include <sys/time.h>

#include <boost/sort/sort.hpp>
#include <map>

inline double wtime() {
    double time[2];
    struct timeval time1;
    gettimeofday(&time1, NULL);

    time[0] = time1.tv_sec;
    time[1] = time1.tv_usec;

    return time[0] + time[1] * 1.0e-6;
}

template <template <typename, template <typename> class> class DeltaSteppingType, typename GraphType, unsigned int num_threads, bool pps = false>
class PeeKAdaptiveWithEdgeSwap : public KSP<GraphType, num_threads, pps> {
    using ParentType = KSP<GraphType, num_threads, pps>;
    using DS = DeltaSteppingType<GraphType, KSPGraph>;
    using ParentType::candidates;
    using ParentType::k;
    using ParentType::orig;
    using ParentType::parallel_ps;
    using ParentType::paths;

    std::unique_ptr<KSPGraph<GraphType>> reverse_graph;
    w_type k_bound;
    SsspTree reverse_sssp_t;
    SsspTree normal_sssp_t;
    unsigned int actual_k;
    const double adaptive_ratio = 0.1;  // adaptive graph compact ratio is 10% of total nodes, use regeneration when less, use status array or edge swap when great
    bool* is_k_bound_nodes;
    NODE_ID* original_2_compact;
    unsigned int k_bound_nodes_num;

   public:
    explicit PeeKAdaptiveWithEdgeSwap(const GraphType& g, const unsigned int k) noexcept : ParentType(g, k), reverse_graph(GraphType::is_directed ? std::make_unique<KSPGraph<GraphType>>(g.get_reverse_graph()) : nullptr), reverse_sssp_t(g.get_num_nodes()), normal_sssp_t(g.get_num_nodes()) {
        assert(GraphType::is_directed);
        assert(GraphType::is_weighted);
        k_bound = INFINITE_DISTANCE;
        actual_k = 0;
        is_k_bound_nodes = new bool[g.get_num_nodes()]();
        original_2_compact = new NODE_ID[g.get_num_nodes()]();
        k_bound_nodes_num = 0;
    }

    ~PeeKAdaptiveWithEdgeSwap() {
        delete[] is_k_bound_nodes;
        delete[] original_2_compact;
    }

    w_type getKBound() { return k_bound; }

    std::vector<Path> compute(const NODE_ID source, const NODE_ID destination) {
        double start_time, end_time;
        std::vector<Path> result;

        // start_time = wtime();
        reverseSSSP(source, destination);
        // std::cout << "," << wtime() - start_time;

        // start_time = wtime();
        normalSSSP(source, destination);
        // std::cout << "," << wtime() - start_time;

        // start_time = wtime();
        identifyLooplessKBound(source, destination);
        // std::cout << "," << wtime() - start_time;

        if (k_bound_nodes_num > adaptive_ratio * orig.get_num_nodes()) {
            // start_time = wtime();
            OptYenEdgeSwap<DeltaSteppingStatic, GraphType, num_threads, pps> optyen(orig.get_original_graph(), k, reverse_sssp_t, paths, reverse_graph, is_k_bound_nodes);

#ifdef TECHNIQUE_BOUND
            // float number could cause one less path in some cases, in order to make sure correct totally, now comment it out or make the k bound is greater
            if (actual_k == k) {
                optyen.setKBound(k_bound * 1.01);
            }
#endif

            result = optyen.compute(source, destination);
            // end_time = wtime();
            // optyen.print_paths();
            // std::cout << ",edge_swap:" << end_time - start_time;
            return result;
        }

        // start_time = wtime();
        const auto compact_g = regenerateGraph();

        SsspTree compact_reverse_sssp_t(compact_g.get_csr_graph()->num_nodes);
        regenerateGraphSsspTree(compact_reverse_sssp_t);

        // transfer original nodes into compact nodes among source destination and paths
        NODE_ID compact_source;
        NODE_ID compact_destination;
        if (is_k_bound_nodes[source] && is_k_bound_nodes[destination]) {
            compact_source = original_2_compact[source];
            compact_destination = original_2_compact[destination];
        } else {
            std::cerr << "cannot find compact_path_node in compact_graph" << std::endl;
            return result;
        }

        for (unsigned int i = 0; i < paths.size(); i++) {
            if (paths[i].length == INFINITE_DISTANCE) continue;

            for (unsigned int j = 0; j < paths[i].p.size(); j++) {
                if (is_k_bound_nodes[paths[i].p[j]]) {
                    paths[i].p[j] = original_2_compact[paths[i].p[j]];
                } else {
                    std::cerr << "cannot find compact_path_node in compact_graph" << std::endl;
                    return result;
                }
            }
        }
        // std::cout << "," << wtime() - start_time;

        // start_time = wtime();
        OptYen<DeltaSteppingStatic, GraphType, num_threads, pps> optyen(compact_g, k, compact_reverse_sssp_t, paths);

#ifdef TECHNIQUE_BOUND
        // float number could cause one less path in some cases, in order to make sure correct totally, now comment it out or make the k bound is greater
        if (actual_k == k) {
            optyen.setKBound(k_bound * 1.01);
        }
#endif

        result = optyen.compute(compact_source, compact_destination, false);
        // end_time = wtime();
        // optyen.print_compact_paths(original_2_compact, is_k_bound_nodes, orig.get_num_nodes());
        // std::cout << "," << end_time - start_time;

        compact_g.get_csr_graph()->destroy();

        return result;
    }

    void normalSSSP(const NODE_ID source, const NODE_ID destination) {
        const NODE_ID num_nodes = orig.get_num_nodes();
        auto csr = orig.get_original_graph().get_csr_graph();
#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID node = 0; node < num_nodes; node++) {
            if (INFINITE_DISTANCE == reverse_sssp_t[node].tent) {
                const EDGE_ID edge_list_end = orig.end(node);
                for (EDGE_ID neighor = orig.begin(node); neighor < edge_list_end; neighor++) {
                    orig.edge_is_removed[neighor] = true;
                }
            } else {
                const EDGE_ID edge_list_end = orig.end(node);
                for (EDGE_ID neighor = orig.begin(node); neighor < edge_list_end; neighor++) {
                    NODE_ID node_to = csr->adj[neighor];
                    if (INFINITE_DISTANCE == reverse_sssp_t[node_to].tent) {
                        orig.edge_is_removed[neighor] = true;
                    }
                }
            }
        }

        DS normal_ds(orig, num_threads);
        normal_ds.template compute_sssp(source);
        normal_sssp_t = normal_ds.get_sssp_tree();
    }

    void reverseSSSP(const NODE_ID source, const NODE_ID destination) {
        DS reverse_ds(*reverse_graph, num_threads);
        reverse_ds.template compute_sssp(destination);
        paths.push_back({reverse_ds.get_distance(source), 0, 0, reverse_ds.get_shortest_path_to(source, true)});
        reverse_sssp_t = reverse_ds.get_sssp_tree();
    }

    void identifyLooplessKBound(const NODE_ID source, const NODE_ID destination) {
        const NODE_ID num_nodes = orig.get_num_nodes();
        nodeWeight* node_weight = new nodeWeight[num_nodes];

        // sort distance
#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID i = 0; i < num_nodes; i++) {
            w_type weight = normal_sssp_t[i].tent + reverse_sssp_t[i].tent;
            node_weight[i].id = i;
            node_weight[i].distance = weight;
        }
        boost::sort::block_indirect_sort(node_weight, node_weight + num_nodes, cmp, num_threads);
        // std::sort(node_weight, node_weight + num_nodes, cmp);

        // calc k nodes
        NODE_ID i;
        w_type path_weight;
        for (i = 0; i < num_nodes; i++) {
            path_weight = node_weight[i].distance;
            if (path_weight == INFINITE_DISTANCE) {
                break;
            }

            // check if the shortest path is loop
            NODE_ID node = node_weight[i].id;
            if (is_k_bound_nodes[node]) {
                continue;
            }

            NODE_ID normal_node = node;
            std::vector<NODE_ID> normal_path;
            normal_path.clear();
            normal_path.push_back(node);
            NODE_ID normal_parent = normal_sssp_t[normal_node].parent;
            while (normal_parent != normal_node) {
                normal_node = normal_parent;
                normal_parent = normal_sssp_t[normal_node].parent;
                normal_path.push_back(normal_node);
            }

            NODE_ID reverse_node = node;
            std::vector<NODE_ID> reverse_path;
            reverse_path.clear();
            std::set<NODE_ID> reverse_path_set;
            NODE_ID reverse_parent = reverse_sssp_t[reverse_node].parent;
            while (reverse_parent != reverse_node) {
                reverse_node = reverse_parent;
                reverse_parent = reverse_sssp_t[reverse_node].parent;
                reverse_path_set.insert(reverse_node);
                reverse_path.push_back(reverse_node);
            }

            volatile bool is_loop = false;
#pragma omp parallel for num_threads(num_threads) shared(is_loop)
            for (NODE_ID i = 0; i < normal_path.size(); ++i) {
                if (!is_loop && reverse_path_set.find(normal_path[i]) != reverse_path_set.end()) {
                    is_loop = true;
                }
            }

            if (!is_loop) {
                // add simple path
#pragma omp parallel for num_threads(num_threads)
                for (NODE_ID i = 0; i < normal_path.size(); ++i) {
                    if (!is_k_bound_nodes[normal_path[i]]) {
                        is_k_bound_nodes[normal_path[i]] = true;
                    }
                }
#pragma omp parallel for num_threads(num_threads)
                for (NODE_ID i = 0; i < reverse_path.size(); ++i) {
                    if (!is_k_bound_nodes[reverse_path[i]]) {
                        is_k_bound_nodes[reverse_path[i]] = true;
                    }
                }

                // only calc k bound when path is loopless
                actual_k++;

                if (actual_k == k) {
                    // add extra same distance nodes
                    while (++i < num_nodes) {
                        NODE_ID new_node = node_weight[i].id;
                        w_type new_path_weight = node_weight[i].distance;
                        if (FLOATEQUAL(new_path_weight, path_weight)) {
                            if (!is_k_bound_nodes[new_node]) {
                                is_k_bound_nodes[new_node] = true;
                            }
                        } else {
                            break;
                        }
                    }
                    // get k difference loopless paths nodes
                    break;
                }
            } else {
                is_k_bound_nodes[node] = true;
            }
        }

        k_bound = path_weight;

        NODE_ID thread_k_bound_nodes[num_threads] = {0};
#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID node = 0; node < num_nodes; node++) {
            if (is_k_bound_nodes[node]) {
                const auto tid = static_cast<unsigned int>(omp_get_thread_num());
                thread_k_bound_nodes[tid]++;
            }
        }

        for (unsigned int i = 0; i < num_threads; i++) {
            k_bound_nodes_num += thread_k_bound_nodes[i];
        }

        delete[] node_weight;
    }

    void regenerateGraphSsspTree(SsspTree& compact_reverse_sssp_t) {
        const NODE_ID num_nodes = orig.get_num_nodes();
#pragma omp parallel for num_threads(num_threads)
        for (NODE_ID original_node = 0; original_node < num_nodes; original_node++) {
            if (is_k_bound_nodes[original_node]) {
                NODE_ID original_parent = reverse_sssp_t[original_node].parent;
                if (NULL_NODE == original_parent) {
                    NODE_ID compact_node = original_2_compact[original_node];
                    
                    compact_reverse_sssp_t[compact_node].tent = INFINITE_DISTANCE;
                    compact_reverse_sssp_t[compact_node].parent = NULL_NODE;
                } else if (is_k_bound_nodes[original_parent]) {
                    w_type tent = reverse_sssp_t[original_node].tent;
                    NODE_ID compact_node = original_2_compact[original_node];
                    NODE_ID compact_parent = original_2_compact[original_parent];

                    compact_reverse_sssp_t[compact_node].tent = tent;
                    compact_reverse_sssp_t[compact_node].parent = compact_parent;
                }
            }
        }
    }

    GraphType regenerateGraph() {
        const NODE_ID num_nodes = orig.get_num_nodes();
        const auto csr = orig.get_original_graph().get_csr_graph();
        NODE_ID thread_num_nodes = num_nodes / num_threads;
        NODE_ID thread_k_bound_nodes[num_threads] = {0};
        NODE_ID thread_k_bound_edges[num_threads] = {0};
        NODE_ID prefix_sum_node[num_threads] = {0};
        EDGE_ID prefix_sum_edge[num_threads] = {0};
        NODE_ID node_index_start[num_threads] = {0};
        NODE_ID node_index_end[num_threads] = {0};
        w_type thread_total_weight[num_threads] = {0.0};
        w_type thread_heaviest_weight[num_threads] = {0.0};
        EDGE_ID thread_max_out_degree[num_threads] = {0};

#pragma omp parallel num_threads(num_threads)
        {
            const auto tid = static_cast<unsigned int>(omp_get_thread_num());

            // step 1: get k_bound_node and k_bound_edge for each thread
            NODE_ID node_start = tid * thread_num_nodes;
            NODE_ID node_end = tid * thread_num_nodes + thread_num_nodes;
            if (tid == (num_threads - 1)) {
                node_end = num_nodes;
            }
            node_index_start[tid] = node_start;
            node_index_end[tid] = node_end;

            for (NODE_ID node = node_start; node < node_end; node++) {
                if (is_k_bound_nodes[node]) {
                    thread_k_bound_nodes[tid]++;
                    EDGE_ID out_degree = 0;
                    const EDGE_ID edge_list_end = orig.end(node);
                    for (EDGE_ID neighor = orig.begin(node); neighor < edge_list_end; neighor++) {
                        NODE_ID node_to = csr->adj[neighor];
                        w_type weight = csr->value[neighor];
                        if (is_k_bound_nodes[node_to]) {
                            if (actual_k == k && weight > k_bound) {
                                continue;
                            }
                            thread_k_bound_edges[tid]++;

                            thread_total_weight[tid] += weight;
                            if (weight > thread_heaviest_weight[tid]) {
                                thread_heaviest_weight[tid] = weight;
                            }
                            out_degree++;
                        }
                    }
                    if (out_degree > thread_max_out_degree[tid]) {
                        thread_max_out_degree[tid] = out_degree;
                    }
                }
            }

            // step 2: get prefix_sum_node and prefix_sum_edge for each thread
#pragma omp barrier
            for (unsigned int i = 0; i < tid; ++i) {
                prefix_sum_node[tid] += thread_k_bound_nodes[i];
                prefix_sum_edge[tid] += thread_k_bound_edges[i];
            }

            // step 3: map original nodes into compact nodes
#pragma omp barrier
            NODE_ID start_node_pos = prefix_sum_node[tid];
            for (NODE_ID node = node_start; node < node_end; node++) {
                if (is_k_bound_nodes[node]) {
                    original_2_compact[node] = start_node_pos++;
                }
            }
        }

        // step 4: allocate new memory
        NODE_ID compact_num_nodes = prefix_sum_node[num_threads - 1] + thread_k_bound_nodes[num_threads - 1];
        EDGE_ID compact_num_edges = prefix_sum_edge[num_threads - 1] + thread_k_bound_edges[num_threads - 1];

        w_type total_weight = 0.0;
        for (unsigned int i = 0; i < num_threads; i++) {
            total_weight += thread_total_weight[i];
        }

        w_type heaviest_weight = 0.0;
        for (unsigned int i = 0; i < num_threads; i++) {
            if (heaviest_weight < thread_heaviest_weight[i]) {
                heaviest_weight = thread_heaviest_weight[i];
            }
        }

        EDGE_ID max_out_degree = 0;
        for (unsigned int i = 0; i < num_threads; i++) {
            if (max_out_degree < thread_max_out_degree[i]) {
                max_out_degree = thread_max_out_degree[i];
            }
        }

        // make light and heavy neighor nodes together separately to use delta-steping algorithm
        NODE_ID* thread_light_node[num_threads];
        w_type* thread_light_weight[num_threads];
        NODE_ID thread_light_num[num_threads] = {0};
        NODE_ID* thread_heavy_node[num_threads];
        w_type* thread_heavy_weight[num_threads];
        NODE_ID thread_heavy_num[num_threads] = {0};
        for (unsigned int i = 0; i < num_threads; i++) {
            thread_light_node[i] = new NODE_ID[max_out_degree]();
            thread_light_weight[i] = new w_type[max_out_degree]();

            thread_heavy_node[i] = new NODE_ID[max_out_degree]();
            thread_heavy_weight[i] = new w_type[max_out_degree]();
        }

        w_type avg_degree = static_cast<w_type>(compact_num_edges) / static_cast<w_type>(compact_num_nodes);
        w_type avg_weight = total_weight / compact_num_edges;
        w_type dist = paths[0].length;
        int l = paths[0].p.size() - 1;
        w_type delta = (dist / l + avg_weight) / avg_degree;
        CSRGraph* csr_compact = new CSRGraph(compact_num_nodes, compact_num_edges);

#pragma omp parallel num_threads(num_threads)
        {
            const auto tid = static_cast<unsigned int>(omp_get_thread_num());

            // step 5: write compact graph data to compact_csr
            NODE_ID node_start = node_index_start[tid];
            NODE_ID node_end = node_index_end[tid];
            EDGE_ID start_edge_pos = prefix_sum_edge[tid];

            for (NODE_ID node = node_start; node < node_end; node++) {
                if (is_k_bound_nodes[node]) {
                    NODE_ID compact_node_from = original_2_compact[node];
                    csr_compact->begin[compact_node_from] = start_edge_pos;
                    thread_light_num[tid] = 0;
                    thread_heavy_num[tid] = 0;
                    const EDGE_ID edge_list_end = orig.end(node);
                    for (EDGE_ID neighor = orig.begin(node); neighor < edge_list_end; neighor++) {
                        NODE_ID node_to = csr->adj[neighor];
                        w_type weight = csr->value[neighor];
                        if (is_k_bound_nodes[node_to]) {
                            if (actual_k == k && weight > k_bound) {
                                continue;
                            }

                            NODE_ID compact_node_to = original_2_compact[node_to];
                            if (weight <= delta) {
                                thread_light_node[tid][thread_light_num[tid]] = compact_node_to;
                                thread_light_weight[tid][thread_light_num[tid]] = weight;
                                thread_light_num[tid]++;
                            } else {
                                thread_heavy_node[tid][thread_heavy_num[tid]] = compact_node_to;
                                thread_heavy_weight[tid][thread_heavy_num[tid]] = weight;
                                thread_heavy_num[tid]++;
                            }
                        }
                    }

                    csr_compact->num_light_edges[compact_node_from] = thread_light_num[tid];
                    for (NODE_ID i = 0; i < thread_light_num[tid]; i++, start_edge_pos++) {
                        csr_compact->adj[start_edge_pos] = thread_light_node[tid][i];
                        csr_compact->value[start_edge_pos] = thread_light_weight[tid][i];
                    }
                    for (NODE_ID i = 0; i < thread_heavy_num[tid]; i++, start_edge_pos++) {
                        csr_compact->adj[start_edge_pos] = thread_heavy_node[tid][i];
                        csr_compact->value[start_edge_pos] = thread_heavy_weight[tid][i];
                    }
                }
            }
        }

        csr_compact->begin[compact_num_nodes] = compact_num_edges;

        for (unsigned int i = 0; i < num_threads; i++) {
            delete[] thread_light_node[i];
            delete[] thread_light_weight[i];

            delete[] thread_heavy_node[i];
            delete[] thread_heavy_weight[i];
        }

        return GraphType(csr_compact, heaviest_weight, delta);
    }
};

#endif  // _PEEK_ADAPTIVE_WITH_EDGE_SWAP_PARALLEL_H
