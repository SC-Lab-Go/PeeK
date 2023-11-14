#ifndef _PEEK_ADAPTIVE_WITH_EDGE_SWAP_H
#define _PEEK_ADAPTIVE_WITH_EDGE_SWAP_H

#include <DeltaSteppingStatic.h>
#include <GraphRW.h>
#include <KSP.h>
#include <MaxHeapNoDuplicate.h>
#include <OptYen.h>
#include <sys/time.h>

#include <map>

#ifndef COMPACT_KSP_NUM_THREADS
#define COMPACT_KSP_NUM_THREADS 1
#endif
#define COMPACT_KSP_PD true

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
    std::vector<NODE_ID> reverse_unreachable_nodes;
    w_type k_bound;
    std::vector<NODE_ID> k_bound_nodes;            // k_bound_nodes[compact_node]=original_node
    std::map<NODE_ID, NODE_ID> k_bound_nodes_map;  // k_bound_nodes_map[original_node]=compact_node
    SsspTree reverse_sssp_t;
    SsspTree normal_sssp_t;
    unsigned int actual_k;
    const double adaptive_ratio = 0.1;  // adaptive graph compact ratio is 10% of total nodes, use regeneration when less, use status array or edge swap when great
    bool* is_k_bound_nodes;
    unsigned int k_bound_nodes_num;

   public:
    explicit PeeKAdaptiveWithEdgeSwap(const GraphType& g, const unsigned int k) noexcept : ParentType(g, k), reverse_graph(GraphType::is_directed ? std::make_unique<KSPGraph<GraphType>>(g.get_reverse_graph()) : nullptr), reverse_sssp_t(g.get_num_nodes()), normal_sssp_t(g.get_num_nodes()) {
        assert(GraphType::is_directed);
        assert(GraphType::is_weighted);
        k_bound = INFINITE_DISTANCE;
        actual_k = 0;
        is_k_bound_nodes = new bool[g.get_num_nodes()]();
        k_bound_nodes_num = 0;
    }

    ~PeeKAdaptiveWithEdgeSwap() { delete[] is_k_bound_nodes; }

    w_type getKBound() { return k_bound; }

    std::vector<Path> compute(const NODE_ID source, const NODE_ID destination) {
        double start_time, end_time;
        std::vector<Path> result;

        // start_time = wtime();
        reverseSSSP(source, destination);
        // std::cout << "reverseSSSP time:" << wtime() - start_time << std::endl;

        // start_time = wtime();
        normalSSSP(source, destination);
        // std::cout << "normalSSSP time:" << wtime() - start_time << std::endl;

        // start_time = wtime();
        identifyLooplessKBound(source, destination);
        // std::cout << "identifyLooplessKBound time:" << wtime() - start_time << std::endl;

        if (k_bound_nodes_num > adaptive_ratio * orig.get_num_nodes()) {
            OptYenEdgeSwap<DeltaSteppingStatic, GraphType, num_threads, COMPACT_KSP_PD> optyen(orig.get_original_graph(), k, reverse_sssp_t, paths, reverse_graph, is_k_bound_nodes);
#ifdef TECHNIQUE_BOUND
            // float number could cause one less path in some cases, in order to make sure correct totally, now comment it out or solution 2, make the k bound is greater
            if (actual_k == k) {
                optyen.setKBound(k_bound * 1.01);
            }
#endif
            result = optyen.compute(source, destination);
            // optyen.print_paths();
            return result;
        }

        // start_time = wtime();
        SsspTree compact_reverse_sssp_t(k_bound_nodes.size());
        regenerateGraphSsspTree(compact_reverse_sssp_t);
        // std::cout << "regenerateGraphSsspTree time:" << wtime() - start_time << std::endl;

        // start_time = wtime();
        const auto compact_g = regenerateGraph();
        // std::cout << "regenerateGraph time:" << wtime() - start_time << std::endl;

        // start_time = wtime();
        // transfer original nodes into compact nodes among source destination and paths
        NODE_ID compact_source;
        NODE_ID compact_destination;
        std::map<NODE_ID, NODE_ID>::iterator it_source = k_bound_nodes_map.find(source);
        std::map<NODE_ID, NODE_ID>::iterator it_destination = k_bound_nodes_map.find(destination);
        if (it_source != k_bound_nodes_map.end() && it_destination != k_bound_nodes_map.end()) {
            compact_source = it_source->second;
            compact_destination = it_destination->second;
        } else {
            std::cerr << "cannot find compact_source or compact_destination in compact_graph" << std::endl;
            return result;
        }

        for (unsigned int i = 0; i < paths.size(); i++) {
            if (paths[i].length == INFINITE_DISTANCE) continue;

            for (unsigned int j = 0; j < paths[i].p.size(); j++) {
                std::map<NODE_ID, NODE_ID>::iterator it = k_bound_nodes_map.find(paths[i].p[j]);
                if (it != k_bound_nodes_map.end()) {
                    paths[i].p[j] = it->second;
                } else {
                    std::cerr << "cannot find compact_path_node in compact_graph" << std::endl;
                    return result;
                }
            }
        }
        // std::cout << "transfer original nodes into compact nodes among source destination and paths time:" << wtime() - start_time << std::endl;

        // start_time = wtime();
        OptYen<DeltaSteppingStatic, GraphType, COMPACT_KSP_NUM_THREADS, COMPACT_KSP_PD> optyen(compact_g, k, compact_reverse_sssp_t, paths);

#ifdef TECHNIQUE_BOUND
        // float number could cause one less path in some cases, in order to make sure correct totally, now comment it out or solution 2, make the k bound is greater
        if (actual_k == k) {
            optyen.setKBound(k_bound * 1.01);
        }
#endif

        result = optyen.compute(compact_source, compact_destination, false);
        // end_time = wtime();
        // optyen.print_compact_paths();
        // std::cout << "OptYen compute time:" << end_time - start_time << std::endl;

        compact_g.get_csr_graph()->destroy();

        return result;
    }

    void normalSSSP(const NODE_ID source, const NODE_ID destination) {
        for (int i = 0; i < reverse_unreachable_nodes.size(); ++i) {
            orig.remove_node(reverse_unreachable_nodes[i]);
        }

        DS normal_ds(orig, num_threads);
        normal_ds.template compute_sssp(source);
        normal_sssp_t = normal_ds.get_sssp_tree();
    }

    void reverseSSSP(const NODE_ID source, const NODE_ID destination) {
        DS reverse_ds(*reverse_graph, num_threads);
        reverse_ds.template compute_sssp(destination);
        paths.push_back({reverse_ds.get_distance(source), 0, 0, reverse_ds.get_shortest_path_to(source, true)});
        reverse_unreachable_nodes = reverse_ds.get_unreachable_nodes();
        reverse_sssp_t = reverse_ds.get_sssp_tree();
    }

    void identifyLooplessKBound(const NODE_ID source, const NODE_ID destination) {
        const NODE_ID num_nodes = orig.get_num_nodes();
        nodeWeight* node_weight = new nodeWeight[num_nodes];

        // sort distance
        for (NODE_ID i = 0; i < num_nodes; i++) {
            w_type weight = normal_sssp_t.at(i).tent + reverse_sssp_t.at(i).tent;
            node_weight[i].id = i;
            node_weight[i].distance = weight;
        }
        std::sort(node_weight, node_weight + num_nodes, cmp);

        // calc k nodes
        NODE_ID i;
        for (i = 0; i < num_nodes; i++) {
            w_type path_weight = node_weight[i].distance;
            if (path_weight == INFINITE_DISTANCE) {
                break;
            }

            // check if the shortest path is loop
            NODE_ID node = node_weight[i].id;
            if (is_k_bound_nodes[node]) {
                continue;
            }

            bool is_loop = false;
            NODE_ID reverse_node = node;
            std::vector<NODE_ID> reverse_path;
            std::set<NODE_ID> normal_path = normal_sssp_t.get_shortest_path_set(node);
            while (reverse_sssp_t.at(reverse_node).parent != reverse_node) {
                reverse_node = reverse_sssp_t.at(reverse_node).parent;
                reverse_path.push_back(reverse_node);
                if (normal_path.find(reverse_node) != normal_path.end()) {
                    is_loop = true;
                    break;
                }
            }

            k_bound = path_weight;
            is_k_bound_nodes[node] = true;
            k_bound_nodes_num++;

            if (!is_loop) {
                // add simple path
                for (auto path_node : normal_path) {
                    if (!is_k_bound_nodes[path_node]) {
                        is_k_bound_nodes[path_node] = true;
                        k_bound_nodes_num++;
                    }
                }
                for (auto path_node : reverse_path) {
                    if (!is_k_bound_nodes[path_node]) {
                        is_k_bound_nodes[path_node] = true;
                        k_bound_nodes_num++;
                    }
                }

                // only calc k bound when path is loopless
                actual_k++;

                if (actual_k == k) {
                    // add extra same distance nodes
                    while (++i < num_nodes) {
                        NODE_ID new_node = node_weight[i].id;
                        w_type new_path_weight = node_weight[i].distance;
                        if (FLOATEQUAL(new_path_weight, k_bound)) {
                            if (!is_k_bound_nodes[new_node]) {
                                is_k_bound_nodes[new_node] = true;
                                k_bound_nodes_num++;
                            }
                        } else {
                            break;
                        }
                    }
                    // get k difference loopless paths nodes
                    break;
                }
            }
        }

        if (k_bound_nodes_num <= adaptive_ratio * num_nodes) {
            // remap original node into compact node
            for (NODE_ID node = 0; node < num_nodes; node++) {
                if (is_k_bound_nodes[node]) {
                    k_bound_nodes_map[node] = k_bound_nodes.size();
                    k_bound_nodes.push_back(node);
                }
            }
        }

        delete[] node_weight;
    }

    void regenerateGraphSsspTree(SsspTree& compact_reverse_sssp_t) {
        NODE_ID num_nodes = k_bound_nodes.size();
        for (NODE_ID compact_node = 0; compact_node < num_nodes; ++compact_node) {
            NODE_ID original_node = k_bound_nodes[compact_node];
            w_type tent = reverse_sssp_t[original_node].tent;
            NODE_ID original_parent = reverse_sssp_t[original_node].parent;
            std::map<NODE_ID, NODE_ID>::iterator it = k_bound_nodes_map.find(original_parent);
            if (it != k_bound_nodes_map.end()) {
                NODE_ID compact_parent = it->second;
                compact_reverse_sssp_t[compact_node].tent = tent;
                compact_reverse_sssp_t[compact_node].parent = compact_parent;
            }
        }
    }

    GraphType regenerateGraph() {
        NODE_ID num_nodes = k_bound_nodes.size();
        EDGE_ID num_edges = 0;
        w_type total_weight = 0;
        std::vector<std::vector<std::pair<NODE_ID, w_type>>> edges_per_node;
        edges_per_node.resize(num_nodes);

        for (NODE_ID compact_node = 0; compact_node < num_nodes; ++compact_node) {
            NODE_ID original_node_from = k_bound_nodes[compact_node];
            for (EDGE_ID neighor = orig.begin(original_node_from); neighor < orig.end(original_node_from); ++neighor) {
                const auto csr = orig.get_original_graph().get_csr_graph();
                w_type weight = csr->value[neighor];
                if (actual_k == k && weight > k_bound) {
                    continue;
                }
                NODE_ID original_node_to = csr->adj[neighor];
                std::map<NODE_ID, NODE_ID>::iterator it = k_bound_nodes_map.find(original_node_to);
                if (it != k_bound_nodes_map.end()) {
                    edges_per_node[compact_node].push_back(std::make_pair(it->second, weight));
                    num_edges++;
                    total_weight += weight;
                }
            }
        }

        w_type heaviest_weight = 0;
        CSRGraph* csr_compact = new CSRGraph(num_nodes, num_edges);
        EDGE_ID edge_count = 0;
        for (NODE_ID node_from = 0; node_from < num_nodes; node_from++) {
            csr_compact->begin[node_from] = edge_count;
            for (int j = 0; j < edges_per_node[node_from].size(); ++j) {
                NODE_ID node_to = edges_per_node[node_from][j].first;
                w_type weight = edges_per_node[node_from][j].second;
                csr_compact->adj[edge_count] = node_to;
                csr_compact->value[edge_count] = weight;
                edge_count++;

                if (weight > heaviest_weight) {
                    heaviest_weight = weight;
                }
            }
        }
        csr_compact->begin[num_nodes] = num_edges;
        csr_compact->k_bound_nodes = std::move(k_bound_nodes);

        w_type avg_degree = static_cast<w_type>(num_edges) / static_cast<w_type>(num_nodes);
        w_type avg_weight = total_weight / num_edges;
        w_type dist = paths[0].length;
        int l = paths[0].p.size() - 1;
        w_type delta = (dist / l + avg_weight) / avg_degree;
        // std::cout << "compact graph delta:" << delta << std::endl;
        GraphRW::partition_edges_by_weight(csr_compact, delta);
        // std::cout << "k_bound:" << k_bound << ",k_bound_nodes:" << num_nodes << ",k_bound_edges:" << num_edges << ",all_nodes:" << orig.get_num_nodes() << ",all_edges:" << orig.get_num_edges() << std::endl;
        return GraphType(csr_compact, heaviest_weight, delta);
    }
};

#endif  // _PEEK_ADAPTIVE_WITH_EDGE_SWAP_H