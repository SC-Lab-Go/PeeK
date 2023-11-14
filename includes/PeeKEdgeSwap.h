#ifndef _PEEK_EDGE_SWAP_H
#define _PEEK_EDGE_SWAP_H

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
class PeeKEdgeSwap : public KSP<GraphType, num_threads, pps> {
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
    SsspTree reverse_sssp_t;
    SsspTree normal_sssp_t;
    unsigned int actual_k;
    bool* is_k_bound_nodes;
    unsigned int k_bound_nodes_num;

   public:
    explicit PeeKEdgeSwap(const GraphType& g, const unsigned int k) noexcept : ParentType(g, k), reverse_graph(GraphType::is_directed ? std::make_unique<KSPGraph<GraphType>>(g.get_reverse_graph()) : nullptr), reverse_sssp_t(g.get_num_nodes()), normal_sssp_t(g.get_num_nodes()) {
        assert(GraphType::is_directed);
        assert(GraphType::is_weighted);
        k_bound = INFINITE_DISTANCE;
        actual_k = 0;
        is_k_bound_nodes = new bool[g.get_num_nodes()]();
        k_bound_nodes_num = 0;
    }

    ~PeeKEdgeSwap() { delete[] is_k_bound_nodes; }

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

        // start_time = wtime();
        OptYenEdgeSwap<DeltaSteppingStatic, GraphType, num_threads, COMPACT_KSP_PD> optyen(orig.get_original_graph(), k, reverse_sssp_t, paths, reverse_graph, is_k_bound_nodes);

#ifdef TECHNIQUE_BOUND
        // float number could cause one less path in some cases, in order to make sure correct totally, now comment it out or solution 2, make the k bound is greater
        if (actual_k == k) {
            optyen.setKBound(k_bound * 1.01);
        }
#endif

        result = optyen.compute(source, destination);
        // end_time = wtime();
        // optyen.print_paths();
        // std::cout << "OptYen compute time:" << end_time - start_time << std::endl;

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

        delete[] node_weight;
    }
};

#endif  // _PEEK_EDGE_SWAP_H
