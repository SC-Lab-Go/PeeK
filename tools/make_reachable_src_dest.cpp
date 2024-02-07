#include <DeltaSteppingStatic.h>
#include <GraphRW.h>
#include <OptYen.h>
#include <Yen.h>
#include <math.h>
#include <sys/time.h>

#include <exception>
#include <iostream>
#include <set>

#define KSP_NUM_THREADS_1 1
#define KSP_NUM_THREADS_2 2
#define KSP_NUM_THREADS_4 4
#define KSP_NUM_THREADS_8 8
#define KSP_NUM_THREADS_16 16
#define KSP_NUM_THREADS_32 32
#define KSP_PD true
#define KSP_PD_L1 false

template <typename GraphType>
std::set<std::pair<NODE_ID, NODE_ID>> getSrcAndDest(const GraphType& G, const unsigned int path_node_num_min = 3, const unsigned int path_node_num_max = 8, const unsigned int num_pairs = 10) {
    KSPGraph<GraphType> K(G);
    srand(static_cast<unsigned int>(time(nullptr)));
    std::set<std::pair<NODE_ID, NODE_ID>> src_dest;
    unsigned int avg = num_pairs / (path_node_num_max - path_node_num_min + 1);
    unsigned int remain = num_pairs % (path_node_num_max - path_node_num_min + 1);
    std::vector<unsigned int> path_node_num(path_node_num_max + 1, avg);
    for (unsigned int i = 0; i < remain; i++) {
        path_node_num[i + path_node_num_min]++;
    }

    for (unsigned int i = path_node_num_min; i <= path_node_num_max; ++i) {
        while (path_node_num[i] > 0) {
            NODE_ID src = rand() % G.get_num_nodes();

            DeltaSteppingStatic<GraphType> dss(K, KSP_NUM_THREADS_32);
            dss.template compute<false>(src);
            std::vector<std::pair<NODE_ID, unsigned int>> reachable_nodes = dss.get_reachable_nodes(i, i);
            if (0 == reachable_nodes.size()) {
                continue;
            }
            std::pair<NODE_ID, unsigned int> elem = reachable_nodes[rand() % reachable_nodes.size()];
            NODE_ID dest = elem.first;
            unsigned int node_num = elem.second;

            if (node_num != i || src == dest || src_dest.find(std::make_pair(src, dest)) != src_dest.end()) {
                continue;
            }

            path_node_num[i]--;
            src_dest.insert(std::make_pair(src, dest));

            std::vector<NODE_ID> path = dss.get_shortest_path_to(dest);
            std::cout << "#src:" << src << ",dest:" << dest << ",path_node_num:" << node_num << ",length:" << dss.get_distance(dest) << ",path is:";
            for (auto elem : path) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
    }
    return src_dest;
}


int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Arguments:" << std::endl << " [1] fw_beg_file" << std::endl << " [2] fw_csr_file" << std::endl << " [3] fw_value_file" << std::endl;
        exit(0);
    }

    char default_delta[] = "0.1";
    using GraphType = BasicGraph<true, true>;
    const auto G = GraphRW::read_graph<true, true, true>(argv[1], argv[2], argv[3], default_delta);
    std::set<std::pair<NODE_ID, NODE_ID>> src_dest = getSrcAndDest(G);
    for (auto elem : src_dest) {
        NODE_ID src = elem.first;
        NODE_ID dest = elem.second;
        std::cout << src << " " << dest << std::endl;
    }
}