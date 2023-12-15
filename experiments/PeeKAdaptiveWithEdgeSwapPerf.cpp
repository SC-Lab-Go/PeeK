#include <DeltaSteppingStatic.h>
#include <GraphRW.h>
#include <KSPGraph.h>
#include <MaxHeapNoDuplicate.h>
#include <OptYen.h>
#include <PeeKAdaptiveWithEdgeSwap.h>
#include <sys/time.h>

#include <numeric>

#ifdef KSP_PARALLEL_DEVIATIONS_L1
#define KSP_PD_L1 true
#else
#define KSP_PD_L1 false
#endif

#ifndef KSP_NUM_THREADS
#define KSP_NUM_THREADS 1
#endif

#if defined(KSP_PARALLEL_DEVIATIONS_L2) || defined(KSP_PARALLEL_DEVIATIONS_L1)
#define KSP_PD true
#else
#define KSP_PD false
#endif

#define KSP_NUM_THREADS_32 32

std::set<std::pair<NODE_ID, NODE_ID>> getPairsOfSrcAndDestFromFile(const char* path = NULL) {
    std::set<std::pair<NODE_ID, NODE_ID>> src_dest;
    src_dest.clear();

    if (path) {
        std::ifstream graph_file(path);
        if (!graph_file.is_open()) throw std::invalid_argument("Could not open the file.");

        std::string line;
        while (std::getline(graph_file, line)) {
            // Skip comment lines
            if (line.empty() || line.at(0) == '#') continue;

            NODE_ID src, dest;
            std::stringstream ss(line);
            ss >> src >> dest;

            src_dest.insert(std::make_pair(src, dest));
        }
    }

    return src_dest;
}

// restore graph because edge swap can change original graph
template <bool directed = true, bool weighted = true>
void restoreGraph(const BasicGraph<directed, weighted>& G, const BasicGraph<directed, weighted>& G_backup) {
    const auto csr_backup = G_backup.get_csr_graph();
    const auto csr = G.get_csr_graph();

    csr->num_nodes = csr_backup->num_nodes;
    csr->num_edges = csr_backup->num_edges;
    csr->k_bound_nodes.clear();
#pragma omp parallel for num_threads(KSP_NUM_THREADS_32)
    for (NODE_ID i = 0; i < csr_backup->num_nodes + 1; i++) {
        csr->begin[i] = csr_backup->begin[i];
        csr->num_light_edges[i] = csr_backup->num_light_edges[i];
    }
#pragma omp parallel for num_threads(KSP_NUM_THREADS_32)
    for (NODE_ID i = 0; i < csr_backup->num_edges; i++) {
        csr->adj[i] = csr_backup->adj[i];
        csr->value[i] = csr_backup->value[i];
    }
}

int main(int argc, char** argv) {
    double start_time, end_time;

    if (argc < 5) {
        std::cout << "Arguments:" << std::endl << " [1] fw_beg_file" << std::endl << " [2] fw_adj_file" << std::endl << " [3] fw_value_file" << std::endl << " [4] src_dest_file" << std::endl;
        exit(0);
    }

    char default_delta[] = "0.1";
    char* src_dest_file = argv[4];

    std::vector<unsigned int> k = {8, 128};

    const unsigned int runs = 1;

    using GraphType = BasicGraph<true, true>;
    start_time = wtime();
    const auto G = GraphRW::read_graph<true, true, true>(argv[1], argv[2], argv[3], default_delta);
    const auto G_backup = GraphRW::read_graph<true, true, true>(argv[1], argv[2], argv[3], default_delta);
    end_time = wtime();
    std::cout << "file read time:" << end_time - start_time << std::endl;
    std::cout << "KSP_NUM_THREADS:" << KSP_NUM_THREADS << ",KSP_PD:" << KSP_PD << ",KSP_PD_L1:" << KSP_PD_L1 << ",COMPACT_KSP_NUM_THREADS:" << COMPACT_KSP_NUM_THREADS << ",COMPACT_KSP_PD:" << COMPACT_KSP_PD << std::endl;

    std::set<std::pair<NODE_ID, NODE_ID>> src_dest = getPairsOfSrcAndDestFromFile(src_dest_file);

    for (auto _k : k) {
        std::vector<double> avg_run_times;
        avg_run_times.clear();
        double avg_run_time_max = 0;
        double avg_run_time_min = INFINITE_DISTANCE;
        for (auto elem : src_dest) {
            NODE_ID src = elem.first;
            NODE_ID dest = elem.second;

            std::vector<double> run_times;
            run_times.clear();
            for (unsigned i = 0; i < runs; i++) {
                PeeKAdaptiveWithEdgeSwap<DeltaSteppingStatic, GraphType, KSP_NUM_THREADS, KSP_PD> peek(G, _k);
                start_time = wtime();
                const auto paths = peek.compute(src, dest);
                end_time = wtime();
                double run_time = end_time - start_time;
                run_times.push_back(run_time);
                restoreGraph<true, true>(G, G_backup);
            }

            double total_run_time = std::accumulate(run_times.begin(), run_times.end(), 0.0);
            double avg_run_time = total_run_time / runs;

            avg_run_times.push_back(avg_run_time);

            if (avg_run_time_max < avg_run_time) {
                avg_run_time_max = avg_run_time;
            }
            if (avg_run_time_min > avg_run_time) {
                avg_run_time_min = avg_run_time;
            }

            std::cout << "k:" << _k << ",src:" << src << ",dest:" << dest << ",avg_run_time:" << avg_run_time << std::endl;
        }

        double total_avg_run_time = std::accumulate(avg_run_times.begin(), avg_run_times.end(), 0.0);
        double avg_run_time = total_avg_run_time / avg_run_times.size();

        std::cout << "average k:" << _k << ",avg_run_time:" << avg_run_time << ",avg_run_time_max:" << avg_run_time_max << ",avg_run_time_min:" << avg_run_time_min << std::endl;
    }

    return 0;
}