#include <DeltaSteppingStatic.h>
#include <GraphRW.h>
#include <PeeKAdaptiveWithEdgeSwap.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

using GraphType = BasicGraph<true, true>;
using PeeKType = PeeKAdaptiveWithEdgeSwap<DeltaSteppingStatic, GraphType, 1, false>;

struct SimplePath {
    std::vector<NODE_ID> nodes;
    w_type length = INFINITE_DISTANCE;
};

GraphType build_graph(NODE_ID num_nodes, const std::vector<std::tuple<NODE_ID, NODE_ID, w_type>>& edges) {
    auto csr = new CSRGraph(num_nodes, edges.size());

    std::vector<std::vector<std::pair<NODE_ID, w_type>>> adjacency(num_nodes);
    for (const auto& edge : edges) {
        const auto [u, v, w] = edge;
        adjacency[u].push_back({v, w});
    }

    EDGE_ID cursor = 0;
    w_type heaviest_weight = std::numeric_limits<w_type>::lowest();
    w_type total_weight = 0.0;
    for (NODE_ID node = 0; node < num_nodes; ++node) {
        csr->begin[node] = cursor;
        for (const auto& [neighbor, weight] : adjacency[node]) {
            csr->adj[cursor] = neighbor;
            csr->value[cursor] = weight;
            heaviest_weight = std::max(heaviest_weight, weight);
            total_weight += weight;
            ++cursor;
        }
    }
    csr->begin[num_nodes] = cursor;

    if (heaviest_weight <= 0) {
        heaviest_weight = 1.0;
    }

    const w_type avg_degree = edges.empty() ? 0.0 : static_cast<w_type>(edges.size()) / static_cast<w_type>(num_nodes);
    const w_type avg_weight = edges.empty() ? 0.0 : total_weight / static_cast<w_type>(edges.size());
    const w_type delta = avg_degree > 0.0 ? std::max<w_type>(std::abs(avg_weight) / std::max<w_type>(avg_degree, 1e-9), 1e-6)
                                          : 1.0;

    GraphRW::partition_edges_by_weight(csr, delta);
    auto reverse_graph = GraphRW::get_reverse_graph<true, true>(csr, heaviest_weight, delta);

    return GraphType(csr, heaviest_weight, delta, std::move(reverse_graph));
}

void enumerate_paths_recursive(const GraphType& graph, NODE_ID current, NODE_ID destination, std::vector<bool>& visited,
                               std::vector<NODE_ID>& stack, w_type distance, std::vector<SimplePath>& results) {
    stack.push_back(current);

    if (current == destination) {
        results.push_back({stack, distance});
        stack.pop_back();
        return;
    }

    visited[current] = true;
    const CSRGraph* csr = graph.get_csr_graph();
    for (EDGE_ID e = csr->begin[current]; e < csr->begin[current + 1]; ++e) {
        const NODE_ID next = csr->adj[e];
        if (visited[next]) {
            continue;
        }

        enumerate_paths_recursive(graph, next, destination, visited, stack, distance + csr->value[e], results);
    }

    visited[current] = false;
    stack.pop_back();
}

std::vector<SimplePath> brute_force_k_shortest(const GraphType& graph, NODE_ID source, NODE_ID destination, std::size_t k) {
    std::vector<SimplePath> all_paths;
    std::vector<bool> visited(graph.get_num_nodes(), false);
    std::vector<NODE_ID> stack;
    enumerate_paths_recursive(graph, source, destination, visited, stack, 0.0, all_paths);

    std::sort(all_paths.begin(), all_paths.end(), [](const SimplePath& lhs, const SimplePath& rhs) {
        if (!FLOATEQUAL(lhs.length, rhs.length)) {
            return lhs.length < rhs.length;
        }
        return lhs.nodes < rhs.nodes;
    });

    if (all_paths.size() > k) {
        all_paths.resize(k);
    }

    return all_paths;
}

void print_paths(const std::vector<Path>& paths, std::size_t limit) {
    std::size_t printed = 0;
    for (const auto& path : paths) {
        if (path.length == INFINITE_DISTANCE) {
            continue;
        }
        std::cout << "  length=" << path.length << " path=";
        for (auto node : path.p) {
            std::cout << node << ' ';
        }
        std::cout << "\n";
        if (++printed >= limit) {
            break;
        }
    }
}

void print_simple_paths(const std::vector<SimplePath>& paths) {
    for (const auto& path : paths) {
        std::cout << "  length=" << path.length << " path=";
        for (auto node : path.nodes) {
            std::cout << node << ' ';
        }
        std::cout << "\n";
    }
}

void test_uniform_grid() {
    std::cout << "== Uniform grid test ==\n";
    const int rows = 4;
    const int cols = 4;
    const NODE_ID num_nodes = rows * cols;
    std::vector<std::tuple<NODE_ID, NODE_ID, w_type>> edges;
    edges.reserve(rows * cols * 2);

    auto node_id = [cols](int r, int c) { return static_cast<NODE_ID>(r * cols + c); };

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const NODE_ID u = node_id(r, c);
            if (c + 1 < cols) {
                edges.emplace_back(u, node_id(r, c + 1), 1.0);
            }
            if (r + 1 < rows) {
                edges.emplace_back(u, node_id(r + 1, c), 1.0);
            }
        }
    }

    auto graph = build_graph(num_nodes, edges);
    PeeKType peek(graph, 20);
    const NODE_ID source = node_id(0, 0);
    const NODE_ID destination = node_id(rows - 1, cols - 1);
    const auto peek_paths = peek.compute(source, destination);

    std::cout << "K bound value: " << peek.getKBound() << "\n";
    std::cout << "Nodes kept after pruning: " << peek.getKBoundNodeCount() << " / " << graph.get_num_nodes() << "\n";
    std::cout << "Used edge-swap mode: " << std::boolalpha << peek.usedEdgeSwapMode() << "\n";
    std::cout << "Paths returned (up to 5):\n";
    print_paths(peek_paths, 5);

    const auto brute_paths = brute_force_k_shortest(graph, source, destination, 5);
    std::cout << "Ground-truth top paths (up to 5):\n";
    print_simple_paths(brute_paths);

    graph.get_csr_graph()->destroy();
}

void test_negative_weight() {
    std::cout << "\n== Negative edge test ==\n";
    const NODE_ID num_nodes = 4;
    std::vector<std::tuple<NODE_ID, NODE_ID, w_type>> edges = {
        {0, 1, 2.0}, {1, 3, 2.0}, {0, 2, 3.0}, {2, 3, 3.0}, {1, 2, -4.0}};

    auto graph = build_graph(num_nodes, edges);
    PeeKType peek(graph, 3);
    const NODE_ID source = 0;
    const NODE_ID destination = 3;

    const auto peek_paths = peek.compute(source, destination);
    const auto brute_paths = brute_force_k_shortest(graph, source, destination, 3);

    std::cout << "PeeK paths:\n";
    print_paths(peek_paths, 3);
    std::cout << "Ground-truth paths:\n";
    print_simple_paths(brute_paths);

    for (std::size_t i = 0; i < brute_paths.size() && i < peek_paths.size(); ++i) {
        const auto& expected = brute_paths[i];
        const auto& actual = peek_paths[i];
        const bool length_mismatch = !FLOATEQUAL(expected.length, actual.length);
        const bool sequence_mismatch = std::vector<NODE_ID>(actual.p.begin(), actual.p.end()) != expected.nodes;
        if (length_mismatch || sequence_mismatch) {
            std::cout << "  Mismatch at rank " << i + 1 << ": PeeK length=" << actual.length
                      << ", expected=" << expected.length << '\n';
            break;
        }
    }

    graph.get_csr_graph()->destroy();
}

int main() {
    test_uniform_grid();
    test_negative_weight();
    return 0;
}