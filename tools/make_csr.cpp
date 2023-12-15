#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using NODE_ID = unsigned int;
using EDGE_ID = size_t;
using w_type = double;

template <typename T>
void writeBinaryFile(std::vector<T> data, std::string output_file) {
    std::ofstream ofs(output_file, std::ios::out | std::ios::binary);
    ofs.write((const char *)&data[0], data.size() * sizeof(T));
    ofs.close();
}

void writeCSR(std::vector<std::vector<std::pair<NODE_ID, w_type> > > &edges_per_node, NODE_ID num_nodes, std::string begin_file, std::string adj_file, std::string value_file) {
    EDGE_ID total_num_edges = 0;
    for (NODE_ID node_from = 0; node_from < num_nodes; node_from++) {
        total_num_edges += edges_per_node[node_from].size();
    }

    std::vector<EDGE_ID> begin(num_nodes + 1); /* adding a sentinel at the end */
    std::vector<NODE_ID> adj(total_num_edges);
    std::vector<w_type> value(total_num_edges);
    EDGE_ID edge_count = 0;
    for (NODE_ID i = 0; i < num_nodes; i++) {
        begin[i] = edge_count;
        for (int j = 0; j < edges_per_node[i].size(); ++j) {
            adj[edge_count] = edges_per_node[i][j].first;
            value[edge_count] = edges_per_node[i][j].second;
            edge_count++;
            
        }
    }
    // set the sentinel
    begin[num_nodes] = total_num_edges;

    writeBinaryFile<EDGE_ID>(begin, begin_file);
    writeBinaryFile<NODE_ID>(adj, adj_file);
    writeBinaryFile<w_type>(value, value_file);
}

int main(int argc, char const *argv[]) {
    if (argc != 5) {
        std::cout << "Arguments:" << std::endl
                  << " [1] graph file"
                  << " [2] begin.bin"
                  << " [3] adj.bin"
                  << " [4] value.bin" << std::endl;
        exit(0);
    }
    std::string begin(argv[2]);
    std::string adj(argv[3]);
    std::string value(argv[4]);

    auto *graph_file_path = const_cast<char *>(argv[1]);

    std::ifstream graph_file(graph_file_path);

    if (!graph_file.is_open()) throw std::invalid_argument("Could not open the file.");

    NODE_ID num_nodes;
    EDGE_ID num_edges;
    std::string line;
    std::getline(graph_file, line);
    std::stringstream s(line);
    s >> num_nodes;

    std::vector<std::vector<std::pair<NODE_ID, w_type> > > edges_per_node_directed;
    edges_per_node_directed.resize(num_nodes);

    while (std::getline(graph_file, line)) {
        // Skip comment lines
        if (line.empty() || line.at(0) == '#') continue;

        NODE_ID src, dest;
        w_type weight;
        std::stringstream ss(line);
        ss >> src >> dest >> weight;

        edges_per_node_directed[src].push_back(std::make_pair(dest, weight));
    }

    graph_file.close();

    writeCSR(edges_per_node_directed, num_nodes, begin, adj, value);

    return (0);
}
