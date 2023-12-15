//
// Created by Alex Schickedanz <alex@ae.cs.uni-frankfurt.de> on 22.05.17.
//

#ifndef _GRAPH_RW_H
#define _GRAPH_RW_H

#include "BasicGraph.h"
#include "SsspTree.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <sstream>
#include <random>
#include <chrono>
#include <tuple>
#include <memory>
#include <sys/stat.h>

#include <parallel/algorithm>

/**
 * Reads and writes graphs from and to files.
 */
class GraphRW
{
public:
    template<bool directed = false, bool weighted = true, bool add_reverse_graph = false>
    static BasicGraph<directed, weighted> read_graph( const char *beg_file, const char *adj_file, const char *value_file,
#ifdef TECHNIQUE_BRIDGE_WCC
        const char *wcc_file, const char *bridge_file, 
#endif
        const char* delta_option )
    {
        w_type total_weight = 0;
        w_type  heaviest_weight = weighted ? 0.0 : 1.0;
        NODE_ID num_nodes   = fsize( beg_file ) / sizeof(EDGE_ID) - 1;
        EDGE_ID num_edges   = fsize( adj_file ) / sizeof(NODE_ID);

        FILE *file = fopen( beg_file, "rb" );
        if ( file == NULL )
        {
            std::cout << beg_file << " cannot open " << beg_file << std::endl;
            exit( -1 );
        }
        EDGE_ID *tmp_begin  = new EDGE_ID[num_nodes + 1];
        int ret     = fread( tmp_begin, sizeof(EDGE_ID), num_nodes + 1, file );
        assert( ret == num_nodes + 1 );
        fclose( file );

        file = fopen( adj_file, "rb" );
        if ( file == NULL )
        {
            std::cout << adj_file << " cannot open " << adj_file << std::endl;
            exit( -1 );
        }
        NODE_ID *tmp_adj = new NODE_ID[num_edges];
        ret = fread( tmp_adj, sizeof(NODE_ID), num_edges, file );
        assert( ret == num_edges );
        fclose( file );

        file = fopen( value_file, "rb" );
        if ( file == NULL )
        {
            std::cout << value_file << " cannot open " << value_file << std::endl;
            exit( -1 );
        }
        w_type *tmp_value = new w_type[num_edges];
        ret = fread( tmp_value, sizeof(w_type), num_edges, file );
        assert( ret == num_edges );
        fclose( file );

        CSRGraph* csr = new CSRGraph( num_nodes, num_edges );

        for ( NODE_ID i = 0; i < num_nodes + 1; ++i )
        {
            csr->begin[i] = (EDGE_ID) tmp_begin[i];
        }

        for ( EDGE_ID i = 0; i < num_edges; ++i )
        {
            csr->adj[i] = (NODE_ID) tmp_adj[i];
        }

        for ( EDGE_ID i = 0; i < num_edges; ++i )
        {
            csr->value[i] = (w_type) tmp_value[i];
            if ( weighted )
            {
                total_weight += csr->value[i];
                if ( csr->value[i] > heaviest_weight )
                    heaviest_weight = csr->value[i];
            }
        }


        delete[] tmp_begin;
        delete[] tmp_adj;
        delete[] tmp_value;

#ifdef TECHNIQUE_BRIDGE_WCC
        //read wcc file
        csr->wcc = new NODE_ID[num_nodes];
        std::ifstream wcc(wcc_file);
        if(!wcc.is_open())
        {
            throw std::invalid_argument("Could not open the wcc file.");
        }
        std::string wcc_line;
        NODE_ID node_count=0;
        while(std::getline(wcc, wcc_line))
        {
            // Skip comment lines
            if(wcc_line.empty() || wcc_line.at(0) == '#')
                continue;

            NODE_ID wcc_id;
            std::stringstream ss(wcc_line);
            ss >> wcc_id ;
            csr->wcc[node_count] = wcc_id;
            node_count++;
        }

        //read bridge file
        csr->bridge = new std::set<std::pair<NODE_ID, NODE_ID>>;
        std::ifstream bridge(bridge_file);
        if(!bridge.is_open())
        {
            throw std::invalid_argument("Could not open the bridge file.");
        }
        std::string bridge_line;
        while(std::getline(bridge, bridge_line))
        {
            // Skip comment lines
            if(bridge_line.empty() || bridge_line.at(0) == '#')
                continue;

            NODE_ID src,dest;
            std::stringstream ss(bridge_line);
            ss >> src >> dest;
            csr->bridge->insert(std::make_pair(src,dest));
        }
#endif
        

        w_type avg_degree = static_cast<w_type>(num_edges) / static_cast<w_type>(num_nodes);
        w_type avg_weight = total_weight / num_edges;
        w_type delta = avg_weight / avg_degree;
        // std::cout << "original graph delta:" << delta << std::endl;
        GraphRW::partition_edges_by_weight( csr, delta );

        
        if(directed && add_reverse_graph)
        {
        auto g_reverse = GraphRW::get_reverse_graph<directed, weighted>(csr, heaviest_weight, delta);
        return BasicGraph<directed, weighted>(csr, heaviest_weight, delta, std::move(g_reverse));
        }

        return(BasicGraph<directed, weighted>( csr, heaviest_weight, delta ) );
    }


public:
    template<bool directed, bool weighted>
    static std::unique_ptr<const BasicGraph<directed, weighted> > get_reverse_graph( CSRGraph* csr_orig, const w_type heaviest_weight, const w_type delta )
    {
        if ( !directed )
            throw std::invalid_argument( "Don't compute the reverse of a undirected graph." );

        CSRGraph * csr = new CSRGraph( csr_orig->num_nodes, csr_orig->num_edges );
        std::vector<std::vector<std::pair<NODE_ID, w_type> > >  edges_per_node( csr_orig->num_nodes );

        for ( NODE_ID node_from = 0; node_from < csr_orig->num_nodes; node_from++ )
        {
            for ( EDGE_ID neighor = csr_orig->begin[node_from]; neighor < csr_orig->begin[node_from + 1]; neighor++ )
            {
                w_type  weight  = csr_orig->value[neighor];
                NODE_ID node_to = csr_orig->adj[neighor];
                edges_per_node[node_to].emplace_back( std::make_pair( node_from, weight ) );
            }
        }

        NODE_ID total_edge_count = 0;
        for ( NODE_ID node_to = 0; node_to < csr->num_nodes; node_to++ )
        {
            csr->begin[node_to] = total_edge_count;
            for ( NODE_ID edge_count = 0; edge_count < edges_per_node[node_to].size(); edge_count++ )
            {
                NODE_ID node_from   = edges_per_node[node_to][edge_count].first;
                w_type  weight      = edges_per_node[node_to][edge_count].second;
                csr->adj[total_edge_count]  = node_from;
                csr->value[total_edge_count]    = weight;
                total_edge_count++;
            }
        }
        csr->begin[csr->num_nodes] = csr->num_edges; /* setting the new sentinel */

        assert( total_edge_count == csr->num_edges );
        GraphRW::partition_edges_by_weight( csr, delta );

        /* the object gets deleted in the destructor of the original graph. */
        return(std::make_unique<const BasicGraph<directed, weighted> >( csr, heaviest_weight, delta ) );
    }

public:
    /**
     * Sets delta according to delta_option. It does not precompute any edge weights.
     * If delta_option is not valid, delta is set to 0.1
     * @param delta_option Takes one of the following
     *                     "-D" for Dijkstra, which sets delta to 1
     *                     "-BF" for Bellman-Ford which sets delta to the heaviest weight
     *                     "<float>" sets delta to this number, if it is >0
     */
    static w_type get_delta(const char* delta_option, w_type heaviest_weight)
    {
        if(delta_option[0] == '-')
        {
            const char dijkstra[] = "-D";
            const char bellman_ford[] = "-BF";

            if(strcmp(delta_option, dijkstra) == 0)
            {
                return 1;  // todo is this correct if edge weights < 1 exist?   this needs to be the lightest weight...
            }

            if(strcmp(delta_option, bellman_ford) == 0)
            {
                return heaviest_weight;
            }

            throw std::invalid_argument("wrong option for delta. Only -D (Dijkstra) and -BF (Bellman Ford) is accepted.");
        }
        else
        {
            w_type delta = strtod(delta_option, nullptr);

            if(delta <= 0)
            {
                std::cerr << "WARNING: delta must be > 0, delta will be set to 0.1" << std::endl;
                return 0.1;
            }

            return delta;
        }
    }


    static off_t fsize( const char *filename )
    {
        struct stat st;
        if ( stat( filename, &st ) == 0 )
            return(st.st_size);
        return(-1);
    }


public:
    template<bool directed = false, bool weighted = true>
    static void partition_edges_by_weight( CSRGraph* csr, const w_type delta )
    {
        if ( !weighted )
        {
            /*
             * Every edge has the same weight. If the first edge is light, every edge is. Then each number of light
             * edges needs to be set to the number of neighbors.
             */
            if ( csr->value && csr->value[0] <= delta )
            {
                for ( NODE_ID node_from = 0; node_from < csr->num_nodes; node_from++ )
                {
                    csr->num_light_edges[node_from] = static_cast<NODE_ID>(csr->begin[node_from + 1] - csr->begin[node_from]);
                }
            }
            /* else, every node has only heavy out edges or zero light edges. Zero is the default value so no changes needed. */

            return;
        }

        for ( NODE_ID node_from = 0; node_from < csr->num_nodes; node_from++ )
        {
            std::vector<std::pair<NODE_ID, w_type> >    light_weight;
            std::vector<std::pair<NODE_ID, w_type> >    heavy_weight;
            for ( EDGE_ID neighor = csr->begin[node_from]; neighor < csr->begin[node_from + 1]; neighor++ )
            {
                w_type  weight  = csr->value[neighor];
                NODE_ID node_to = csr->adj[neighor];
                if ( weight <= delta )
                {
                    light_weight.push_back( std::make_pair( node_to, weight ) );
                }else{
                    heavy_weight.push_back( std::make_pair( node_to, weight ) );
                }
            }

            csr->num_light_edges[node_from] = light_weight.size();
            if ( csr->num_light_edges[node_from] > 0 )
            {
                for ( int i = 0; i < light_weight.size(); i++ )
                {
                    csr->value[csr->begin[node_from] + i]   = light_weight[i].second;
                    csr->adj[csr->begin[node_from] + i] = light_weight[i].first;
                }
                for ( int i = 0; i < heavy_weight.size(); i++ )
                {
                    csr->value[csr->begin[node_from] + light_weight.size() + i] = heavy_weight[i].second;
                    csr->adj[csr->begin[node_from] + light_weight.size() + i]   = heavy_weight[i].first;
                }
            }
        }
    }
};

#endif  // _GRAPH_RW_H
