#ifndef _NODE_H
#define _NODE_H

#include <ostream>
#include <limits>
#include <vector>
#include <set>

// Node types
using NODE_ID = unsigned int;
constexpr NODE_ID NULL_NODE = std::numeric_limits<NODE_ID>::max();

// Edge types
using EDGE_ID = size_t;
constexpr EDGE_ID NULL_EDGE = std::numeric_limits<EDGE_ID>::max();

using w_type = double;      // Can be double or float
static_assert(std::numeric_limits<w_type>::has_infinity, "w_type needs to support numeric_limits<w_type>::infinity()");

using path_type = std::vector<NODE_ID>;

// distances
constexpr w_type INFINITE_DISTANCE = std::numeric_limits<w_type>::infinity();


struct CSRGraph
{
    EDGE_ID *begin          = NULL;
    NODE_ID *adj            = NULL;
    w_type  *value          = NULL;
    NODE_ID *num_light_edges    = NULL;
    NODE_ID *valid_light_num    = NULL;
    NODE_ID *valid_heavy_num    = NULL;

#ifdef TECHNIQUE_BRIDGE_WCC
    NODE_ID *wcc    = NULL;
    std::set<std::pair<NODE_ID, NODE_ID>> *bridge   = NULL;
#endif

    std::vector<NODE_ID> k_bound_nodes;

    NODE_ID num_nodes   = 0;
    EDGE_ID num_edges   = 0;

    CSRGraph() = default;
    ~CSRGraph()
    {
        destroy();
    }


    void destroy()
    {
        if ( begin )
        {
            delete[] begin;
            begin = NULL;
        }
        if ( adj )
        {
            delete[] adj;
            adj = NULL;
        }
        if ( value )
        {
            delete[] value;
            value = NULL;
        }
        if ( num_light_edges )
        {
            delete[] num_light_edges;
            num_light_edges = NULL;
        }
        if ( valid_light_num )
        {
            delete[] valid_light_num;
            valid_light_num = NULL;
        }
        if ( valid_heavy_num )
        {
            delete[] valid_heavy_num;
            valid_heavy_num = NULL;
        }
#ifdef TECHNIQUE_BRIDGE_WCC
        if ( wcc )
        {
            delete[] wcc;
            wcc = NULL;
        }
        if ( bridge )
        {
            delete bridge;
            bridge = NULL;
        }
#endif
        num_nodes   = 0;
        num_edges   = 0;
    }

    explicit CSRGraph( NODE_ID num_nodes, NODE_ID num_edges ) noexcept :
        num_nodes( num_nodes ),
        num_edges( num_edges )
    {
        begin       = new EDGE_ID[num_nodes + 1];
        adj     = new NODE_ID[num_edges];
        value       = new w_type[num_edges];
        num_light_edges = new NODE_ID[num_nodes + 1];
        for ( NODE_ID i = 0; i < num_nodes + 1; i++ )
        {
            num_light_edges[i] = 0;
        }
    }

    void initEdgeSwap()
    {
        if ( valid_light_num )
        {
            delete[] valid_light_num;
            valid_light_num = NULL;
        }
        if ( valid_heavy_num )
        {
            delete[] valid_heavy_num;
            valid_heavy_num = NULL;
        }
        valid_light_num = new NODE_ID[num_nodes];
        valid_heavy_num = new NODE_ID[num_nodes];
    }

    friend std::ostream &operator<<( std::ostream & output, const CSRGraph & csr ) noexcept
    {
        output << "csr_graph:[num_nodes = " << csr.num_nodes << "; num_edges = " << csr.num_edges << "]" << " [begin:" << csr.begin << " adj:" << csr.adj << " value:" << csr.value << "]" << std::endl;
        if ( csr.num_nodes < 20 && csr.num_edges < 100 )
        {
            for ( NODE_ID node_from = 0; node_from < csr.num_nodes; node_from++ )
            {
                for ( EDGE_ID neighor = csr.begin[node_from]; neighor < csr.begin[node_from + 1]; neighor++ )
                {
                    w_type  weight  = csr.value[neighor];
                    NODE_ID node_to = csr.adj[neighor];
                    output << node_from << "->" << node_to << " " << weight << " " << csr.num_light_edges[node_from] << std::endl;
                }
            }
        }
        return(output);
    }
};



#endif  // _NODE_H
