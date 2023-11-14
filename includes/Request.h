#ifndef _REQUEST_H
#define _REQUEST_H

#include <KSPGraph.h>

/**
 * Stores a request to improve the distance of node to the source by the edge (parent_node, node) to the distance proposed_tent.
 */
struct UpdateRequest
{
    NODE_ID parent_node;
    NODE_ID node;
    w_type proposed_tent;

    UpdateRequest() noexcept
        : parent_node(NULL_NODE), node(NULL_NODE),
          proposed_tent(INFINITE_DISTANCE)
    {}

    UpdateRequest(const NODE_ID v, const NODE_ID node, const w_type proposed_tent) noexcept
        : parent_node(v), node(node),
          proposed_tent(proposed_tent)
    {}

    friend std::ostream &operator<<(std::ostream& o, const UpdateRequest& r)
    {
        o << "[" << r.parent_node << "," << r.node << "," << r.proposed_tent << "]";
        return o;
    }
};

#endif  // _REQUEST_H