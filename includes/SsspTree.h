#ifndef _SSSP_TREE_H
#define _SSSP_TREE_H

#include <KSPGraph.h>
#include <algorithm>
#include <cassert>

#ifdef DO_STATISTICS
#include <exception>

class NoPathException : public std::exception
{
    const char* what() const noexcept override
    {
        return "Source and target are not connected";
    }
};

#endif

class SsspTree
{
public:
    enum class path_direction { s2t, t2s };   // todo is this really needed?

private:
    struct tree_node;   // is defined at the end

    const unsigned int num_nodes;
    std::vector<tree_node> the_tree;
    path_direction path_type;

    NODE_ID destination = NULL_NODE;    // this should only be set if the SSSP tree is only used for SPSP computation.

public:
    /**
     * @param num_nodes The number of nodes in the tree
     */
    explicit SsspTree(const unsigned int num_nodes, const path_direction pt = path_direction::s2t) noexcept
        : num_nodes(num_nodes), the_tree(num_nodes), path_type(pt)
    {}

    std::vector<std::pair<NODE_ID, unsigned int>> get_reachable_nodes(unsigned int path_node_num_min = 3, unsigned int path_node_num_max = 8) {
        std::vector<std::pair<NODE_ID, unsigned int>> reachable_nodes;
        for (NODE_ID node = 0; node < the_tree.size(); node++) {
            if (INFINITE_DISTANCE != the_tree[node].tent) {
                std::vector<NODE_ID> path = get_shortest_path(node);
                unsigned int node_num = path.size();
                if (path_node_num_min <= node_num && node_num <= path_node_num_max) {
                    reachable_nodes.push_back(std::make_pair(node, node_num));
                }
            }
        }

        return reachable_nodes;
    }

    /**
     * Set a destination if this SSSP tree is used for SPSP computation.
     * @param destination
     */
    void set_destination(const NODE_ID destination) noexcept
    {
        this->destination = destination;
    }

    tree_node& operator[](const unsigned int i) noexcept
    {
        assert(i < the_tree.size());
        return the_tree[i];
    }

    const tree_node& at(const unsigned int i) const noexcept
    {
        assert(i < the_tree.size());
        return the_tree[i];
    }

    SsspTree& operator=(const SsspTree& to_copy) noexcept
    {
        the_tree = to_copy.the_tree;
        return *this;
    }

    void reset() noexcept
    {
        destination = NULL_NODE;

        for(unsigned int i = 0; i < num_nodes; i++)
        {
            the_tree[i].reset();
        }
    }

    void print_path(const NODE_ID src, const NODE_ID dest) const noexcept
    {
        std::vector<NODE_ID> path;
        path.push_back(dest);

        NODE_ID p1 = dest;

        while(p1 != src)
        {
            const NODE_ID p2 = the_tree[p1].parent;
            if(p2 == NULL_NODE)
            {
                std::cout << "no shortest path in this component!";
                return;
            }

            path.push_back(p2);
            p1 = p2;
        }

        std::cout << "distance = " << the_tree[dest].tent << std::endl;
        std::cout << "shortest path: ";
        for(size_t i = path.size(); i > 0; i--)
        {
            std::cout << path[i-1] << " ";
        }

        std::cout << std::endl;
    }

    void print_tree() const
    {
        for(unsigned int i = 0; i < the_tree.size(); i++)
        {
            std::cout << "node " << i << ", distance " <<  the_tree[i].tent << ", parent " << the_tree[i].parent << std::endl;
        }
    }

    std::vector<NODE_ID> get_shortest_path(const NODE_ID dest, const bool get_reverse_path = false) const
    {
        assert(destination == NULL_NODE || destination == dest);    // avoid to ask for shortest paths on incomplete trees

        if(get_distance(dest) == INFINITE_DISTANCE)
        {
#ifdef DO_STATISTICS
            throw NoPathException();        // this is not thread safe!
#else
            std::cerr << "ERROR: no parent for node " << dest << std::endl;
            return std::vector<NODE_ID>();
#endif
        }

        std::vector<NODE_ID> path;
        path.push_back(dest);
        NODE_ID node = dest;

        while(the_tree[node].parent != node)
        {
            node = the_tree[node].parent;
            path.push_back(node);
        }

        if(!get_reverse_path)
            std::reverse(path.begin(), path.end());

        return path;
    }
    
    std::set<NODE_ID> get_shortest_path_set(const NODE_ID dest, const bool include_dest = false) const {
        if (get_distance(dest) == INFINITE_DISTANCE) {
            return std::set<NODE_ID>();
        }

        std::set<NODE_ID> path;
        if (include_dest) {
            path.insert(dest);
        }
        NODE_ID node = dest;

        while (the_tree[node].parent != node) {
            node = the_tree[node].parent;
            path.insert(node);
        }

        return path;
    }

    /**
     * Returns the distance from the source node to the specified target.
     * @param target
     * @return The distance from the source node to the specified target.
     */
    w_type get_distance(const NODE_ID target) const noexcept
    {
        assert(destination == NULL_NODE ||
               destination == target);    // avoid to ask for shortest paths on incomplete trees

        return the_tree[target].tent;
    }

    std::vector<NODE_ID> get_unreachable_nodes(){
        std::vector<NODE_ID> unreachable_nodes;
        for(unsigned int i = 0; i < the_tree.size(); i++)
        {
            if (INFINITE_DISTANCE==the_tree[i].tent)
            {
                unreachable_nodes.push_back(i);
            }
        }

        return unreachable_nodes;
    }


    void get_all_distance(w_type *distance_from_src_2_nodes){
        for(unsigned int i = 0; i < the_tree.size(); i++)
        {
            distance_from_src_2_nodes[i]=the_tree[i].tent;
        }
    }


private:
    struct tree_node
    {
        w_type tent;        // distance from source to this node
        NODE_ID parent;     // parent node ID

        tree_node() noexcept
            : tent(INFINITE_DISTANCE),
              parent(NULL_NODE)
        {}

        void reset() noexcept
        {
            tent = INFINITE_DISTANCE;
            parent = NULL_NODE;
        }

        friend std::ostream& operator<<(std::ostream& o, const tree_node& t)
        {
            o << "tent:" << t.tent << " parent:" << t.parent;
            return o;
        }
    };

};

#endif  // _SSSP_TREE_H
