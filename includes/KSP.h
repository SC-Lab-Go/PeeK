#ifndef _KSP_H
#define _KSP_H

#include <KSPGraph.h>
#include <Profiling.h>
#include <limits>
#include <deque>
#include <set>
#include <DeltaSteppingStatic.h>

#ifndef NO_SSSP_EARLY_STOPPING
    #define SSSP_EARLY_STOPPING true
#else
    #define SSSP_EARLY_STOPPING false
#endif

#ifndef KSP_PD_L1
    #define KSP_PD_L1 false
#endif

constexpr unsigned int INVALID_PATH_PARENT = std::numeric_limits<unsigned int>::max();
constexpr unsigned int INVALID_PATH_ALPHA  = std::numeric_limits<unsigned int>::max();

struct Path
{
    //! The total weight of the path
    w_type length = INFINITE_DISTANCE;
    //! Id of the path, this path is deviated from
    unsigned int parent_index = INVALID_PATH_PARENT;
    //! branching_node_index
    unsigned int alpha = INVALID_PATH_ALPHA;
    //! The actual path as a list of nodes
    std::vector<NODE_ID> p;

    Path() = default;
    Path(const Path& p) = default;
    Path(Path&& p) = default;
    Path& operator=(Path&&) = default;

    Path(const w_type length, const unsigned int parent_index, const unsigned int alpha, std::vector<NODE_ID> p)
        : length(length),
          parent_index(parent_index),
          alpha(alpha),
          p(std::move(p))
    {}

    bool operator<(const Path& r) const
    {
        const auto p_size = p.size();
        const auto r_p_size = r.p.size();
        return std::tie(length, p_size, r.alpha) < std::tie(r.length, r_p_size, alpha);
    }

    void print_path()
    {
        std::cout << "path length: " << length << ", parent index: " << parent_index << ", alpha: " << alpha << std::endl;
        std::cout << "path:";
        for(unsigned int i : p)
            std::cout << " " << i;
        std::cout << std::endl;
    }

#ifdef TECHNIQUE_BOUND
    std::vector<NODE_ID> getSubPath(unsigned int index)
    {
        std::vector<NODE_ID> sub_path;
        for (int i = 0; i < index; i++)
        {
            sub_path.push_back(p[i]);
        }
        return sub_path;
    }
#endif
};

/**
 * Maintains a sorted set of at most k candidates.
 */
class Candidates
{
private:
    unsigned int paths_left;

    std::multiset<Path> candidates;

#ifdef TECHNIQUE_BOUND
    w_type k_bound = INFINITE_DISTANCE;
#endif

public:
    explicit Candidates(const unsigned int k) noexcept
        : paths_left(k)
    {}

    /**
     * Checks if the candidate list is empty.
     * @return True if empty, false otherwise.
     */
    bool empty() const noexcept
    {
        return candidates.empty();
    }

    /**
     * Returns the length of the current k-th candidate.
     * @return The length of the current k-th candidate.
     */
    w_type get_length_threshold() const noexcept
    {
        return candidates.size() < paths_left ? INFINITE_DISTANCE : candidates.crbegin()->length;
    }

#ifdef TECHNIQUE_BOUND
    w_type getKBound() const noexcept
    {
        return candidates.size() != paths_left ? k_bound : candidates.crbegin()->length;
    }

    void setKBound(w_type k_bound) noexcept
    {
        this->k_bound=k_bound;
    }
#endif

    /**
     * Adds path p as a new candidate if it is at least as good as the current k path.
     * The candidate list is than shorted to only contain a many candidates as needed.
     * @param p
     */
    void insert(Path&& p) noexcept
    {
        if(candidates.size() < paths_left || p < *candidates.rbegin())
        {
            if(candidates.size() == paths_left)
                candidates.erase(--candidates.end());

            candidates.insert(std::move(p));
        }
    }

    void insert(std::vector<Path>& paths) noexcept
    {
        for(auto& p : paths)
            insert(std::move(p));
    }

    /**
     * Returns the best candidate and removes it from the list.
     * @return The best candidate from the list.
     */
    Path get_best() noexcept
    {
        assert(!candidates.empty());
        Path temp = *candidates.begin();
        candidates.erase(candidates.begin());
        paths_left--;
        return temp;
    }
};

template<typename GraphType, unsigned int num_threads, bool pps = false>
class KSP
{
protected:
#ifdef PROF
    ProfilerStatistics profiler_statistics;
#endif
    KSPGraph<GraphType> orig;
    std::vector<Path> paths;
    Candidates candidates;
    const unsigned int k;
    static constexpr bool parallel_ps = num_threads > 1 && pps;

public:
    explicit KSP(const GraphType& g, const unsigned int k) noexcept
            : orig(g), paths(0), candidates(k), k(k)
    {}

#ifdef TECHNIQUE_BOUND
    void setKBound(w_type k_bound)
    {
        candidates.setKBound(k_bound);
    }
#endif

    void calcNodesAndEdges(std::set<NODE_ID> &nodes,std::set<std::pair<NODE_ID,NODE_ID>> &edges)
    {
        for(unsigned int i = 0; i < paths.size(); i++)
        {
            if(paths[i].length == INFINITE_DISTANCE)
                break;
            
            for(unsigned int j = 0; j < paths[i].p.size(); j++){
                nodes.insert(paths[i].p[j]);
                if(j<(paths[i].p.size()-1)){
                    edges.insert(std::make_pair(paths[i].p[j],paths[i].p[j+1]));
                }
                
            }
        }
    }

    double get_remain_edge_ratio() { return orig.get_remain_edge_ratio(); }
    
    void print_paths()
    {
        std::cout << "printing paths" << std::endl;
        for(unsigned int i = 0; i < paths.size(); i++)
        {
            if(paths[i].length == INFINITE_DISTANCE)
                break;

            std::cout << "Path " << i << ":length " << paths[i].length << ":";

            for(unsigned int j = 0; j < paths[i].p.size(); j++)
                std::cout << " " << paths[i].p[j];

            std::cout << std::endl;
        }
    }

    void print_compact_paths()
    {
        CSRGraph* csr_compact=orig.get_original_graph().get_csr_graph();
        std::cout << "printing paths" << std::endl;
        for(unsigned int i = 0; i < paths.size(); i++)
        {
            if(paths[i].length == INFINITE_DISTANCE)
                break;

            std::cout << "Path " << i << ":length " << paths[i].length << ":";

            for(unsigned int j = 0; j < paths[i].p.size(); j++)
                std::cout << " " << csr_compact->k_bound_nodes[paths[i].p[j]];

            std::cout << std::endl;
        }
    }

    void print_compact_paths(NODE_ID* original_2_compact, bool* is_k_bound_nodes, const NODE_ID original_num_nodes) {
        std::cout << "printing paths" << std::endl;
        for (unsigned int i = 0; i < paths.size(); i++) {
            if (paths[i].length == INFINITE_DISTANCE) break;

            std::cout << "Path " << i << ":length " << paths[i].length << ":";

            for (unsigned int j = 0; j < paths[i].p.size(); j++) {
                for (NODE_ID node = 0; node < original_num_nodes; node++) {
                    if (is_k_bound_nodes[node] && original_2_compact[node] == paths[i].p[j]) {
                        std::cout << " " << node;
                    }
                }
            }

            std::cout << std::endl;
        }
    }

    double get_avg_path_size()
    {
        double avg_path_size = 0;

        if(paths.size() > 0)
        {
            for(unsigned int i = 0; i < paths.size(); i++)
            {
                avg_path_size += paths[i].p.size();
            }

            avg_path_size = avg_path_size / paths.size();
        }

        return avg_path_size;
    }


#ifdef PROF
    std::string return_stats()
    {
        return profiler_statistics.return_stats();
    }
#endif
};


template<typename GraphType, unsigned int num_threads, bool pps = false>
class KSPEdgeSwap
{
protected:
#ifdef PROF
    ProfilerStatistics profiler_statistics;
#endif
    KSPGraphEdgeSwap<GraphType> orig;
    std::vector<Path> paths;
    Candidates candidates;
    const unsigned int k;
    static constexpr bool parallel_ps = num_threads > 1 && pps;

public:
    explicit KSPEdgeSwap(const GraphType& g, const unsigned int k) noexcept
            : orig(g), paths(0), candidates(k), k(k)
    {}

#ifdef TECHNIQUE_BOUND
    void setKBound(w_type k_bound)
    {
        candidates.setKBound(k_bound);
    }
#endif

    void print_paths()
    {
        std::cout << "printing paths" << std::endl;
        for(unsigned int i = 0; i < paths.size(); i++)
        {
            if(paths[i].length == INFINITE_DISTANCE)
                break;

            std::cout << "Path " << i << ":length " << paths[i].length << ":";

            for(unsigned int j = 0; j < paths[i].p.size(); j++)
                std::cout << " " << paths[i].p[j];

            std::cout << std::endl;
        }
    }

#ifdef PROF
    std::string return_stats()
    {
        return profiler_statistics.return_stats();
    }
#endif
};

#endif  // _KSP_H
