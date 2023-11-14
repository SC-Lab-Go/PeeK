#ifndef MAX_HEAP_NO_DUPLICATE_H
#define MAX_HEAP_NO_DUPLICATE_H

#include <algorithm>
#include <iterator>
#include <map>
#include <queue>
#include <set>
#include <vector>

#define EPSILON 0.0000001
#define FLOATEQUAL(a, b) (((a) >= ((b)-EPSILON)) && ((a) <= ((b) + EPSILON)))

template <typename T>
class MaxHeapNoDuplicate {
   private:
    std::priority_queue<T> pq;
    std::set<T> unique_data;
    unsigned int k;

   public:
    MaxHeapNoDuplicate(const unsigned int k) : k(k) {}

    void push(T t) {
        if (pq.size() < k) {
            if (!isExist(t)) {
                pq.push(t);
                unique_data.insert(t);
            }
        } else {
            T top_data = pq.top();
            if (t < top_data && !isExist(t)) {
                // unique_data.erase(top_data);
                pq.pop();

                pq.push(t);
                unique_data.insert(t);
            }
        }
    }

    void push2(T t) {
        if (pq.size() < k) {
            if (unique_data.find(t) == unique_data.end()) {
                pq.push(t);
                unique_data.insert(t);
            }
        } else {
            T top_data = pq.top();
            if (t < top_data && unique_data.find(t) == unique_data.end()) {
                // unique_data.erase(top_data);
                pq.pop();

                pq.push(t);
                unique_data.insert(t);
            }
        }
    }

    void pop() { pq.pop(); }

    T top() { return pq.top(); }

    unsigned int size() { return pq.size(); }

    bool empty() { return pq.empty(); }

    friend std::ostream &operator<<(std::ostream &output, const MaxHeapNoDuplicate &max_heap) {
        output << "max_heap:[";
        for (auto elem : max_heap.unique_data) {
            output << elem << " ";
        }
        output << "]";
        return (output);
    }

   private:
    bool isExist(T t) {
        for (auto elem : unique_data) {
            if (FLOATEQUAL(elem, t)) {
                return true;
            }
        }
        return false;
    }
};

template <typename T>
class KBoundPath {
   private:
    std::priority_queue<T> pq;
    std::map<w_type, std::set<std::set<NODE_ID>>> weight_path;
    unsigned int k;

   public:
    KBoundPath(const unsigned int k) : k(k) {}

    void push(T t, std::set<NODE_ID> &path) {
        if (getK() < k) {
            if (!isPathExist(t, path)) {
                pq.push(t);
                weight_path[t].insert(path);
            }
        } else {
            T top_data = pq.top();
            if (t < top_data && !isPathExist(t, path)) {
                weight_path[top_data].erase(weight_path[top_data].begin());
                if (0 == weight_path[top_data].size()) {
                    weight_path.erase(top_data);
                }
                pq.pop();

                pq.push(t);
                weight_path[t].insert(path);
            }
        }
    }

    std::set<NODE_ID> getKBoundPathNodes(unsigned int k, w_type &k_bound, unsigned int &actual_k) {
        std::set<std::set<NODE_ID>> paths;
        paths.clear();
        std::set<NODE_ID> nodes;
        bool is_k_bound_node = false;
        for (auto elem : weight_path) {
            k_bound = elem.first;
            for (auto path : elem.second) {
                paths.insert(path);
                nodes.insert(path.begin(), path.end());
                if (paths.size() == k) {
                    is_k_bound_node = true;
                    break;
                }
            }
            if (is_k_bound_node) {
                break;
            }
        }
        actual_k = paths.size();
        return nodes;
    }

    void pop() { pq.pop(); }

    T top() { return pq.top(); }

    unsigned int size() { return pq.size(); }

    bool empty() { return pq.empty(); }

    friend std::ostream &operator<<(std::ostream &output, const KBoundPath &max_heap) {
        output << "max_heap:[" << std::endl;

        for (auto elem : max_heap.weight_path) {
            output << "weight:" << elem.first << ",path:";

            for (auto path : elem.second) {
                output << "[";
                for (auto node : path) {
                    output << node << " ";
                }
                output << "] ";
            }
            output << std::endl;
        }

        output << "]";

        return (output);
    }

   private:
    bool isPathExist(T t, std::set<NODE_ID> &path) {
        for (auto elem : weight_path[t]) {
            if (elem == path) {
                return true;
            }
        }
        return false;
    }

    unsigned int getK() {
        std::set<std::set<NODE_ID>> paths;
        paths.clear();
        for (auto elem : weight_path) {
            for (auto path : elem.second) {
                paths.insert(path);
            }
        }
        return paths.size();
    }
};

struct nodeWeight {
    NODE_ID id;
    w_type distance;
};

bool cmp(nodeWeight a, nodeWeight b) {
    if (a.distance < b.distance) {
        return true;
    }
    return false;
};

#endif
