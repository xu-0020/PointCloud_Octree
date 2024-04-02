#include <vector>
#include <future>

#include "HRTree.cpp"
typedef RTree<Point, float, 3> RTreePoints;
typedef Point ValueType;

struct OctreeNode {
    Bounds bound;                               // Bounding box of the node
    vector<Point> points;                       // Points contained in this node
    OctreeNode* children[8] = {nullptr};        // Pointers to octants
    HRTree* hrtree = nullptr;               // R-tree in this leaf node

    OctreeNode() {}

    OctreeNode(const Bounds& b) : bound(b) {}

    ~OctreeNode() {
        for (auto& child : children) {
            delete child;
            child = nullptr;
        }
        delete hrtree; 
        hrtree = nullptr;
    }

    bool isLeaf() const {
        for (const auto& child : children) {
            if (child != nullptr) {
                return false; 
            }
        }
        return true; 
    }


    void convertToLeaf() {
        for (int i = 0; i < 8; i++) {
            if (children[i]) {
                delete children[i];
                children[i] = nullptr;
            }
        }
    }
};