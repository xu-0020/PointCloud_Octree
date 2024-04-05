#include <vector>
#include <future>

#include "Point.h"
#include "Bound.h"
#include "RTree.h"

typedef RTree<Point, float, 3> RTreePoints;
typedef Point ValueType;

struct OctreeNode {
    Bounds bound;                               // Bounding box of the node
    vector<int> points;                         // indexs of Points contained in this node
    OctreeNode* children[8] = {nullptr};        // Pointers to octants
    RTreePoints* rtree = nullptr;               // R-tree in this leaf node

    OctreeNode() {}

    OctreeNode(const Bounds& b) : bound(b) {}

    ~OctreeNode() {
        for (auto& child : children) {
            delete child;
            child = nullptr;
        }
        delete rtree; 
        rtree = nullptr;
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