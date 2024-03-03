#include <vector>
#include "Point.h"
#include "Bound.h"


struct OctreeNode {
    Bounds bound;                               // Bounding box of the node
    vector<Point> points;                       // Points contained in this node
    OctreeNode* children[8] = {nullptr};        // Pointers to octants

    OctreeNode(const Bounds& b) : bound(b) {}

    ~OctreeNode() {
        for (auto& child : children) {
            delete child;
            child = nullptr;
        }
    }

    bool isLeaf() const {
        return children[0] == nullptr;
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