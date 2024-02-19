#include <vector>
#include "Point.h"


struct OctreeNode {
    Point origin;                               // Center point of the node
    float size;                                 // Length of the side of the cubic node
    vector<Point> points;                       // Points contained in this node
    OctreeNode* children[8] = {nullptr};        // Pointers to octants

    OctreeNode(Point o, float s) : origin(o), size(s) {}

    ~OctreeNode() {
        for (auto& child : children) {
            delete child;
        }
    }

    bool isLeaf() const {
        return children[0] == nullptr;
    }
};