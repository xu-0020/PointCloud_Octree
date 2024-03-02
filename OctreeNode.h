#include <vector>
#include "Point.h"


struct OctreeNode {
    Point origin;                               // Centroid point
    float size;                                 // Length of the side of the cubic node
    vector<Point> points;                       // Points contained in this node
    OctreeNode* children[8] = {nullptr};        // Pointers to octants

    OctreeNode(Point o, float s) : origin(o), size(s) {}

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