#include <vector>
#include "OctreeNode.h"
using namespace std;


class Octree {
private:
    OctreeNode* root;
    int maxDepth;
    int maxPointsPerNode;


    // Function to determine the child index for a point
    int getOctant(const Point& origin, Point& point);

    // Recursive function to insert a point into the tree
    void insert(OctreeNode* node, Point& point, int depth);
    
    void subdivideAndInsert(OctreeNode* node, Point& point, int depth);

public:
    Octree(Point origin, float initialSize, int maxDepth, int maxPoints) : maxDepth(maxDepth), maxPointsPerNode(maxPoints) {
        root = new OctreeNode(origin, initialSize);
    }

    ~Octree() {
        delete root;
    }

    void insert(Point& point) {
        insert(root, point, 0);
    }

};