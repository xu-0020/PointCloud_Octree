#include <vector>
#include "OctreeNode.h"
using namespace std;


class Octree {
private:
    OctreeNode* root;
    int maxDepth;
    int maxPointsPerNode;
    int depthAdjustmentFactor;

    // Function to determine the child index for a point
    int getOctant(const Point& origin, Point& point);

    // Recursive function to insert a point into the tree
    void insert(OctreeNode* node, Point& point, int depth);
    
    void subdivideAndInsert(OctreeNode* node, Point& point, int depth);

    // Function to visualize built node in the Octree
    void visualizeNode(OctreeNode* node, int level);

public:
    Octree(Point origin, float initialSize, int maxDepth, int maxPoints, int depthFactor) : maxDepth(maxDepth), maxPointsPerNode(maxPoints), depthAdjustmentFactor(depthFactor) {
        root = new OctreeNode(origin, initialSize);
    }

    ~Octree() {
        delete root;
    }

    void insert(Point& point) {
        insert(root, point, 0);
    }

    // Function to visualize octree
    void visualize() {
        visualizeNode(root, 0);
    }
};