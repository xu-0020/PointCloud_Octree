#include <fstream>
#include "OctreeNode.h"
using namespace std;


class Octree {
private:
    OctreeNode* root;
    int maxDepth;
    int maxPointsPerNode;
    int minPointsPerNode;
    int depthAdjustmentFactor;

    // Function to determine the child index for a point
    int getOctant(const Point& origin, Point& point);

    // Recursive function to insert a point into the tree
    void insert(OctreeNode* node, Point& point, int depth);
    
    void subdivideAndInsert(OctreeNode* node, Point& point, int depth);

    // Function to merge underpopulated leaf nodes
    void mergeUnderpopulatedNodes(OctreeNode* node, int depth, int startDepth);

    // Function to visualize built node in the Octree
    void visualizeNode(OctreeNode* node, int level, ofstream& outFile);

public:
    Octree(Point origin, float initialSize, int maxDepth, int maxPoints, int minPoints, int depthFactor) : maxDepth(maxDepth), maxPointsPerNode(maxPoints), minPointsPerNode(minPoints), depthAdjustmentFactor(depthFactor) {
        root = new OctreeNode(origin, initialSize);
    }

    ~Octree() {
        delete root;
    }

    void insert(Point& point) {
        insert(root, point, 0);
    }

    void trim(int startDepth) {
        mergeUnderpopulatedNodes(root, 0, startDepth);
    }

    // Function to visualize octree
    void visualize(string trialNum) {
        ofstream outFile(trialNum + " octree.txt");
        if (!outFile.is_open()) {
            cerr << "Failed to open file for writing.\n";
            return;
        }

        visualizeNode(root, 0, outFile);
        outFile.close();
    }
};