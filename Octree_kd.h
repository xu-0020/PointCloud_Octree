
#ifndef OCTREEKD_H
#define OCTREEKD_H

#include <fstream>
#include "OctreeNode_kd.h"
using namespace std;


class Octree_kd {
private:
    OctreeNode_kd* root;
    int maxDepth;
    int maxPointsPerNode;
    int minPointsPerNode;

    // Function to determine the child index for a point
    int getOctant(const Point& origin, Point& point);

    // Recursive function to insert a point into the tree
    void insert(OctreeNode_kd* node, Point& point, int depth);
    
    void subdivideAndInsert(OctreeNode_kd* node, Point& point, int depth);

    // Function to merge underpopulated leaf nodes
    void mergeUnderpopulatedNodes(OctreeNode_kd* node, int depth, int startDepth);

    // Function to visualize built node in the Octree
    void visualizeNode(OctreeNode_kd* node, int level, ofstream& outFile);

    // Function to calculate bounding box for each octants
    Bounds calculateChildBounds(Bounds& parentBounds, int octant);

    // Range Query function
    void rangeQuery(Bounds& queryRange, vector<Point>& results, OctreeNode_kd* node, int depth = 0);

    // Function to create R-trees for each leaf node
    void initializeRTrees(OctreeNode_kd* node);

    // Function to create Kd-trees for each leaf node
    void initializeKdTrees(OctreeNode_kd* node);

public:
    Octree_kd(Bounds bound, int maxDepth, int maxPoints, int minPoints) : maxDepth(maxDepth), maxPointsPerNode(maxPoints), minPointsPerNode(minPoints) {
        root = new OctreeNode_kd(bound);
    };

    ~Octree_kd() {
        delete root;
    };

    void insert(Point& point) {
        insert(root, point, 0);
    };

    // Function to optimize octree 
    void trim(int startDepth) {
        mergeUnderpopulatedNodes(root, 0, startDepth);
    };

    // Function to visualize octree
    void visualize(string trialNum) {
        ofstream outFile(trialNum + " octree.txt");
        if (!outFile.is_open()) {
            cerr << "Failed to open file for writing.\n";
            return;
        }

        visualizeNode(root, 0, outFile);
        outFile.close();
    };

    void executeRangeQuery(Bounds& queryRange, vector<Point>& results) {
        rangeQuery(queryRange, results, root);
    };
    
    // Search leaf nodes and build R-trees
    void buildRtrees() {
        initializeRTrees(root);
    }

    void buildKdtrees() {
        initializeKdTrees(root);
    }
};



#endif  // BOUNDS_H