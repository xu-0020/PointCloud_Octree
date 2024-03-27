#include <fstream>
#include "OctreeNode.h"
using namespace std;


class Octree {
private:
    OctreeNode* root;
    int maxDepth;
    int maxPointsPerNode;


    // Function to merge underpopulated leaf nodes
    void mergeUnderpopulatedNodes(OctreeNode* node, int depth, const int startDepth);

    // Function to visualize built node in the Octree
    void visualizeNode(OctreeNode* node, int level, ofstream& outFile);

    // Function to calculate bounding box for each octants
    Bounds calculateChildBounds(Bounds& parentBounds, int octant);

    // Function to recalculate bounds for a set of points (Regenerate for optimization)
    Bounds calculateBoundsForPoints(const vector<Point>& points);

    // Range Query function
    void rangeQuery(Bounds& queryRange, vector<Point>& results, OctreeNode* node, int depth = 0);

    // Function to create R-trees for each leaf node
    void initializeRTrees(OctreeNode* node, vector<future<void>>& futures);

public:
    Octree(Bounds bound, int maxDepth, int maxPoints) : maxDepth(maxDepth), maxPointsPerNode(maxPoints) {
        root = new OctreeNode(bound);
    };

    ~Octree() {
        delete root;
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
        vector<future<void>> futures;
        initializeRTrees(root, futures);

        // Wait for all threads to complete
        for (auto& fut : futures) {
            fut.get();
        }
    }
};