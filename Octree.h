#include <fstream>
#include "OctreeNode.h"
using namespace std;


class Octree {
private:
    OctreeNode* root;
    int maxDepth;
    int maxPointsPerNode;

    // Function to initialize the leaf nodes
    vector<OctreeNode*> createLeafNodes(const vector<Point>& dataPoints, const vector<pair<int, uint64_t>>& mortonCode);

    // Function to build the tree from leaf nodes
    void buildFromLeafNodes(vector<OctreeNode*>& leafNodes, int depth);

    // Function to merge underpopulated leaf nodes
    void mergeUnderpopulatedNodes(OctreeNode* node, int depth, const vector<Point>& dataPoints);

    // Function to divide overpopulated leaf nodes
    void divideOverpopulatedNodes(OctreeNode* node, int depth, const vector<Point>& dataPoints);

    // Function to visualize built node in the Octree
    void visualizeNode(OctreeNode* node, int level, ofstream& outFile);

    // Function to recalculate bounds for a set of points (Regenerate for optimization)
    Bounds calculateBoundsForPoints(const vector<Point>& dataPoints, const vector<int>& pointsIndices);

    // Range Query function
    void rangeQuery(Bounds& queryRange, vector<Point>& results, OctreeNode* node, int depth = 0);

    // Function to create R-trees for each leaf node
    void initializeRTrees(OctreeNode* node, vector<future<void>>& futures);

    Bounds calculateChildBounds(Bounds& parentBounds, int octant);

    void subdivideAndInsert(OctreeNode* node, int pointIdx, int depth, const vector<Point>& dataPoints);

    // Function to determine the child index for a point
    int getOctant(const Point& origin, int pointIdx, const vector<Point>& dataPoints);

    // Recursive function to insert a point into the tree
    void insert(OctreeNode* node, int pointIdx, int depth, const vector<Point>& dataPoints);

public:
    Octree(Bounds bound, int maxDepth, int maxPoints) : maxDepth(maxDepth), maxPointsPerNode(maxPoints) {
        root = new OctreeNode(bound);
    };

    ~Octree() {
        delete root;
    };

    void constructOctree(const vector<Point>& dataPoints, const vector<pair<int, uint64_t>>& mortonCode) {
        // Initialize bottom level leaf nodes and count for each grid
        vector<OctreeNode*> leafNodes = createLeafNodes(dataPoints, mortonCode);

        // Recursively build the tree from leaf nodes
        buildFromLeafNodes(leafNodes, maxDepth);
    }

    // Function to rebalance the octree (divide the cluster at the bottom)
    void rebalance(const vector<Point>& dataPoints) {
        divideOverpopulatedNodes(root, 0, dataPoints);
    };

    // Function to optimize octree 
    void trim(const vector<Point>& dataPoints) {
        mergeUnderpopulatedNodes(root, 0, dataPoints);
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