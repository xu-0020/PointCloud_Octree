#include <fstream>
#include "OctreeNode.h"
using namespace std;


class Octree {
private:
    OctreeNode* root;
    int maxDepth;
    int maxPointsPerNode;

    // Function to create the leaf nodes multithreaded
    vector<OctreeNode*> createLeafNodesMultithreaded(const vector<Point>& dataPoints, const vector<pair<int, uint64_t>>& mortonCode);

    // Helper Function to create the leaf nodes for the thread
    vector<OctreeNode*> createLeafNodesSegment(const vector<Point>& dataPoints, const vector<pair<int, uint64_t>>& mortonCode, size_t start, size_t end, uint64_t gridSize);

    // Function to build the tree from leaf nodes
    void buildFromLeafNodes(vector<OctreeNode*>& leafNodes);

    // Function to merge underpopulated leaf nodes
    void mergeUnderpopulatedNodes(OctreeNode* node, int depth, const vector<Point>& dataPoints);

    // Function to divide overpopulated leaf nodes multithreaded
    void divideOverpopulatedNodes(OctreeNode* node, vector<future<void>>& futures, int depth, const vector<Point>& dataPoints);

    // Function to visualize built node in the Octree
    void visualizeNode(OctreeNode* node, int level, ofstream& outFile);

    // Helper Function to calculate bounds
    Bounds calculateBoundsForPoints(const vector<Point>& dataPoints, const vector<int>& pointsIndices);

    // Range Query function
    void rangeQuery(Bounds& queryRange, vector<Point>& results, OctreeNode* node);

    // Function to create R-trees for each leaf node
    void initializeRTrees(OctreeNode* node, vector<future<void>>& futures, const vector<Point>& dataPoints);

    void initializeKdTrees(OctreeNode* node, vector<future<void>>& futures, const vector<Point>& dataPoints);

    Bounds calculateChildBounds(Bounds& parentBounds, int octant);

    void subdivideAndInsert(OctreeNode* node, int pointIdx, int depth, const vector<Point>& dataPoints);

    // Function to determine the child index for a point
    int getOctant(const Point& origin, int pointIdx, const vector<Point>& dataPoints);

    // Recursive function to insert a point into the tree
    void insert(OctreeNode* node, int pointIdx, int depth, const vector<Point>& dataPoints);

public:
    Octree(int maxDepth, int maxPoints) : maxDepth(maxDepth), maxPointsPerNode(maxPoints) {
        root = new OctreeNode();
    };

    ~Octree() {
        delete root;
    };

    void constructOctree(const vector<Point>& dataPoints, const vector<pair<int, uint64_t>>& mortonCode) {
        // Initialize bottom level leaf nodes and count for each grid
        vector<OctreeNode*> leafNodes = createLeafNodesMultithreaded(dataPoints, mortonCode);

        // Recursively build the tree from leaf nodes
        buildFromLeafNodes(leafNodes);
    }

    // Function to rebalance the octree (divide the cluster at the bottom)
    void rebalance(const vector<Point>& dataPoints) {
        vector<future<void>> futures;
        divideOverpopulatedNodes(root, futures, 0, dataPoints);

        // Wait for all threads to complete
        for (auto& fut : futures) {
            fut.get();
        }
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
    void buildRtrees(const vector<Point>& dataPoints) {
        vector<future<void>> futures;
        initializeRTrees(root, futures, dataPoints);

        // Wait for all threads to complete
        for (auto& fut : futures) {
            fut.get();
        }
    }
    void buildKdtrees(const vector<Point>& dataPoints) {
        vector<future<void>> futures;
        initializeKdTrees(root, futures, dataPoints);

        // Wait for all threads to complete
        for (auto& fut : futures) {
            fut.get();
        }
    }
};