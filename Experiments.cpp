#include "Octree.h"
#include "MortonEncoding.h"
#include "RadixSort.h"
#include "CreateBox.h"


// Build the octree + rtree
void buildTree(vector<Point> dataPoints, const vector<string>& filenames) {
    const Bounds bounds = computeTotalBoundingBox(filenames);   // Compute the bounding box of the entire dataset

    // Morton encoding, corresponding to each dataPoints, (index, mortoncode)
    vector<pair<int, uint64_t>> mortonCode = computeMortonCodesMultithreaded(dataPoints, bounds);
    
    // Sort based on Morton code
    radixSort(mortonCode);

    // Tree Variable Settings
    int maxDepth = 10;
    int maxPointsPerLeaf = 6000;

    Octree octree(maxDepth, maxPointsPerLeaf);
    octree.constructOctree(dataPoints, mortonCode);
    
    octree.rebalance(dataPoints);
    octree.trim(dataPoints);
    
    octree.buildRtrees(dataPoints);
}

void runQuery(Bounds& queryRange, vector<Point>& results, Octree octree) {
    octree.executeRangeQuery(queryRange, results);
}


