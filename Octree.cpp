#include <iostream>
#include <algorithm>

#include "Octree.h"
#include <string.h>

const size_t max_concurrent_tasks = 30;     // concurrency control for building r-trees


// R-tree functions
void RInsert(RTreePoints* tree, const vector<Point>& points){
    if (points.size()==0){
        return;
    }
    for (const Point& point : points){
        // Point* pointer = new Point();
        tree->Insert(point.arr, point.arr, point);
    }
    // auto list = tree->ListTree();
    // int counter = 0;
    // for (auto aabb : list) {
    // cout << "TreeList [" << counter++ << "]: "
    //     << aabb.m_min[0] << ", "
    //     << aabb.m_min[1] << ", "
    //     << aabb.m_min[2] << "; "
    //     << aabb.m_max[0] << ", "
    //     << aabb.m_max[1] << ", "
    //     << aabb.m_max[2] << endl;
    // }
}

bool MySearchCallback(ValueType point)
{
//   point.print_point();
    return true; // keep going
}

void RSearch(RTreePoints* tree, std::vector<Point>& results, Bounds& queryRange){
    auto list = tree->ListTree();
    int counter = 0;
    // for (auto aabb : list) {
    // cout << "TreeList [" << counter++ << "]: "
    //     << aabb.m_min[0] << ", "
    //     << aabb.m_min[1] << ", "
    //     << aabb.m_min[2] << "; "
    //     << aabb.m_max[0] << ", "
    //     << aabb.m_max[1] << ", "
    //     << aabb.m_max[2] << endl;
    // }
    // queryRange.min.print_point();
    // queryRange.max.print_point();
    // cout<<endl<<endl;
    tree->Search(queryRange.min.arr, queryRange.max.arr, results, MySearchCallback);
    // for (Point point : results){
    //     point.print_point();
    // }
    return;
}
// End of R-tree functions



vector<OctreeNode*> Octree::createLeafNodesMultithreaded(const vector<Point>& dataPoints, const vector<pair<int, uint64_t>>& mortonCode) {
    const size_t numThreads = 32;
    const size_t segmentSize = ceil(mortonCode.size() / static_cast<double>(numThreads));

    uint64_t numberOfGrid = pow(8, maxDepth);    // Number of grids at the maximum depth
    uint64_t gridSize = (pow(2, 32) - 1) / numberOfGrid;     // Grid size for the current depth

    vector<future<vector<OctreeNode*>>> futures;

    for (size_t i = 0; i < numThreads && i * segmentSize < mortonCode.size(); ++i) {
        size_t start = i * segmentSize;
        size_t end = min(start + segmentSize, mortonCode.size());

        futures.push_back(async(launch::async, [&, start, end] {
            return createLeafNodesSegment(dataPoints, mortonCode, start, end, gridSize);
        }));
    }

    vector<OctreeNode*> leafNodes;
    for (auto& fut : futures) {
        auto segmentLeafNodes = fut.get();
        leafNodes.insert(leafNodes.end(), segmentLeafNodes.begin(), segmentLeafNodes.end());
    }

    return leafNodes;
}


vector<OctreeNode*> Octree::createLeafNodesSegment(const vector<Point>& dataPoints, const vector<pair<int, uint64_t>>& mortonCode, size_t start, size_t end, uint64_t gridSize) {
    vector<OctreeNode*> segmentLeafNodes;

    while (start < end) {
        // Morton code at the start pointer represents the beginning of the current grid
        uint64_t startGridCode = mortonCode[start].second - (mortonCode[start].second % gridSize);

        // Collect indices of points in the current grid
        vector<int> pointIndices;

        while (start < end && (mortonCode[start].second - (mortonCode[start].second % gridSize)) == startGridCode) {
            pointIndices.push_back(mortonCode[start].first);
            start++;  
        }

        // Create a leaf node for the current grid if it contains points
        if (!pointIndices.empty()) {
            Bounds bounds = calculateBoundsForPoints(dataPoints, pointIndices);

            // Create a new leaf node
            OctreeNode* leafNode = new OctreeNode(bounds);
            leafNode->points.swap(pointIndices); 
            segmentLeafNodes.push_back(leafNode);
        }

    }

    return segmentLeafNodes;
}



void Octree::buildFromLeafNodes(vector<OctreeNode*>& nodes, int currentDepth) {
    if (nodes.size() <= 1 || currentDepth == 0) {
        // If there's only one node left or we've reached the maximum depth, this is the root
        if (!nodes.empty()) {
            OctreeNode* head = nodes.front();
            Bounds bounds = root->bound;
            delete root;
            root = head;
            head->bound = bounds;
        }
        return;
    }

    vector<OctreeNode*> parentNodes;
    uint32_t parentCount = ceil(nodes.size() / 8.0);  // Calculate the number of parent nodes needed

    for (int i = 0; i < parentCount; i++) {     // Loop over each parent node
        Bounds parentBounds;
        vector<int> parentPointIndices;
        bool firstNode = true;
        bool isInternal = false;    // Flag to check if any child is an internal node

        for (int j = 0; j < 8 && i * 8 + j < nodes.size(); j++) {  // Each parent can have up to 8 children
            int childIndex = i * 8 + j;
            OctreeNode* childNode = nodes[childIndex];

            // Update parent bounds to include child bounds
            if (firstNode) {
                parentBounds = childNode->bound;
                firstNode = false;
            } else {
                parentBounds.update(childNode->bound.min);
                parentBounds.update(childNode->bound.max);
            }

            // Check if any child is an internal node
            if (!childNode->isLeaf()) {
                isInternal = true;
            } else {
                // Merge child points into parent
                parentPointIndices.insert(parentPointIndices.end(), childNode->points.begin(), childNode->points.end());
            }
        }

        // Create a new parent node
        OctreeNode* parentNode = new OctreeNode(parentBounds);

        // Internal: Assign children to the parent node if any child is internal or total size exceeds the limit
        if (isInternal || parentPointIndices.size() > maxPointsPerNode) {
            for (int j = 0; j < 8 && i * 8 + j < nodes.size(); j++) {
                parentNode->children[j] = nodes[i * 8 + j];
            }
        }  
        // Leaf: Assign points to the parent node and delete the child nodes
        else {
            parentNode->points.swap(parentPointIndices);
            parentNode->convertToLeaf();

            // Delete the child nodes that have been merged
            for (int j = 0; j < 8 && i * 8 + j < nodes.size(); j++) {
                delete nodes[i * 8 + j];  // Delete the child node
                nodes[i * 8 + j] = nullptr; 
            }
        }

        parentNodes.push_back(parentNode);
    }
    
    // Recursively build the tree
    buildFromLeafNodes(parentNodes, currentDepth - 1);
}







void Octree::visualizeNode(OctreeNode* node, int level, ofstream& outFile) {
    if (!node) return;

    for (int i=0; i<level; i++) outFile << "  ";   // Indentation

    if (node->isLeaf()) {
        outFile << "Level " << level << ": Leaf node with " << node->points.size() << " points\n";
    }
    else {
        outFile << "Level " << level << ": Internal node\n";
        for (int i=0; i<8; i++) {
            visualizeNode(node->children[i], level+1, outFile);
        }
    }
}


void Octree::divideOverpopulatedNodes(OctreeNode* node, vector<future<void>>& futures, int depth, const vector<Point>& dataPoints) {
    if (!node) return;

    if (node->isLeaf()) {
        // Check if reached the maximum number of concurrent tasks
        if (futures.size() >= max_concurrent_tasks) {
            // Wait for at least one task to complete
            bool taskCompleted = false;
            while (!taskCompleted) {
                for (auto it = futures.begin(); it != futures.end(); ) {
                    auto& fut = *it;
                    if (fut.wait_for(chrono::seconds(0)) == future_status::ready) {
                        fut.get();              // Get the result to clear any stored exception
                        it = futures.erase(it);     // Remove the completed future
                        taskCompleted = true;
                        break;           // Break the loop as we only need one task to complete
                    } else {
                        it++;
                    }
                }
            }
        }

        futures.push_back(async(launch::async, [this, node, depth, &dataPoints]() {
            if (node->isLeaf() && node->points.size() > maxPointsPerNode * 1.2) {
            // Handle Overpopulated leaf nodes at bottom level (further divide)
                for (int index : node->points) {
                    subdivideAndInsert(node, index, depth, dataPoints);  // Reinsert existing points in the current node
                }
                node->points.clear();  // Remove points as they are moved to a lower level
            }}));
    } 
    else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                divideOverpopulatedNodes(node->children[i], futures, depth + 1, dataPoints);
            }
        }
    }
}




void Octree::mergeUnderpopulatedNodes(OctreeNode* node, int depth, const vector<Point>& dataPoints) {
    if (!node || node->isLeaf()) return;

    // Perform bottom-up merge
    for (int i = 0; i < 8; i++) {
        mergeUnderpopulatedNodes(node->children[i], depth + 1, dataPoints);
    }
    
    bool allChildrenAreLeaves = true;
    int totalPoints = 0;
    for (int i = 0; i < 8; i++) {   // check if children are all leaves
        if (node->children[i]) {
            if (!node->children[i]->isLeaf()) {
                allChildrenAreLeaves = false;
                break;  // break if one child is not leaf (already checked-exceed max)
            } else {
                totalPoints += node->children[i]->points.size();
            }
        }
    }

    if (allChildrenAreLeaves && totalPoints <= maxPointsPerNode * 1.5) {    // if all children combined have less than 1.5*max threshold, merge
        vector<int> mergedPointIndices;
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                mergedPointIndices.insert(mergedPointIndices.end(), node->children[i]->points.begin(), node->children[i]->points.end());
                delete node->children[i];
                node->children[i] = nullptr;
            }
        }
        node->points.swap(mergedPointIndices);  
        node->convertToLeaf();  
    }
    /*
    // Not Improving the performance
    else if (allChildrenAreLeaves) {    // if some children combined have less than max threshold, combine them
        // sort children
        vector<pair<int, int>> childPointCounts;
        childPointCounts.reserve(8);    // Reserve memory to avoid reallocations
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                childPointCounts.emplace_back(node->children[i]->points.size(), i);;
            }
        }   
        sort(childPointCounts.begin(), childPointCounts.end()); // Sort the vector by point count in ascending order

        vector<int> mergedPointIndices;
        vector<OctreeNode*> newLeafNodes;

        for (const auto& [pointCount, index] : childPointCounts) {
            if (!mergedPointIndices.empty() && mergedPointIndices.size() + pointCount > maxPointsPerNode) {    // Newly merged exceed the max limit
                OctreeNode* newLeaf = new OctreeNode();
                newLeaf->points.swap(mergedPointIndices);  // Move the merged points to the new leaf node, use swap to save memory
                newLeaf->convertToLeaf();
                newLeafNodes.push_back(newLeaf);    // Add node
            }

            auto& childPointIndices = node->children[index]->points;
            mergedPointIndices.insert(mergedPointIndices.end(), childPointIndices.begin(), childPointIndices.end()); 
            delete node->children[index]; 
            node->children[index] = nullptr;
        }

        // Handle remaining merged points
        if (!mergedPointIndices.empty()) {
            OctreeNode* newLeaf = new OctreeNode();
            newLeaf->points.swap(mergedPointIndices);
            newLeaf->convertToLeaf();
            newLeafNodes.push_back(newLeaf);
        }

        // Reinsert any remaining children after the merged one
        for (size_t i = 0; i < newLeafNodes.size(); i++) {
            node->children[i] = newLeafNodes[i];
        }

        // Clear any remaining child pointers
        for (size_t i = newLeafNodes.size(); i < 8; i++) {
            node->children[i] = nullptr;
        }
    }
    
    else {
        vector<OctreeNode*> leafNodes;
        vector<OctreeNode*> internalNodes;
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]) {
                if (node->children[i]->isLeaf()) {
                    leafNodes.push_back(node->children[i]);
                } else {
                    internalNodes.push_back(node->children[i]);
                }
                node->children[i] = nullptr;  // Clear the child pointer
            }
        }
        
        // Sort leaf nodes by their point count in ascending order
        sort(leafNodes.begin(), leafNodes.end(), [](const OctreeNode* a, const OctreeNode* b) {
            return a->points.size() < b->points.size();
        });

        vector<OctreeNode*> nodesToReattach = internalNodes; 

        while (!leafNodes.empty()) {
            OctreeNode* mergeBase = leafNodes.front();
            leafNodes.erase(leafNodes.begin());
            bool merged = false;

            for (auto it = leafNodes.begin(); it != leafNodes.end(); ) {
                int potentialPointCount = mergeBase->points.size() + (*it)->points.size();
                if (potentialPointCount <= maxPointsPerNode) {
                    mergeBase->points.insert(mergeBase->points.end(), (*it)->points.begin(), (*it)->points.end());
                    delete *it;  // Delete the absorbed node
                    it = leafNodes.erase(it);  // Remove from leafNodes
                    merged = true;
                } else {
                    it++;
                }
            }

            if (merged) {
                mergeBase->convertToLeaf();  // Ensure mergeBase is a leaf if any merging occurred
                nodesToReattach.push_back(mergeBase);  // Add mergeBase to the list of nodes to reattach
            } else {
                nodesToReattach.push_back(mergeBase);  // Reattach mergeBase if no merging occurred
            }
        }

        // Clear and rebuild node->children with reattached nodes
        memset(node->children, 0, sizeof(node->children));      // Clear existing children pointers
        for (int i = 0; i < nodesToReattach.size() && i < 8; i++) {
            node->children[i] = nodesToReattach[i];  // Reattach nodes
        }

    }
    */
}

void Octree::subdivideAndInsert(OctreeNode* node, int pointIdx, int depth, const vector<Point>& dataPoints) {
    int octant = getOctant(node->bound.getCenter(), pointIdx, dataPoints);
    
    if (node->children[octant] == nullptr) {
        Bounds childBounds = calculateChildBounds(node->bound, octant);
        node->children[octant] = new OctreeNode(childBounds);
    }
    insert(node->children[octant], pointIdx, depth + 1, dataPoints);
}

void Octree::insert(OctreeNode* node, int pointIdx, int depth, const vector<Point>& dataPoints) {
    if (node->points.size() < maxPointsPerNode * 1.5 || depth >= maxDepth * 1.5) {
        node->points.push_back(pointIdx);
        return;
    }
    // Subdivide the node if it exceeds capacity and within depth limit
    for (int indices : node->points) {
        subdivideAndInsert(node, indices, depth, dataPoints);     // Reinsert existing points in the current node
    }
    node->points.clear();
    subdivideAndInsert(node, pointIdx, depth, dataPoints);     // Insert target point

}




void Octree::rangeQuery(Bounds& queryRange, vector<Point>& results, OctreeNode* node) {

    if (node->isLeaf()) {
        // If it's a leaf node, query the R-tree
        RSearch(node->rtree, results, queryRange);
    } 
    else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                if (node->children[i]->bound.intersects(queryRange)) {  // Check if the node's bounds intersect with the query range
                    rangeQuery(queryRange, results, node->children[i]);
                }
            }
        }
    }
}



// Function to create R-trees for each leaf node with thread limitation
void Octree::initializeRTrees(OctreeNode* node, vector<future<void>>& futures, const vector<Point>& dataPoints) {
    if (node->isLeaf()) {

        // Check if reached the maximum number of concurrent tasks
        if (futures.size() >= max_concurrent_tasks) {
            // Wait for at least one task to complete
            bool taskCompleted = false;
            while (!taskCompleted) {
                for (auto it = futures.begin(); it != futures.end(); ) {
                    auto& fut = *it;
                    if (fut.wait_for(chrono::seconds(0)) == future_status::ready) {
                        fut.get();              // Get the result to clear any stored exception
                        it = futures.erase(it);     // Remove the completed future
                        taskCompleted = true;
                        break;           // Break the loop as we only need one task to complete
                    } else {
                        it++;
                    }
                }
            }
        }

        // Launch a new task for R-tree construction in the leaf node
        futures.push_back(async(launch::async, [this, node, &dataPoints]() {
            // Regenerate bounds for the leaf node to ensure tight fitting
            node->bound = calculateBoundsForPoints(dataPoints, node->points);

            node->rtree = new RTreePoints();
            vector<Point> insertPoints;         // Retrieve points to be inserted in to the Rtree.
            for (int i=0; i<node->points.size(); i++) {
                insertPoints.push_back(dataPoints[node->points[i]]);
            }
            node->points.clear();
            RInsert(node->rtree, insertPoints);
        }));
    } else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                initializeRTrees(node->children[i], futures, dataPoints);
            }
        }
    }
}



Bounds Octree::calculateBoundsForPoints(const vector<Point>& dataPoints, const vector<int>& pointsIndices) {
    if (pointsIndices.empty()) return Bounds();  // Return default bounds if no points

    Point min = dataPoints[pointsIndices[0]], max = dataPoints[pointsIndices[0]];
    for (int idx : pointsIndices) {
        const Point& point = dataPoints[idx];

        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);

        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);
    }
    return Bounds(min, max);
}

Bounds Octree::calculateChildBounds(Bounds& parentBounds, int octant) {
    Point center = parentBounds.getCenter();
    Point min = parentBounds.min;
    Point max = parentBounds.max;
    Point childMin, childMax;

    // Calculate min and max points for the child bounds based on the octant
    if (octant & 4) {
        childMin.x = center.x;
        childMax.x = max.x;
    } else {
        childMin.x = min.x;
        childMax.x = center.x;
    }
    if (octant & 2) {
        childMin.y = center.y;
        childMax.y = max.y;
    } else {
        childMin.y = min.y;
        childMax.y = center.y;
    }
    if (octant & 1) {
        childMin.z = center.z;
        childMax.z = max.z;
    } else {
        childMin.z = min.z;
        childMax.z = center.z;
    }

    return Bounds(childMin, childMax);
}


int Octree::getOctant(const Point& origin, int pointIdx, const vector<Point>& dataPoints) {  
    int octant = 0;
    const Point& point = dataPoints[pointIdx];
    if (point.x >= origin.x) octant |= 4;       // The third bit (from the right, 0-indexed) of octant is set 
    if (point.y >= origin.y) octant |= 2;       // The second bit (from the right, 0-indexed) of octant is set 
    if (point.z >= origin.z) octant |= 1;       // The first bit (from the right, 0-indexed) of octant is set 
    return octant;      // 8 possible results from 3 bits 000 to 111 representing 8 octants
}