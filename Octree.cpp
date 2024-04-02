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


int Octree::getOctant(const Point& origin, Point& point) {  
    int octant = 0;
    if (point.x >= origin.x) octant |= 4;       // The third bit (from the right, 0-indexed) of octant is set 
    if (point.y >= origin.y) octant |= 2;       // The second bit (from the right, 0-indexed) of octant is set 
    if (point.z >= origin.z) octant |= 1;       // The first bit (from the right, 0-indexed) of octant is set 
    return octant;      // 8 possible results from 3 bits 000 to 111 representing 8 octants
}

void Octree::insert(OctreeNode* node, Point& point, int depth) {
    if (node->isLeaf()) {
        if (node->points.size() < maxPointsPerNode || depth >= maxDepth) {
            node->points.push_back(point);
            return;
        }
        // Subdivide the node if it exceeds capacity and within depth limit
        for (Point p : node->points) {
            subdivideAndInsert(node, p, depth);     // Reinsert existing points in the current node
        }
        node->points.clear();
        subdivideAndInsert(node, point, depth);     // Insert target point
    }
    else {
        subdivideAndInsert(node, point, depth);
    }
}

void Octree::subdivideAndInsert(OctreeNode* node, Point& point, int depth) {
    int octant = getOctant(node->bound.getCenter(), point);
    
    if (node->children[octant] == nullptr) {
        Bounds childBounds = calculateChildBounds(node->bound, octant);
        node->children[octant] = new OctreeNode(childBounds);
    }
    insert(node->children[octant], point, depth + 1);
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


void Octree::mergeUnderpopulatedNodes(OctreeNode* node, int depth, int startDepth) {
    if (!node || node->isLeaf()) return; 

    // bottom-up merge
    for (int i = 0; i < 8; i++) {
        mergeUnderpopulatedNodes(node->children[i], depth + 1, startDepth);
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

    if (allChildrenAreLeaves && totalPoints <= maxPointsPerNode * 1.2) {    // if all children combined have less than 1.5*max threshold, merge
        vector<Point> mergedPoints;
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                mergedPoints.insert(mergedPoints.end(), node->children[i]->points.begin(), node->children[i]->points.end());
                delete node->children[i];
                node->children[i] = nullptr;
            }
        }
        node->points.swap(mergedPoints);  
        node->convertToLeaf();  
    }
    
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

        vector<Point> mergedPoints;
        vector<OctreeNode*> newLeafNodes;

        for (const auto& [pointCount, index] : childPointCounts) {
            if (!mergedPoints.empty() && mergedPoints.size() + pointCount > maxPointsPerNode) {    // Newly merged exceed the max limit
                OctreeNode* newLeaf = new OctreeNode();
                newLeaf->points.swap(mergedPoints);  // Move the merged points to the new leaf node, use swap to save memory
                newLeaf->convertToLeaf();
                newLeafNodes.push_back(newLeaf);    // Add node
            }

            auto& childPoints = node->children[index]->points;
            mergedPoints.insert(mergedPoints.end(),
                        make_move_iterator(childPoints.begin()),
                        make_move_iterator(childPoints.end()));     // Use move iterators to save memory
            delete node->children[index]; 
            node->children[index] = nullptr;
        }

        // Handle remaining merged points
        if (!mergedPoints.empty()) {
            OctreeNode* newLeaf = new OctreeNode();
            newLeaf->points.swap(mergedPoints);
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

    /*
    // Not improving efficiency
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
                    mergeBase->points.insert(mergeBase->points.end(),
                                            make_move_iterator((*it)->points.begin()),
                                            make_move_iterator((*it)->points.end()));

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


void Octree::rangeQuery(Bounds& queryRange, vector<Point>& results, OctreeNode* node, int depth, float start=-1, float end=-1) {
    // Check if the current node's bounds intersect with the query range
    if (!node->bound.intersects(queryRange)) {
        return;
    }

    if (node->isLeaf()) {
        // If it's a leaf node, query the R-tree
        if (node->hrtree) {
            node->hrtree->HRTreeSearch(start, end, results, queryRange);
        } 
    } 
    else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                rangeQuery(queryRange, results, node->children[i], depth + 1);
            }
        }
    }
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



// Function to create R-trees for each leaf node with thread limitation
void Octree::initializeRTrees(OctreeNode* node, vector<future<void>>& futures) {
    if (node->isLeaf()) {

        // Check if reached the maximum number of concurrent tasks
        if (futures.size() >= max_concurrent_tasks) {
            // Wait for at least one task to complete
            bool taskCompleted = false;
            while (!taskCompleted) {
                for (auto it = futures.begin(); it != futures.end(); ) {
                    auto& fut = *it;
                    if (fut.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
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
        futures.push_back(async(launch::async, [this, node]() {
            // Regenerate bounds for the leaf node to ensure tight fitting
            node->bound = calculateBoundsForPoints(node->points);

            node->hrtree = new HRTree(node->points);
            node->points.clear(); // Clear points after moving them to R-tree to save memory
        }));
    } else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                initializeRTrees(node->children[i], futures);
            }
        }
    }
}

Bounds Octree::calculateBoundsForPoints(const vector<Point>& points) {
    if (points.empty()) return Bounds();  // Return default bounds if no points

    Point min = points[0], max = points[0];
    for (const Point& point : points) {
        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);

        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);
    }
    return Bounds(min, max);
}
