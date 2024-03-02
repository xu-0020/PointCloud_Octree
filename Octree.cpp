#include <iostream>
#include <algorithm>

#include "Octree.h"


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
    int octant = getOctant(node->origin, point);
    
    if (node->children[octant] == nullptr) {
        float newSize = node->size / 2.0f;
        // Move origin to the determined octant
        Point newOrigin = node->origin;     
        if (octant & 4) newOrigin.x += newSize / 2.0f;
        else newOrigin.x -= newSize / 2.0f;
        if (octant & 2) newOrigin.y += newSize / 2.0f;
        else newOrigin.y -= newSize / 2.0f;
        if (octant & 1) newOrigin.z += newSize / 2.0f;
        else newOrigin.z -= newSize / 2.0f;
        
        node->children[octant] = new OctreeNode(newOrigin, newSize);
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
    
    if (depth >= startDepth) {
        bool allChildrenAreLeaves = true;
        int totalPoints = 0;
        for (int i = 0; i < 8; ++i) {   // check if children are all leaves
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
            vector<Point> mergedPoints;
            for (int i = 0; i < 8; ++i) {
                if (node->children[i]) {
                    mergedPoints.insert(mergedPoints.end(), node->children[i]->points.begin(), node->children[i]->points.end());
                    delete node->children[i];
                    node->children[i] = nullptr;
                }
            }
            node->points = mergedPoints;  
            node->convertToLeaf();  
        }
        
        else if (allChildrenAreLeaves) {    // if some children combined have less than max threshold, combine them
            // sort children
            vector<pair<int, int>> childPointCounts;
            for (int i = 0; i < 8; ++i) {
                if (node->children[i]) {
                    childPointCounts.push_back(pair(node->children[i]->points.size(), i));
                }
            }   
            sort(childPointCounts.begin(), childPointCounts.end()); // Sort the vector by point count in ascending order

            vector<Point> mergedPoints;
            vector<OctreeNode*> newLeafNodes;

            for (auto& [pointCount, index] : childPointCounts) {
                if (mergedPoints.size() + pointCount > maxPointsPerNode * 1.2) {    // Newly merged exceed the max limit
                    OctreeNode* newLeaf = new OctreeNode(node->origin, node->size);
                    newLeaf->points = mergedPoints;  // Move the merged points to the new leaf node
                    newLeaf->convertToLeaf();
                    newLeafNodes.push_back(newLeaf);    // Add node
                    mergedPoints.clear();
                }
                mergedPoints.insert(mergedPoints.end(), node->children[index]->points.begin(), node->children[index]->points.end());
                delete node->children[index]; 
                node->children[index] = nullptr;
            }

            // Handle remaining merged points
            if (!mergedPoints.empty()) {
                OctreeNode* newLeaf = new OctreeNode(node->origin, node->size);
                newLeaf->points = mergedPoints;
                newLeaf->convertToLeaf();
                newLeafNodes.push_back(newLeaf);
            }

            // Reinsert any remaining children after the merged one
            for (size_t i = 0; i < newLeafNodes.size(); i++) {
                node->children[i] = newLeafNodes[i];
            }

            // Clear any remaining child pointers
            for (size_t i = newLeafNodes.size(); i < 8; ++i) {
                node->children[i] = nullptr;
            }
        }
    }

}