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
        // Subdivide the node if it exceeds capacity and depth limit
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