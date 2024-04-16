#include "Octree.h"
#include "MortonEncoding.h"
#include "RadixSort.h"
#include "CreateBox.h"
#include <list>
#include <filesystem>
#include <unordered_map>
#include <iostream>
#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include <experimental/filesystem>
#include <random>
#include "KdTree.cpp"
using namespace std;
// #include "Octree_kd.cpp"
namespace fs = std::experimental::filesystem;
void write_result(string filename, list<string> results){
    std::ofstream outputFile(filename, std::ios::trunc);
    if (outputFile.is_open()) { // Check if the file is successfully opened
        for (const auto& str : results) {
            outputFile << str << std::endl; // Write each string to the file
        }

        outputFile.close(); // Close the file
        std::cout << filename << " File write successful." << std::endl;
    } else {
        std::cout << "Failed to open the file." << std::endl;
    }

    return;
}
// summary of the points(max and min of x, y, z, time)
std::unordered_map<std::string, float> Summary(vector<Point> points){
    std::unordered_map<std::string, float> summary;
    auto maxXPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.x < b.x;
        }
    );
    auto minXPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.x > b.x;
        }
    );
    auto maxYPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.x < b.x;
        }
    );
    auto minYPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.y > b.y;
        }
    );
    auto maxZPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.z < b.z;
        }
    );
    auto minZPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.z > b.z;
        }
    );
    summary["maxX"] = maxXPoint->x;
    summary["minX"] = minXPoint->x;
    summary["maxY"] = maxYPoint->y;
    summary["minY"] = minYPoint->y;
    summary["maxZ"] = maxZPoint->z;
    summary["minZ"] = minZPoint->z;
    return summary;


}

std::vector<Point> getRandomPoints(const std::vector<Point>& data, int count) {
    std::vector<Point> randomPoints;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, data.size() - 1);

    for (int i = 0; i < count; ++i) {
        int index = dis(gen);
        randomPoints.push_back(data[index]);
    }

    return randomPoints;
}

vector<Point> GetDataFromSingleCSV(const string& filename){
    bool firstPoint = true;
    vector<Point> points;
    ifstream file(filename);
    string line;
    if (!getline(file, line)) { // Skip header
        cerr << "Error reading the file or the file is empty: " << filename << endl;
        return points; // Return empty bounds if file is empty or can't be read
    }

    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        vector<string> values;

        while (getline(ss, value, ',')) {
            values.push_back(value);
        }

        float x = stof(values[0]);
        float y = stof(values[1]);
        float z = stof(values[2]);
        float time = stof(values[3]);
        Point point = Point(x, y, z, -1, -1, -1, filename, floor(time));
        points.push_back(point);
    }
 
    return points;


}


vector<Point> GetDataFromSingleFolder(const std::string& folderPath) {
    vector<Point> points;
    vector<Point> tmp;
    DIR* directory = opendir(folderPath.c_str());
    if (directory == nullptr) {
        std::cout << "Failed to open directory: " << folderPath << std::endl;
        return points;
    }

    dirent* entry;
    while ((entry = readdir(directory)) != nullptr) {
        std::string fileName = entry->d_name;
        if (fileName != "." && fileName != "..") {
            std::string absolutePath = folderPath + "/" + fileName;
            std::cout << "File: " << absolutePath << std::endl;
            tmp = GetDataFromSingleCSV(absolutePath);
            points.insert(points.end(), tmp.begin(), tmp.end());
        }
    }

    closedir(directory);
    return points;
}

// Build the octree + rtree
Octree buildTree(vector<Point> dataPoints,const Bounds bounds, int maxPointsPerLeaf) {
    // const Bounds bounds = computeTotalBoundingBox(filenames);   // Compute the bounding box of the entire dataset

    // Morton encoding, corresponding to each dataPoints, (index, mortoncode)
    vector<pair<int, uint64_t>> mortonCode = computeMortonCodesMultithreaded(dataPoints, bounds);
    
    // Sort based on Morton code
    radixSort(mortonCode);

    // Tree Variable Settings
    int maxDepth = 10;

    Octree octree(maxDepth, maxPointsPerLeaf);
    octree.constructOctree(dataPoints, mortonCode);
    
    octree.rebalance(dataPoints);
    octree.trim(dataPoints);
    
    octree.buildRtrees(dataPoints);
    return octree;
}


Octree buildTree_kd(vector<Point> dataPoints,const Bounds bounds, int maxPointsPerLeaf) {
    // const Bounds bounds = computeTotalBoundingBox(filenames);   // Compute the bounding box of the entire dataset

    // Morton encoding, corresponding to each dataPoints, (index, mortoncode)
    vector<pair<int, uint64_t>> mortonCode = computeMortonCodesMultithreaded(dataPoints, bounds);
    
    // Sort based on Morton code
    radixSort(mortonCode);

    // Tree Variable Settings
    int maxDepth = 10;

    Octree octree(maxDepth, maxPointsPerLeaf);
    octree.constructOctree(dataPoints, mortonCode);
    
    octree.rebalance(dataPoints);
    octree.trim(dataPoints);
    
    octree.buildKdtrees(dataPoints);
    return octree;
}




void runQuery(Bounds& queryRange, vector<Point>& results, Octree& octree) {
    octree.executeRangeQuery(queryRange, results);
}


void runQuery_kd(Bounds& queryRange, vector<Point>& results, Octree& octree) {
    octree.executeRangeQuery(queryRange, results);
}
//************************************************
//code from Octree.cpp
//************************************************
// R-tree functions
void RInsert(RTreePoints* tree, const vector<Point>& points){
    if (points.size()==0){
        return;
    }
    for (const Point& point : points){
        // Point* pointer = new Point();
        tree->Insert(point.arr, point.arr, point);
    }
}

bool MySearchCallback(ValueType point)
{
    return true; // keep going
}

void RSearch(RTreePoints* tree, std::vector<Point>& results, Bounds& queryRange){
    auto list = tree->ListTree();
    int counter = 0;
    tree->Search(queryRange.min.arr, queryRange.max.arr, results, MySearchCallback);
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



void KdInsert(KdTree* tree, vector<Point>& points){
    if (points.size() == 0){
        return;
    }
    tree->buildTree(points);
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






void Octree::buildFromLeafNodes(vector<OctreeNode*>& nodes) {
    if (nodes.size() <= 1) {
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
    uint32_t parentCount = ceil(static_cast<double>(nodes.size()) / 8.0);  // Calculate the number of parent nodes needed

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
    buildFromLeafNodes(parentNodes);
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

void KdSearch(KdTree* kdtree, std::vector<Point>& results, Bounds& queryRange){
    kdtree->rangeQuery(results,queryRange);
    return;
}

void Octree::rangeQuery(Bounds& queryRange, vector<Point>& results, OctreeNode* node) {

    if (node->isLeaf()) {
        if (node->rtree)
            RSearch(node->rtree, results, queryRange);
        else if (node->kdtree){
            vector<Point> tmp;
            KdSearch(node->kdtree, tmp, queryRange);
            results.insert(results.end(), tmp.begin(), tmp.end());
        }
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
void Octree::initializeRTrees(OctreeNode* node, vector<future<void>>& futures, const vector<Point>& dataPoints, RTreeStats& stats) {
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
        futures.push_back(async(launch::async, [this, node, &dataPoints, &stats]() {
            // Regenerate bounds for the leaf node to ensure tight fitting
            node->bound = calculateBoundsForPoints(dataPoints, node->points);

            node->rtree = new RTreePoints();
            vector<Point> insertPoints;         // Retrieve points to be inserted in to the Rtree.
            for (int i=0; i<node->points.size(); i++) {
                insertPoints.push_back(dataPoints[node->points[i]]);
            }

            // Update RTreeStats
            {
                lock_guard<mutex> guard(stats.mtx);
                stats.totalLeafNodes += 1;
                stats.totalPoints += insertPoints.size();
                stats.minPoints = min(stats.minPoints, insertPoints.size());
                stats.maxPoints = max(stats.maxPoints, insertPoints.size());
            }

            node->points.clear();
            RInsert(node->rtree, insertPoints);
        }));
    } else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                initializeRTrees(node->children[i], futures, dataPoints, stats);
            }
        }
    }
}



void Octree::initializeKdTrees(OctreeNode* node, vector<future<void>>& futures, const vector<Point>& dataPoints) {
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

            node->kdtree = new KdTree();
            vector<Point> insertPoints;         // Retrieve points to be inserted in to the Rtree.
            for (int i=0; i<node->points.size(); i++) {
                insertPoints.push_back(dataPoints[node->points[i]]);
            }
            node->points.clear();
            KdInsert(node->kdtree, insertPoints);
        }));
    } else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                initializeKdTrees(node->children[i], futures, dataPoints);
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


int main() {
    
    // Read in data.
    auto startTime = std::chrono::high_resolution_clock::now();
    cout << "----------------------------------Reading data---------------------------------" << endl; 
    string folderPath = "/export/project/hjingaa/PointCloud_Octree/data_csv/HK_large"; // 替换为实际的文件夹路径
    vector<Point> points = GetDataFromSingleFolder(folderPath);
    // vector<Point> points = GetDataFromSingleCSV("/export/project/hjingaa/PointCloud_Octree/data_csv/HK/demo.csv");
    cout<< "The number of points is " << points.size() << endl;
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    cout<< "Reading data cost                        "  <<duration.count()<< "ms"<<endl;
    // Get the summary data.
    cout << "----------------------------------Process data summary---------------------------------" << endl; 
    int count = points.size();
    std::unordered_map<std::string, float> map = Summary(points);
    // Print the summary data.
    cout << "maxX = " << map["maxX"] << " " << "minX = " << map["minX"] << " " << "maxY = " << map["maxY"] << " " << "minY = " << map["minY"] << " " << "maxZ = " << map["maxZ"] << " " << "minZ = " << map["minZ"] << " " << endl;
    Bounds maxBound = Bounds(Point(map["minX"], map["minY"], map["minZ"]), Point(map["maxX"], map["maxY"], map["maxZ"]));

    // // Build the Octree+Rtree
    // cout << "----------------------------------Build the Octree+Rtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // Octree tree = buildTree(points, maxBound, count / 100);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Building the Octree+Rtree cost            "  <<duration.count()<< "ms"<<endl;

    // // RTree test
    // cout << "----------------------------------Build the Rtree---------------------------------" << endl;
    // startTime = std::chrono::high_resolution_clock::now();
    // RTreePoints* Rtree = new RTreePoints;
    // RInsert(Rtree, points);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Buiding the Rtree cost                        "  <<duration.count()<< "ms"<<endl;

    // OcTree+KDTree test
    cout << "----------------------------------Build the Octree+KDtree---------------------------------" << endl; 
    startTime = std::chrono::high_resolution_clock::now();
    Octree oktree = buildTree_kd(points, maxBound, count/100);
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    cout<< "Buiding the Octree+KDtree cost                        "  <<duration.count()<< "ms"<<endl;

    points.clear();

    cout << "----------------------------------Begin the experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost_or;
    // std::list<string> result_time_cost_r;
    std::list<string> result_time_cost_ok;
    for (float i = 1.0; i <= 20.0; i++){
        cout << endl;
        cout << "**********Run the " << i << "th experiment.**********" << endl;
        std::vector<Point> results_or;
        std::vector<Point> results_r;
        std::vector<Point> results_ok;


        // Compoute the bound in this loop.
        float tmp_minX = (143.0+i * 57.0/20.0) / 200.0 * map["minX"] + (57.0-i * 57.0/20.0) / 200.0 * map["maxX"];
        float tmp_maxX = (57.0-i * 57.0/20.0) / 200.0 * map["minX"] + (143.0+i * 57.0/20.0) / 200.0 * map["maxX"];
        float tmp_minY = (143.0+i * 57.0/20.0) / 200.0 * map["minY"] + (57.0-i * 57.0/20.0) / 200.0 * map["maxY"];
        float tmp_maxY = (57.0-i * 57.0/20.0) / 200.0 * map["minY"] + (143.0+i * 57.0/20.0) / 200.0 * map["maxY"];
        float tmp_minZ = (143.0+i * 57.0/20.0) / 200.0 * map["minZ"] + (57.0-i * 57.0/20.0) / 200.0 * map["maxZ"];
        float tmp_maxZ = (57.0-i * 57.0/20.0) / 200.0 * map["minZ"] + (143.0+i * 57.0/20.0) / 200.0 * map["maxZ"];
        Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
        Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
        Bounds bound = Bounds(pmin, pmax);


        // // Get the random points.
        // std::vector<Point> query_points = getRandomPoints(points, i * points.size() / 10000.0);
        // std::unordered_map<std::string, float> query_map = Summary(query_points);
        // float tmp_minX = query_map["minX"];
        // float tmp_maxX = query_map["maxX"];
        // float tmp_minY = query_map["minY"];
        // float tmp_maxY = query_map["maxY"];
        // float tmp_minZ = query_map["minZ"];
        // float tmp_maxZ = query_map["maxZ"];
        // Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
        // Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
        // Bounds bound = Bounds(pmin, pmax);

        // // get a certain point
        // std::vector<Point> query_points = getRandomPoints(points, 1);
        // std::unordered_map<std::string, float> query_map = Summary(query_points);
        // float tmp_minX = query_map["minX"];
        // float tmp_maxX = query_map["maxX"];
        // float tmp_minY = query_map["minY"];
        // float tmp_maxY = query_map["maxY"];
        // float tmp_minZ = query_map["minZ"];
        // float tmp_maxZ = query_map["maxZ"];
        // Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
        // Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
        // Bounds bound = Bounds(pmin, pmax);


        // auto startTime_or = std::chrono::high_resolution_clock::now();
        // runQuery(bound, results_or, tree);
        // auto endTime_or = std::chrono::high_resolution_clock::now();
        // auto duration_or = std::chrono::duration_cast<std::chrono::microseconds>(endTime_or - startTime_or);
        // // Make the bounds words
        // string bound_words_or = "ORORORThe bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration_or.count())+" microseconds." + "find " + to_string(results_or.size()) + " points.";
        // result_time_cost_or.push_back(bound_words_or);
        // cout << bound_words_or << endl;


        // auto startTime_r = std::chrono::high_resolution_clock::now();
        // RSearch(Rtree, results_r, bound);
        // auto endTime_r = std::chrono::high_resolution_clock::now();
        // auto duration_r = std::chrono::duration_cast<std::chrono::microseconds>(endTime_r - startTime_r);
        // // Make the bounds words
        // string bound_words_r = "RRRThe bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration_r.count())+" microseconds." + "find " + to_string(results_r.size()) + " points.";
        // result_time_cost_r.push_back(bound_words_r);
        // cout << bound_words_r << endl;


        auto startTime_ok = std::chrono::high_resolution_clock::now();
        runQuery_kd(bound, results_ok, oktree);
        auto endTime_ok = std::chrono::high_resolution_clock::now();
        auto duration_ok = std::chrono::duration_cast<std::chrono::microseconds>(endTime_ok - startTime_ok);
        // Make the bounds words
        string bound_words_ok = "OKOKOKThe bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration_ok.count())+" microseconds." + "find " + to_string(results_ok.size()) + " points.";
        result_time_cost_ok.push_back(bound_words_ok);
        cout << bound_words_ok << endl;
        cout << endl;
    }
    string file_path = "result_hkl/";
    // write_result(file_path + "Result_OR.txt", result_time_cost_or);
    // write_result(file_path + "Result_R.txt", result_time_cost_r);
    write_result(file_path + "Result_OK.txt", result_time_cost_ok);



    // // Build the Octree+Rtree
    // cout << "----------------------------------Build the Octree+Rtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // Octree tree4000 = buildTree(points, maxBound, 4000);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Building the Octree+Rtree cost            "  <<duration.count()<< "ms"<<endl;

    // // Build the Octree+Rtree
    // cout << "----------------------------------Build the Octree+Rtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // Octree tree6000 = buildTree(points, maxBound, 6000);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Building the Octree+Rtree cost            "  <<duration.count()<< "ms"<<endl;

    // // Build the Octree+Rtree
    // cout << "----------------------------------Build the Octree+Rtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // Octree tree8000 = buildTree(points, maxBound, 8000);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Building the Octree+Rtree cost            "  <<duration.count()<< "ms"<<endl;

    // // Build the Octree+Rtree
    // cout << "----------------------------------Build the Octree+Rtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // Octree tree10000 = buildTree(points, maxBound, 10000);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Building the Octree+Rtree cost            "  <<duration.count()<< "ms"<<endl;

    // points.clear();

    // // Octree+Rtree Coordinate test.
    // cout << "----------------------------------Begin the experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost_ok;
    // for (float i = 1.0; i < 2.0; i++){
    //     cout << endl;
    //     cout << "**********Run the " << i << "th experiment.**********" << endl;
    //     std::vector<Point> results;


    //     // Compoute the bound in this loop.
    //     float tmp_minX = (143.0+i * 57.0/50.0) / 200.0 * map["minX"] + (57.0-i * 57.0/50.0) / 200.0 * map["maxX"];
    //     float tmp_maxX = (57.0-i * 57.0/50.0) / 200.0 * map["minX"] + (143.0+i * 57.0/50.0) / 200.0 * map["maxX"];
    //     float tmp_minY = (143.0+i * 57.0/50.0) / 200.0 * map["minY"] + (57.0-i * 57.0/50.0) / 200.0 * map["maxY"];
    //     float tmp_maxY = (57.0-i * 57.0/50.0) / 200.0 * map["minY"] + (143.0+i * 57.0/50.0) / 200.0 * map["maxY"];
    //     float tmp_minZ = (143.0+i * 57.0/50.0) / 200.0 * map["minZ"] + (57.0-i * 57.0/50.0) / 200.0 * map["maxZ"];
    //     float tmp_maxZ = (57.0-i * 57.0/50.0) / 200.0 * map["minZ"] + (143.0+i * 57.0/50.0) / 200.0 * map["maxZ"];
    //     Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
    //     Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
    //     Bounds bound = Bounds(pmin, pmax);


    //     auto startTime_ok4000 = std::chrono::high_resolution_clock::now();
    //     runQuery(bound, results, tree4000);
    //     auto endTime_ok4000 = std::chrono::high_resolution_clock::now();
    //     auto duration_ok4000 = std::chrono::duration_cast<std::chrono::microseconds>(endTime_ok4000 - startTime_ok4000);
    //     // Make the bounds words
    //     string bound_words_ok4000 = "4000ORORORThe bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration_ok4000.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     cout << bound_words_ok4000 << endl;
    //     cout << endl;

    //     auto startTime_ok6000 = std::chrono::high_resolution_clock::now();
    //     runQuery(bound, results, tree6000);
    //     auto endTime_ok6000 = std::chrono::high_resolution_clock::now();
    //     auto duration_ok6000 = std::chrono::duration_cast<std::chrono::microseconds>(endTime_ok6000 - startTime_ok6000);
    //     // Make the bounds words
    //     string bound_words_ok6000 = "6000ORORORThe bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration_ok6000.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     cout << bound_words_ok6000 << endl;
    //     cout << endl;

    //     auto startTime_ok8000 = std::chrono::high_resolution_clock::now();
    //     runQuery(bound, results, tree8000);
    //     auto endTime_ok8000 = std::chrono::high_resolution_clock::now();
    //     auto duration_ok8000 = std::chrono::duration_cast<std::chrono::microseconds>(endTime_ok8000 - startTime_ok8000);
    //     // Make the bounds words
    //     string bound_words_ok8000 = "8000ORORORThe bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration_ok8000.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     cout << bound_words_ok8000 << endl;
    //     cout << endl;

    //     auto startTime_ok10000 = std::chrono::high_resolution_clock::now();
    //     runQuery(bound, results, tree10000);
    //     auto endTime_ok10000 = std::chrono::high_resolution_clock::now();
    //     auto duration_ok10000 = std::chrono::duration_cast<std::chrono::microseconds>(endTime_ok10000 - startTime_ok10000);
    //     // Make the bounds words
    //     string bound_words_ok10000 = "10000ORORORThe bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration_ok10000.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     cout << bound_words_ok10000 << endl;
    //     cout << endl;
    // }

    return 0;
}