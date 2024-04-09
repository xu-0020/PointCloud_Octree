#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include "KdTree.cpp"
using namespace std;
#include "Octree_kd.cpp"
#include <experimental/filesystem>


// Octree functions
Bounds computeBoundingBoxFromCSV(const vector<string>& filenames) {
    Bounds bounds;
    bool firstPoint = true;

    for (const auto& filename : filenames) {
        ifstream file(filename);
        string line;

        // Skip the header
        if (!getline(file, line)) {
            cerr << "Error reading the file or the file is empty: " << filename << endl;
            continue;   // Skip to the next file
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

            if (firstPoint) {
                bounds.min = Point(x, y, z);
                bounds.max = Point(x, y, z);
                firstPoint = false;
            }
            else {
                bounds.update(Point(x, y, z));
            }
        }
    }
    return bounds;
}

void buildOctreeFromCSV(const vector<string>& filenames, Octree_kd& octree) {
    for (const auto& filename : filenames) {
        ifstream file(filename);
        string line;

        // Skip header
        getline(file, line);
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
            int r = stoi(values[3]);
            int g = stoi(values[4]);
            int b = stoi(values[5]);
            string source = values[6];
            Point point(x, y, z, r, g, b, source);
            octree.insert(point);
        }
    }
}
// End of Octree functions


int main() {
    vector<string> filenames;
    filenames.push_back("/export/project/hjingaa/PointCloud_Octree/data_csv/Montreal/TileO01_2013.csv");



    // Compute the bounding box 
    Bounds bounds = computeBoundingBoxFromCSV(filenames);
    // bounds.min.print_point();
    // bounds.max.print_point();
    Point origin = bounds.getCenter();       
    float initialSize = bounds.getSize();   


    // Normal Process
    /*
    
    // Variable settings
    int maxDepth = 7;
    int maxPointsPerNode = 10000;
    int minPointsPerNode = 10000;
    Octree octree(bounds, maxDepth, maxPointsPerNode, minPointsPerNode);

    buildOctreeFromCSV(filenames, octree);

    // Trim octree
    octree.trim(maxDepth/2);

    octree.buildRtrees();

    // octree.visualize("Octree Structure");

    // Range query
    Bounds queryRange;
    queryRange.min = Point(bounds.getCenter().x - 15, bounds.getCenter().y - 15, bounds.getCenter().z - 15); 
    queryRange.max = Point(bounds.getCenter().x + 15, bounds.getCenter().y + 15, bounds.getCenter().z + 15);  
    vector<Point> queryResults;

    octree.executeRangeQuery(queryRange, queryResults);

    cout << "Query Result: " << endl;
    for (const Point& point : queryResults) {
        cout << "Point: " << point.x << " " << point.y << " " << point.z << endl;
    }

    */



    // Code for testing the time for query and construrction

    // Variable settings arrays
    vector<int> multimMaxDepths = {3, 5, 7, 10};
    vector<int> multiMaxPointsPerNodes = {10000, 11000, 12000, 13000};
    vector<int> multiMinPointsPerNodes = {6000, 7000, 8000, 9000};

    // Open CSV file for output in append mode
    ofstream csvFile("octree_timing_results.csv", ios_base::app); 
    csvFile << "MaxDepth,MaxPointsPerNode,MinPointsPerNode,ConstructionTime(ms),RangeQueryTime(ms)\n";

    int fixed = 0;

    for (int i=0; i<multimMaxDepths.size(); i++) {
        Octree_kd octree(bounds, multimMaxDepths[i], multiMaxPointsPerNodes[fixed], multiMinPointsPerNodes[fixed]);

        // Measure construction time
        auto start1 = chrono::high_resolution_clock::now();
        buildOctreeFromCSV(filenames, octree);
        octree.trim(multimMaxDepths[i]/ 2);
        octree.buildKdtrees();
        auto stop1 = chrono::high_resolution_clock::now();

        // Measure range query time
        Bounds queryRange;
        vector<Point> queryResults;
        int searchSize = 50;
        queryRange.min = Point(bounds.getCenter().x - searchSize, bounds.getCenter().y - searchSize, bounds.getCenter().z - searchSize); 
        queryRange.max = Point(bounds.getCenter().x + searchSize, bounds.getCenter().y + searchSize, bounds.getCenter().z + searchSize);  

        auto start2 = chrono::high_resolution_clock::now();
        octree.executeRangeQuery(queryRange, queryResults);
        auto stop2 = chrono::high_resolution_clock::now();
        // Calculate durations
        auto constructionDuration = chrono::duration_cast<chrono::milliseconds>(stop1 - start1);
        auto queryDuration = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);

        // Write to CSV
        csvFile << multimMaxDepths[i] << "," << multiMaxPointsPerNodes[fixed] << "," << multiMaxPointsPerNodes[fixed] << ","
                << constructionDuration.count() << "," << queryDuration.count() << "\n";

    }
    
    csvFile.close();





    return 0;
}

