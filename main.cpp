#include <iostream>
#include <fstream>
#include <sstream>

#include "Bound.h"
#include "Octree.h"
using namespace std;



Bounds computeBoundingBoxFromCSV(const string& filename) {
    ifstream file(filename);
    string line;
    Bounds bounds;

    // Skip the header
    if (!getline(file, line)) {
        cerr << "Error reading the file or the file is empty" << endl;
        exit(0);
    }

    bool firstPoint = true;
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
    return bounds;
}

void buildOctreeFromCSV(const string& filename, Octree& octree) {
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

int main() {
    string filename = "Montreal-PointCloud/Montreal1.csv";
    // Compute the bounding box from the CSV file
    Bounds bounds = computeBoundingBoxFromCSV(filename);

    Point origin = bounds.getCenter();       
    float initialSize = bounds.getSize();   

    int maxDepth = 100;
    int maxPointsPerNode = 1000;

    Octree octree(origin, initialSize, maxDepth, maxPointsPerNode);
    buildOctreeFromCSV(filename, octree);

    return 0;
}
