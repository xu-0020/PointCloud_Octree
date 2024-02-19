#include <iostream>
#include <fstream>
#include <sstream>

#include "Octree.h"
using namespace std;

void buildOctreeFromCSV(const std::string& filename, Octree& octree) {
    ifstream file(filename);
    string line;
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

    Point origin = {0, 0, 0};       
    float initialSize = 100.0f;      
    int maxDepth = 100;
    int maxPointsPerNode = 100;

    Octree octree(origin, initialSize, maxDepth, maxPointsPerNode);

    string filename = "Montreal-PointCloud/Montreal1.csv";
    buildOctreeFromCSV(filename, octree);

    return 0;
}
