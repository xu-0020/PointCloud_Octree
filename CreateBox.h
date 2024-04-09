#include <fstream>
#include <sstream>
#include <future>

#include "Bound.h"

Bounds computeBoundingBoxFromSingleCSV(const string& filename) {
    Bounds bounds;
    bool firstPoint = true;

    ifstream file(filename);
    string line;

    if (!getline(file, line)) { // Skip header
        cerr << "Error reading the file or the file is empty: " << filename << endl;
        return bounds; // Return empty bounds if file is empty or can't be read
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
        } else {
            bounds.update(Point(x, y, z));
        }
    }
    return bounds;
}


Bounds computeTotalBoundingBox(const vector<string>& filenames) {
    // Acquiring the bounding box of entire datasets.
    vector<future<Bounds>> futureBounds;
    // Launch a task for each file
    for (const auto& filename : filenames) {
        futureBounds.push_back(async(launch::async, computeBoundingBoxFromSingleCSV, filename));
    }
    Bounds bounds;
    bool firstBounds = true;

    // Wait for all tasks to complete and merge the results
    for (auto& fb : futureBounds) {
        Bounds b = fb.get();    // Block until the future is ready
        if (firstBounds) {
            bounds = b;
            firstBounds = false;
        } else {
            bounds.update(b.min);
            bounds.update(b.max);
        }
    }
    // Boundingbox acquired.
    return bounds;
}