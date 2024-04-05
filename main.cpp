#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>

#include "Octree.h"
#include "RadixSort.h"
using namespace std;

#include <filesystem>
namespace fs = filesystem;

const size_t max_concurrent_tasks = 30;     // concurrency control



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




// Single file processing
void processFileAndInsertPoints(const string& filename, vector<Point>& dataPoints, mutex& dataMutex) {
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
        //int r = stoi(values[3]);
        //int g = stoi(values[4]);
        //int b = stoi(values[5]);
        //string source = values[6];
        //Point point(x, y, z, r, g, b, source);

        Point point(x, y, z);

        // Synchronize access
        lock_guard<mutex> guard(dataMutex);
        dataPoints.push_back(point);
    }
}

void readFromCSV(const vector<string>& filenames, vector<Point>& dataPoints) {
    mutex dataMutex;
    vector<future<void>> futures;

    for (const auto& filename : filenames) {
        // Check if reached the maximum number of concurrent tasks
        if (futures.size() >= max_concurrent_tasks) {
            // Wait for at least one task to complete
            bool taskCompleted = false;
            while (!taskCompleted) {
                for (auto it = futures.begin(); it != futures.end(); ) {
                    auto& fut = *it;
                    if (fut.wait_for(chrono::seconds(0)) == future_status::ready) {
                        fut.get();      // Get the result to clear any stored exception
                        it = futures.erase(it); // Remove the completed future
                        taskCompleted = true;
                        break;      // Break the loop as we only need one task to complete
                    } else {
                        it++;
                    }
                }
            }
        }

        futures.push_back(async(launch::async, processFileAndInsertPoints, filename, ref(dataPoints), ref(dataMutex)));
    }

    for (auto& fut : futures) {
        fut.get();      // waits for the tasks to finish
    }
}

// Morton encoding Functions
vector<pair<int, uint64_t>> computeMortonCodesSegment(const vector<Point>& dataPoints, const Bounds& bounds, size_t start, size_t end) {
    vector<pair<int, uint64_t>> segmentMortonCodes;
    for (size_t i = start; i < end; i++) {
        const Point& point = dataPoints[i];
        uint64_t code = point.mortonCode(point.x, point.y, point.z, bounds.min.x, bounds.min.y, bounds.min.z, bounds.max.x, bounds.max.y, bounds.max.z);
        segmentMortonCodes.push_back(make_pair(i, code));
    }
    return segmentMortonCodes;
}

vector<pair<int, uint64_t>> computeMortonCodesMultithreaded(const vector<Point>& dataPoints, const Bounds& bounds) {
    const size_t segmentSize = ceil(dataPoints.size() / static_cast<double>(max_concurrent_tasks));

    vector<future<vector<pair<int, uint64_t>>>> futures;

    for (size_t i = 0; i < max_concurrent_tasks && i * segmentSize < dataPoints.size(); i++) {
        size_t start = i * segmentSize;
        size_t end = min(start + segmentSize, dataPoints.size());

        futures.push_back(async(launch::async, [&dataPoints, &bounds, start, end] {
            return computeMortonCodesSegment(dataPoints, bounds, start, end);
        }));
    }

    vector<pair<int, uint64_t>> mortonCode;
    for (auto& fut : futures) {
        auto segmentMortonCodes = fut.get();
        mortonCode.insert(mortonCode.end(), segmentMortonCodes.begin(), segmentMortonCodes.end());
    }

    return mortonCode;
}




int main() {
    string choice;
    cout << "Please select input file by folder name(0) or path(1):" << endl;
    getline(cin, choice);

    vector<string> filenames;   // Pre-allocate for optimization

    if (choice == "0") {
        string folderPath;
        cout << "Enter the path to the folder containing CSV files:" << endl;
        getline(cin, folderPath);

        for (auto& entry : fs::directory_iterator(folderPath)) {
            if (entry.path().extension() == ".csv") {
                filenames.push_back(entry.path().string());
            }
        }
    }

    else if (choice == "1") {
        string filePath;
        while (true) {
            cout << "Enter CSV file paths (Enter 'done' when finished):" << endl;
            getline(cin, filePath);
            if (filePath == "done") { 
                break;
            }
            filenames.push_back(filePath);
        }
    }


    if (filenames.empty()) {
        cerr << "No filenames provided." << endl;
        return 1; 
    }

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



    // Read from CSV data
    vector<Point> dataPoints;
    readFromCSV(filenames, dataPoints);
    

    auto start10 = chrono::high_resolution_clock::now();

    // Morton encoding, corresponding to each dataPoints, (index, mortoncode)
    vector<pair<int, uint64_t>> mortonCode = computeMortonCodesMultithreaded(dataPoints, bounds);
    
    auto stop10 = chrono::high_resolution_clock::now();

    // Sort based on Morton code
    radixSort(mortonCode);


    // Tree Variable Settings
    int maxDepth = 10;
    int maxPointsPerLeaf = 6000;


    Octree octree(bounds, maxDepth, maxPointsPerLeaf);

    auto start1 = chrono::high_resolution_clock::now();
    octree.constructOctree(dataPoints, mortonCode);
    auto stop1 = chrono::high_resolution_clock::now();
    
    //octree.visualize("Test 1");

    auto start2 = chrono::high_resolution_clock::now();
    octree.rebalance(dataPoints);
    octree.trim(dataPoints);
    auto stop2 = chrono::high_resolution_clock::now();

    octree.visualize("Test 2");
    
    auto start3 = chrono::high_resolution_clock::now();
    octree.buildRtrees(dataPoints);
    auto stop3 = chrono::high_resolution_clock::now();

    // Measure range query time
    Bounds queryRange;
    vector<Point> queryResults;
    int searchSize = 50;
    queryRange.min = Point(bounds.getCenter().x - searchSize, bounds.getCenter().y - searchSize, bounds.getCenter().z - searchSize); 
    queryRange.max = Point(bounds.getCenter().x + searchSize, bounds.getCenter().y + searchSize, bounds.getCenter().z + searchSize);  

    auto start4 = chrono::high_resolution_clock::now();
    octree.executeRangeQuery(queryRange, queryResults);
    auto stop4 = chrono::high_resolution_clock::now();

    

    // Calculate durations
    auto constructionDuration = chrono::duration_cast<chrono::milliseconds>(stop1 - start1);
    auto rebalanceDuration = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
    auto buildRDuration = chrono::duration_cast<chrono::milliseconds>(stop3 - start3);
    auto queryDuration = chrono::duration_cast<chrono::milliseconds>(stop4 - start4);

    auto mortonDuration = chrono::duration_cast<chrono::milliseconds>(stop10 - start10);



    // Open CSV file for output in append mode
    ofstream csvFile("octree_timing_results_bottom.csv", ios_base::app); 
    csvFile << "ConstructionTime(ms),RebalanceTime(ms),RtreeTime(ms),RangeQueryTime(ms),MortonEncodingTime(ms)\n";
    csvFile << constructionDuration.count() << "," << rebalanceDuration.count() << "," << buildRDuration.count() << "," << queryDuration.count() << "," << mortonDuration.count() << "\n";
    
    //csvFile << constructionDuration.count() << "," << rebalanceDuration.count() << "," << mortonDuration.count() << "\n";


    csvFile.close();

    
    return 0;
}

