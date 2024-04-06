#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>

#include "Octree.h"
#include "MortonEncoding.h"
#include "RadixSort.h"
#include "CreateBox.h"
using namespace std;

#include <filesystem>
namespace fs = filesystem;



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
    

    const Bounds bounds = computeTotalBoundingBox(filenames);   // Compute the bounding box of the entire dataset

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


    Octree octree(maxDepth, maxPointsPerLeaf);

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

