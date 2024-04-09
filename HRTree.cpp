#include "RTreePoints.h"
#include "RTreeT.h"
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <list>
using namespace std;
#include <unordered_map>
#include <random>
#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include <experimental/filesystem>

void write_result(string filename, list<string> results){
    std::ofstream outputFile(filename, std::ios::trunc);
    if (outputFile.is_open()) { // Check if the file is successfully opened
        for (const auto& str : results) {
            outputFile << str << std::endl; // Write each string to the file
        }

        outputFile.close(); // Close the file
        std::cout << "File write successful." << std::endl;
    } else {
        std::cout << "Failed to open the file." << std::endl;
    }

    return;
}
// summary of the points(max and min of x, y, z, time)
std::unordered_map<std::string, float> Summary(vector<Point> points){
    std::unordered_map<std::string, float> summary;
    auto maxTPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.time < b.time;
        }
    );
    auto minTPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.time > b.time;
        }
    );
    auto maxXPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.time < b.time;
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
    summary["maxTime"] = maxTPoint->time;
    summary["minTime"] = minTPoint->time;
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



// summary of the points(max and min of x, y, z, time)
std::unordered_map<std::string, float> Summary_e(vector<Point> points){
    std::unordered_map<std::string, float> summary;
    auto maxXPoint = std::max_element(
        points.begin(),
        points.end(),
        [](Point& a, Point& b) {
            return a.time < b.time;
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


class RTreeNode {
public: 
    float time;
    RTreePoints* tree;
    RTreeNode(vector<Point> points): time(points.front().time) {
        this->tree = new RTreePoints;
        RInsert(this->tree, points);
    }
    ~RTreeNode(){
        delete this->tree;
        this->tree = nullptr;
    }
    // void HRTreeNodeInsert(Point point){
    //     RInsertPoint(this->tree, point);
    // }
    void HRSearch(std::vector<Point>& results, Bounds& queryRange){
        RSearch(this->tree, results, queryRange);
        return;
    }
    void print_tree_node(){
        cout<< "Tree node with the time "<< this->time<<endl;
        auto list = tree->ListTree();
        int counter = 0;
        for (auto aabb : list) {
        cout << "TreeList [" << counter++ << "]: "
            << aabb.m_min[0] << ", "
            << aabb.m_min[1] << ", "
            << aabb.m_min[2] << "; "
            << aabb.m_max[0] << ", "
            << aabb.m_max[1] << ", "
            << aabb.m_max[2] << endl;
        }
        return;
    }
};

bool compareByAttribute(RTreeNode& obj1, RTreeNode& obj2) {
    return obj1.time < obj2.time;
}
bool compareByPoints(Point& obj1, Point& obj2) {
    return obj1.time < obj2.time;
}

class HRTree {
private:
    std::map<float, RTreeNode*> trees;
public: 

    HRTree(vector<Point> points){
        if (points.empty()){
            return;
        }
        std::sort(points.begin(), points.end(), compareByPoints);
         // Group the elements based on 'attribute'
        std::map<float, std::vector<Point>> groupedMap;
        for (const Point& point : points) {
            groupedMap[point.time].push_back(point);
        }
        for (const auto& entry : groupedMap) {
            // cout<< entry.first;
            this->trees[entry.first] = new RTreeNode(groupedMap[entry.first]);
        }
        return;
    }
    ~HRTree() {
        this->trees.clear();
    }
    void print_tree(){
        for (auto& entry : this->trees) {
            entry.second->print_tree_node();
        }
        return;
    }
    void HRTreeSearch(float start, float end, std::vector<Point>& results, Bounds& queryRange){
        for (auto& entry : this->trees) {
            if (entry.first >= start && entry.first <= end){
                entry.second->HRSearch(results, queryRange);
            }
        }
        return;
    }
};


int main(){
    


    // // Read in data.
    // auto startTime = std::chrono::high_resolution_clock::now();
    // cout << "----------------------------------Reading data---------------------------------" << endl; 
    // string folderPath = "/export/project/hjingaa/PointCloud_Octree/data_csv/Montreal"; // 替换为实际的文件夹路径
    // vector<Point> points = GetDataFromSingleFolder(folderPath);
    // auto endTime = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Reading data cost                        "  <<duration.count()<< "ms"<<endl;
    // // Get the summary data.
    // cout << "----------------------------------Process data summary---------------------------------" << endl; 
    // std::unordered_map<std::string, float> map = Summary_e(points);
    // // Print the summary data.
    // cout << "maxX = " << map["maxX"] << " " << "minX = " << map["minX"] << " " << "maxY = " << map["maxY"] << " " << "minY = " << map["minY"] << " " << "maxZ = " << map["maxZ"] << " " << "minZ = " << map["minZ"] << " " << endl;
    // Bounds maxBound = Bounds(Point(map["minX"], map["minY"], map["minZ"]), Point(map["maxX"], map["maxY"], map["maxZ"]));

    // // RTree test
    // cout << "----------------------------------Build the Rtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // RTreePoints* Rtree = new RTreePoints;
    // RInsert(Rtree, points);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Buiding the Rtree cost                        "  <<duration.count()<< "ms"<<endl;

    // // Octree+Rtree Coordinate test.
    // cout << "----------------------------------Begin the experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost;
    // for (float i = 1.0; i < 100.0; i++){
    //     std::vector<Point> results;
    //     // Get the random points.
    //     std::vector<Point> query_points = getRandomPoints(points, i * points.size() / 10000.0);
    //     std::unordered_map<std::string, float> query_map = Summary(query_points);

    //     // Compoute the bound in this loop.
    //     float tmp_minX = query_map["minX"];
    //     float tmp_maxX = query_map["maxX"];
    //     float tmp_minY = query_map["minY"];
    //     float tmp_maxY = query_map["maxY"];
    //     float tmp_minZ = query_map["minZ"];
    //     float tmp_maxZ = query_map["maxZ"];
    //     Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
    //     Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
    //     Bounds bound = Bounds(pmin, pmax);
    //     auto startTime = std::chrono::high_resolution_clock::now();
    //     RSearch(Rtree, results, bound);
    //     auto endTime = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    //     // Make the bounds words
    //     string bound_words = "The bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ")." + " Time cost = " + to_string(duration.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     result_time_cost.push_back(bound_words);
    //     cout << bound_words << endl;
    // }
    // write_result("result_montreal/Result_R.txt", result_time_cost);




    // // Read in data.
    // auto startTime = std::chrono::high_resolution_clock::now();
    // cout << "----------------------------------Reading data---------------------------------" << endl; 
    // vector<Point> points = GetDataFromSingleCSV("data_csv/whampoa_0521.csv");
    // auto endTime = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Reading data cost                        "  <<duration.count()<< "ms"<<endl;
    // // Get the summary data.
    // cout << "----------------------------------Process data summary---------------------------------" << endl; 
    // std::unordered_map<std::string, float> map = Summary(points);
    // // Print the summary data.
    // cout<< "maxTime = " << map["maxTime"] << " " << "minTime = " << map["minTime"] << " " << "maxX = " << map["maxX"] << " " << "minX = " << map["minX"] << " " << "maxY = " << map["maxY"] << " " << "minY = " << map["minY"] << " " << "maxZ = " << map["maxZ"] << " " << "minZ = " << map["minZ"] << " " << endl;
    // // compute the fixed bound.
    // Point Cmin = Point(7./10. *map["minX"] + 3./10. * map["maxX"],7./10. *map["minY"] + 3./10. * map["maxY"],7./10. *map["minZ"] + 3./10. * map["maxZ"]);
    // Point Cmax = Point(3./10. *map["minX"] + 7./10. * map["maxX"],3./10. *map["minY"] + 7./10. * map["maxY"],3./10. *map["minZ"] + 7./10. * map["maxZ"]);
    // Bounds fixedBound = Bounds(Cmin, Cmax);
    // // Compute the fixed time query
    // float Tmin = 7./10. *map["minTime"] + 3./10. * map["maxTime"];
    // float Tmax = 3./10. *map["minTime"] + 7./10. * map["maxTime"];

    // // Build hrtree
    // cout << "----------------------------------Build the hrtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // HRTree* tree = new HRTree(points);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Buiding the hrtree cost                        "  <<duration.count()<< "ms"<<endl;

    // // HRTree Coordinate test.
    // cout << "----------------------------------Begin the Coordinate experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost;
    // for (float i = 1.0; i < 100; i++){
    //     std::vector<Point> results;
        
    //     // Compoute the bound in this loop.
    //     float tmp_minX = (143.0+i * 57.0/100.0) / 200.0 * map["minX"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxX"];
    //     float tmp_maxX = (57.0-i * 57.0/100.0) / 200.0 * map["minX"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxX"];
    //     float tmp_minY = (143.0+i * 57.0/100.0) / 200.0 * map["minY"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxY"];
    //     float tmp_maxY = (57.0-i * 57.0/100.0) / 200.0 * map["minY"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxY"];
    //     float tmp_minZ = (143.0+i * 57.0/100.0) / 200.0 * map["minZ"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxZ"];
    //     float tmp_maxZ = (57.0-i * 57.0/100.0) / 200.0 * map["minZ"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxZ"];
    //     Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
    //     Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
    //     Bounds bound = Bounds(pmin, pmax);
    //     auto startTime = std::chrono::high_resolution_clock::now();
    //     tree->HRTreeSearch(Tmin, Tmax, results, bound);
    //     auto endTime = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    //     // Make the bounds words
    //     string bound_words = "The bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ") with time bound (" + to_string(Tmin) + ", " + to_string(Tmax) + ")." + " Time cost = " + to_string(duration.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     result_time_cost.push_back(bound_words);
    //     cout << bound_words << endl;
    // }
    // write_result("result/FtimeRVbound_HR.txt", result_time_cost);


    // // HRTree Time test 
    // cout << "----------------------------------Begin the Time experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost;
    // for (float i = 1.0; i < 100; i++){
    //     std::vector<Point> results;
    //     float tmp_minT = (100.0+i) / 200.0 * map["minTime"] + (100.0-i) / 200.0 * map["maxTime"];
    //     float tmp_maxT = (100.0-i) / 200.0 * map["minTime"] + (100.0+i) / 200.0 * map["maxTime"];
    //     // Compoute the bound in this loop.
    //     auto startTime = std::chrono::high_resolution_clock::now();
    //     tree->HRTreeSearch(tmp_minT, tmp_maxT, results, fixedBound);
    //     auto endTime = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    //     // Make the bounds words
    //     string bound_words = "The bound is (" + to_string(fixedBound.min.x) + ", " + to_string(fixedBound.max.x) + "), (" + to_string(fixedBound.min.y) + ", " + to_string(fixedBound.max.x) + "), (" + to_string(fixedBound.min.z) + ", " + to_string(fixedBound.max.x) + ") with time bound (" + to_string(tmp_minT) + ", " + to_string(tmp_maxT) + ")." + " Time cost = " + to_string(duration.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     result_time_cost.push_back(bound_words);
    //     cout << bound_words << endl;
    // }
    // write_result("result/FboundRVtime_HR.txt", result_time_cost);

    // // RTree test
    // cout << "----------------------------------Build the Rtree---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // RTreePoints* Rtree = new RTreePoints;
    // RInsert(Rtree, points);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Buiding the Rtree cost                        "  <<duration.count()<< "ms"<<endl;
    // // Coordinate test.
    // cout << "----------------------------------Begin the Coordinate experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost;
    // for (float i = 1.0; i < 100; i++){
    //     std::vector<Point> results;
        
    //     // Compoute the bound in this loop.
    //     float tmp_minX = (143.0+i * 57.0/100.0) / 200.0 * map["minX"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxX"];
    //     float tmp_maxX = (57.0-i * 57.0/100.0) / 200.0 * map["minX"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxX"];
    //     float tmp_minY = (143.0+i * 57.0/100.0) / 200.0 * map["minY"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxY"];
    //     float tmp_maxY = (57.0-i * 57.0/100.0) / 200.0 * map["minY"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxY"];
    //     float tmp_minZ = (143.0+i * 57.0/100.0) / 200.0 * map["minZ"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxZ"];
    //     float tmp_maxZ = (57.0-i * 57.0/100.0) / 200.0 * map["minZ"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxZ"];
    //     Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
    //     Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
    //     Bounds bound = Bounds(pmin, pmax);
    //     auto startTime = std::chrono::high_resolution_clock::now();
    //     RSearch(Rtree, results, bound);
    //     auto endTime = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    //     // Make the bounds words
    //     string bound_words = "The bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ") Time cost = " + to_string(duration.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     result_time_cost.push_back(bound_words);
    //     cout << bound_words << endl;
    // }
    // write_result("result/RVbound_R.txt", result_time_cost);

    // // RTreeT Coordinate test
    // cout << "----------------------------------Build the RtreeT---------------------------------" << endl; 
    // startTime = std::chrono::high_resolution_clock::now();
    // RTreeT* RtreeT = new RTreeT;
    // RTInsert(RtreeT, points);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    // cout<< "Buiding the RtreeT cost                        "  <<duration.count()<< "ms"<<endl;


    // // Coordinate test.
    // cout << "----------------------------------Begin the Coordinate experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost;
    // for (float i = 1.0; i < 100; i++){
    //     std::vector<Point> results;
        
    //     // Compoute the bound in this loop.
    //     float tmp_minX = (143.0+i * 57.0/100.0) / 200.0 * map["minX"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxX"];
    //     float tmp_maxX = (57.0-i * 57.0/100.0) / 200.0 * map["minX"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxX"];
    //     float tmp_minY = (143.0+i * 57.0/100.0) / 200.0 * map["minY"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxY"];
    //     float tmp_maxY = (57.0-i * 57.0/100.0) / 200.0 * map["minY"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxY"];
    //     float tmp_minZ = (143.0+i * 57.0/100.0) / 200.0 * map["minZ"] + (57.0-i * 57.0/100.0) / 200.0 * map["maxZ"];
    //     float tmp_maxZ = (57.0-i * 57.0/100.0) / 200.0 * map["minZ"] + (143.0+i * 57.0/100.0) / 200.0 * map["maxZ"];
    //     Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ, -1, -1, -1, "", Tmin);
    //     Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ, -1, -1, -1, "", Tmax);
    //     Bounds bound = Bounds(pmin, pmax);
    //     auto startTime = std::chrono::high_resolution_clock::now();
    //     RTSearch(RtreeT, results, bound);
    //     auto endTime = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    //     // Make the bounds words
    //     string bound_words = "The bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ") with time bound (" + to_string(Tmin) + ", " + to_string(Tmax) + ")." + " Time cost = " + to_string(duration.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     result_time_cost.push_back(bound_words);
    //     cout << bound_words << endl;
    // }
    // write_result("result/FtimeRVbound_RT.txt", result_time_cost);

    // // Time test.
    // // HRTree Time test 
    // cout << "----------------------------------Begin the Time experiments---------------------------------" << endl; 
    // std::list<string> result_time_cost;
    // for (float i = 1.0; i < 100; i++){
    //     std::vector<Point> results;
    //     float tmp_minT = (100.0+i) / 200.0 * map["minTime"] + (100.0-i) / 200.0 * map["maxTime"];
    //     float tmp_maxT = (100.0-i) / 200.0 * map["minTime"] + (100.0+i) / 200.0 * map["maxTime"];
    //     // Compoute the bound in this loop.
    //     Point pmin = Point(map["minX"],map["minY"],map["minZ"], -1, -1, -1, "", tmp_minT);
    //     Point pmax = Point(map["maxX"],map["maxY"],map["maxZ"], -1, -1, -1, "", tmp_maxT);
    //     Bounds bound = Bounds(pmin, pmax);
    //     auto startTime = std::chrono::high_resolution_clock::now();
    //     RTSearch(RtreeT, results, bound);
    //     auto endTime = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    //     // Make the bounds words
    //     string bound_words = "The bound is (" + to_string(fixedBound.min.x) + ", " + to_string(fixedBound.max.x) + "), (" + to_string(fixedBound.min.y) + ", " + to_string(fixedBound.max.x) + "), (" + to_string(fixedBound.min.z) + ", " + to_string(fixedBound.max.x) + ") with time bound (" + to_string(tmp_minT) + ", " + to_string(tmp_maxT) + ")." + " Time cost = " + to_string(duration.count())+" microseconds." + "find " + to_string(results.size()) + " points.";
    //     result_time_cost.push_back(bound_words);
    //     cout << bound_words << endl;
    // }
    // write_result("result/FboundRVtime_RT.txt", result_time_cost);











    // // Some test.
    // HRTree* tree = new HRTree(points);
    // // // tree->print_tree();
    // std::vector<Point> results;
    // Point min = Point(0,0,0);
    // Point max = Point(15,15,15);
    // Bounds bound = Bounds(min, max);
    // tree->HRTreeSearch(1.2, 1.2, results, bound);
    // for (Point point : results) point.print_point();




    return 0;
}