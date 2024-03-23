#include "RTreePoints.h"
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
using namespace std;
#include <filesystem>
namespace fs = filesystem;


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
    // Read in data.
    cout << "----------------------------------Reading data---------------------------------" << endl; 
    vector<Point> points = GetDataFromSingleCSV("/home/hjingaa/github/PointCloud_Octree/data_csv/whampoa_1223.csv");
    
    // Get the summary data.
    cout << "----------------------------------Process data summary---------------------------------" << endl; 
    std::unordered_map<std::string, float> map = Summary(points);
    // Print the summary data.
    cout<< "maxTime = " << map["maxTime"] << " " << "minTime = " << map["minTime"] << " " << "maxX = " << map["maxX"] << " " << "minX = " << map["minX"] << " " << "maxY = " << map["maxY"] << " " << "minY = " << map["minY"] << " " << "maxZ = " << map["maxZ"] << " " << "minZ = " << map["minZ"] << " " << endl;
    // compute the fixed bound.
    Point Cmin = Point(5.0 / 8 * map["minX"] + 3.0 / 8 * map["maxX"], 5.0 / 8 * map["minY"] + 3.0 / 8 * map["minY"], 5.0 / 8 * map["minZ"] + 3.0 / 8 * map["minZ"]);
    Point Cmax = Point(3.0 / 8 * map["minX"] + 5.0 / 8 * map["maxX"], 3.0 / 8 * map["minY"] + 5.0 / 8 * map["minY"], 3.0 / 8 * map["minZ"] + 5.0 / 8 * map["minZ"]);
    Bounds fixedBound = Bounds(Cmin, Cmax);
    // Compute the fixed time query
    float Tmin = 5.0 / 8 * map["minTime"] + 3.0 / 8 * map["maxTime"];
    float Tmax = 3.0 / 8 * map["minTime"] + 5.0 / 8 * map["maxTime"];


    // vector<Point> points; 
    // Point tmp = Point(1.1, 2, 3, 1.1, 1.1, 1.1, "", 1.1);
    // points.push_back(tmp);
    // tmp = Point(1.2, 2.5, 3.1, 1.2, 1.2, 1.2, "", 1.2);
    // points.push_back(tmp);
    // tmp = Point(1.3,2.2,3.3,1.3 ,1.3, 1.3, "", 1.3);
    // points.push_back(tmp);
    // // for (Point point : points){
    // //     cout<< point.time;
    // // }



    // Coodinate test.
    HRTree* tree = new HRTree(points);
    cout << "----------------------------------Begin the experiments---------------------------------" << endl; 
    std::map<string, float> time_cost_map;
    for (float i = 1.0; i < 100; i++){
        std::vector<Point> results;
        
        // Compoute the bound in this loop.
        float tmp_minX = (100.0+i) / 200.0 * map["minX"] + (100.0-i) / 200.0 * map["maxX"];
        float tmp_maxX = (100.0-i) / 200.0 * map["minX"] + (100.0+i) / 200.0 * map["maxX"];
        float tmp_minY = (100.0+i) / 200.0 * map["minY"] + (100.0-i) / 200.0 * map["maxY"];
        float tmp_maxY = (100.0-i) / 200.0 * map["minY"] + (100.0+i) / 200.0 * map["maxY"];
        float tmp_minZ = (100.0+i) / 200.0 * map["minZ"] + (100.0-i) / 200.0 * map["maxZ"];
        float tmp_maxZ = (100.0-i) / 200.0 * map["minZ"] + (100.0+i) / 200.0 * map["maxZ"];
        Point pmin = Point(tmp_minX, tmp_minY, tmp_minZ);
        Point pmax = Point(tmp_maxX, tmp_maxY, tmp_maxZ);
        Bounds bound = Bounds(pmin, pmax);
        cout << "The bound is (" << tmp_minX << ", " << tmp_maxX << "), (" << tmp_minY << ", " << tmp_maxY << "), (" << tmp_minZ << ", " << tmp_maxZ << ") " << endl;
        auto startTime = std::chrono::high_resolution_clock::now();
        tree->HRTreeSearch(Tmin, Tmax, results, bound);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        // Make the bounds words
        string bound_words = "The bound is (" + to_string(tmp_minX) + ", " + to_string(tmp_maxX) + "), (" + to_string(tmp_minY) + ", " + to_string(tmp_maxY) + "), (" + to_string(tmp_minZ) + ", " + to_string(tmp_maxZ) + ") with time bound (" + to_string(Tmin) + ", " + to_string(Tmax) + ").";
        time_cost_map[bound_words] = duration.count();
    }
    


    // Time query test.



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