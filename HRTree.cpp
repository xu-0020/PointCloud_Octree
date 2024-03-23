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

    vector<Point> points = GetDataFromSingleCSV("/home/hjingaa/github/PointCloud_Octree/data_csv/whampoa_1223.csv");
    std::unordered_map<std::string, float> map = Summary(points);
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
    HRTree* tree = new HRTree(points);
    // // tree->print_tree();
    std::vector<Point> results;
    Point min = Point(0,0,0);
    Point max = Point(15,15,15);
    Bounds bound = Bounds(min, max);
    tree->HRTreeSearch(1.2, 1.2, results, bound);
    for (Point point : results) point.print_point();




    return 0;
}