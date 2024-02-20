#pragma once
#include <string>
using namespace std;

// Define a point structure for standardized data points
struct Point {
    float x, y, z;              // Coordinates
    int r, g, b;                // RGB color data combined into a single integer
    string source;              // Record source file


    Point() : x(0), y(0), z(0), r(-1), g(-1), b(-1), source("") {}

    Point(float x, float y, float z, int r = -1, int g = -1, int b = -1, string source = "") : x(x), y(y), z(z), r(r), g(g), b(b) {
        this->source = source;
    }
};