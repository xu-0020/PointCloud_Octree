#pragma once
#include <string>
#include <iostream>
#include <cstdint>


using namespace std;

// Define a point structure for standardized data points
struct Point {
    float x, y, z;              // Coordinates
    int r, g, b;                // RGB color data combined into a single integer
    string source;              // Record source file
    float arr[3] = {0.,0.,0.};
    float arrT[4] = {0.,0.,0.,0.};
    float time;
    Point() : x(0), y(0), z(0), r(-1), g(-1), b(-1), source(""), time(-1) {}

    Point(float x, float y, float z, int r = -1, int g = -1, int b = -1, string source = "", float time = -1) : x(x), y(y), z(z), r(r), g(g), b(b), time(time) {
        this->source = source;
        
        this->arr[0] = x;
        this->arr[1] = y;
        this->arr[2] = z;
        this->arrT[0] = x;
        this->arrT[1] = y;
        this->arrT[2] = z;
        this->arrT[3] = time;
    }

    void print_point(){
        if (time==-1){
            cout << "x = " <<x << " y = " << y << " z = " << z << " r = " << r << " g = " << g << " b = " << b << endl;
        }
        else{
            cout << "x = " <<x << " y = " << y << " z = " << z << " time = " << time << endl;
        }
        
    }

    // Morton Enconding functions
    uint64_t mortonCode(float x, float y, float z, float minX, float minY, float minZ, float maxX, float maxY, float maxZ) const {
        uint32_t normX = normalizeCoordinate(x, minX, maxX);
        uint32_t normY = normalizeCoordinate(y, minY, maxY);
        uint32_t normZ = normalizeCoordinate(z, minZ, maxZ);
        return interleaveBits(normX, normY, normZ);
    }

private:

    // Normalize the coordinate to [0, 2^32)
    uint32_t normalizeCoordinate(float coord, float minCoord, float maxCoord) const {
        const uint32_t maxRange = 0xFFFFFFFF; // 2^32 - 1
        return static_cast<uint32_t>((coord - minCoord) / (maxCoord - minCoord) * maxRange);    // Min-max normalization
    }

    uint64_t interleaveBits(uint32_t x, uint32_t y, uint32_t z) const {
        uint64_t result = 0;
        for (uint64_t i = 0; i < 21; i++) {     // Limit to 21 bits to avoid overflow
            result |= ((x & ((uint32_t)1 << i)) << (2 * i)) |
                    ((y & ((uint32_t)1 << i)) << (2 * i + 1)) |
                    ((z & ((uint32_t)1 << i)) << (2 * i + 2));
        }
        return result;
    }


};

