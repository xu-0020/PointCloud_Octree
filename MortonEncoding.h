#include <math.h>
#include <future>

#include "Bound.h"

const size_t max_concurrent_tasks = 30;     // concurrency control


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

vector<pair<int, uint64_t>> computeMortonCodesMultithreaded(const vector<Point>& dataPoints, const Bounds bounds) {
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