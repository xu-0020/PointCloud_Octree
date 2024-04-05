#include <vector>
#include <algorithm>
#include <cstdint>

using namespace std;


// Function to get the maximum value in the array of Morton codes
uint64_t getMax(const vector<pair<int, uint64_t>>& mortonCode) {
    uint64_t maxVal = 0;
    for (const auto& elem : mortonCode) {
        if (elem.second > maxVal) maxVal = elem.second;
    }
    return maxVal;
}

// Function to count sort according to the digit represented by exp.
void countSort(vector<pair<int, uint64_t>>& mortonCode, uint64_t exp) {
    vector<pair<int, uint64_t>> output(mortonCode.size());
    int count[10] = {0};

    // Store count of occurrences in count[]
    for (int i = 0; i < mortonCode.size(); i++) {
        count[(mortonCode[i].second / exp) % 10]++;
    }

    // Change count[i] so that count[i] now contains the actual position of this digit in output[]
    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }

    // Build the output array
    for (int i = mortonCode.size() - 1; i >= 0; i--) {
        output[count[(mortonCode[i].second / exp) % 10] - 1] = mortonCode[i];
        count[(mortonCode[i].second / exp) % 10]--;
    }

    // Copy the output, now contains sorted numbers according to the current digit
    for (int i = 0; i < mortonCode.size(); i++) {
        mortonCode[i] = output[i];
    }
}

// Radix Sort
void radixSort(vector<pair<int, uint64_t>>& mortonCode) {
    // Find the maximum number to know the number of digits
    uint64_t m = getMax(mortonCode);

    // Do counting sort for every digit. Note that instead of passing digit number, exp is passed. exp is 10^i where i is the current digit number
    for (uint64_t exp = 1; m / exp > 0; exp *= 10) {
        countSort(mortonCode, exp);
    }
}
