#include "RTree.h"
#include "Point.h" 
#include <vector>
#include <iostream>
typedef RTree<Point, float, 4> RTreeT;
typedef Point ValueType;


void RTInsert(RTreeT* tree, vector<Point> points){
    if (points.size()==0){
        return;
    }
    for (Point point : points){
        // Point* pointer = new Point();
        tree->Insert(point.arrT, point.arrT, point);
    }
}


void RTSearch(RTreeT* tree, std::vector<Point>& results, Bounds& queryRange){
    
    tree->Search(queryRange.min.arrT, queryRange.max.arrT, results, MySearchCallback);
    return;
}

