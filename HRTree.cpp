#include "RTreePoints.h"

class RTreeNode {
public: 
    float time;
    RTreePoints* tree;
    RTreeNode* next;
    RTreeNode(Point point, RTreeNode* next=nullptr): time(point.time), next(next){
        this->tree = new RTreePoints;
        RInsertPoint(this->tree, point);

        
    }
    ~RTreeNode(){
        delete this->tree;
        this->tree = nullptr;
        delete this->next;
        this->next = nullptr;
    }
    void HRTreeNodeInsert(Point point){
        RInsertPoint(this->tree, point);
        cout<<"buggggggggg"<<endl;
    }
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
class HRTree {
private:
    RTreeNode* head;
    RTreeNode* findPointTree(float time){
        //TO DO
        if (this->head == nullptr) {
            return nullptr;  // 链表为空
        }
        RTreeNode* current = this->head;
        while (current != nullptr && current->time != time) {
            current = current->next;
        }
        return current;
    }
    RTreeNode* findLastNode() {
        if (this->head == nullptr) {
            return nullptr;  // 链表为空
        }
        RTreeNode* current = this->head;
        while (current->next != nullptr) {
            current = current->next;
        }
        return current;
    }

public: 

    HRTree(): head(nullptr){}
    HRTree(vector<Point> points){
        if (points.empty()){
            head = nullptr;
            return;
        }
        for (Point point : points){
            if (this->head==nullptr){
                this->head = new RTreeNode(point);
            }
            else{
                RTreeNode* tree = findPointTree(point.time);
                if (tree) {tree->HRTreeNodeInsert(point);}
                else{
                    RTreeNode* current = findLastNode();
                    current->next = new RTreeNode(point);
                }
            }
            
        }
    }
    void HRTreeInsert(Point point){
        if (this->head==nullptr){
            this->head = new RTreeNode(point);
        }
        else{
            RTreeNode* tree = findPointTree(point.time);
            if (tree) tree->HRTreeNodeInsert(point);
            else{
                RTreeNode* current = findLastNode();
                current->next = new RTreeNode(point);
            }
        }
        return;
    }
    ~HRTree() {
        RTreeNode* current = this->head;
        while (current != nullptr) {
            RTreeNode* temp = current;
            current = current->next;
            delete temp;
        }
    }
    void print_tree(){
        if (this->head == nullptr) {
            cout<<"empty tree"<<endl;
            return;  // 链表为空
        }
        RTreeNode* current = this->head;
        while (current != nullptr) {
            current->print_tree_node();
            current = current->next;
        }
        return;
    }
    void HRTreeSearch(float time, std::vector<Point>& results, Bounds& queryRange){
        RTreeNode* tree = findPointTree(time);
        if(!tree) return;
        tree->HRSearch(results, queryRange);
        return;
    }
};


int main(){
    // points.push_back(tmp);
    // tmp = Point(1.2,2.5,3.1,1.2);
    // points.push_back(tmp);
    // tmp = Point(1.3,2.2,3.3,1.3);
    // points.push_back(tmp);
    // // for (Point point : points){
    // //     cout<< point.time;
    // // }
    // HRTree* tree = new HRTree(points);
    // // tree->print_tree();
    // std::vector<Point> results;
    // Point min = Point(0,0,0);
    // Point max = Point(15,15,15);
    // Bounds bound = Bounds(min, max);
    // tree->HRTreeSearch(1.2, results, bound);
    // for (Point point : results) point.print_point();




    return 0;
}