#ifndef A_STAR_FOR_FLIPDISTANCE_OUTERBICONNECTEDCOMPONENTFINDER_HPP
#define A_STAR_FOR_FLIPDISTANCE_OUTERBICONNECTEDCOMPONENTFINDER_HPP


#include <iostream>
#include <vector>
#include <stack>
#include "BasicDataStructures.hpp"
#include <float.h>

struct Node2D {
    Point_2 _point;
    std::vector<int> neighbors;
    int dfsNum;
    int low;
    bool isArticulation;

    Node2D(Point_2 p) : _point(p), dfsNum(-1), low(-1), isArticulation(false) {}
};

class OuterBiConnectedComponentFinder {
    std::vector<Node2D> nodes;
    std::vector<Edge> _edges;
    int dfsCounter;
    std::stack<int> dfsStack;
    std::vector<std::vector<int>> twoConnectedComponents;

public:
    OuterBiConnectedComponentFinder(const std::vector<Point_2>& points, const std::vector<Edge>& edges);

    void findTwoConnectedComponents();

    std::vector<int> find_outer_component();

    std::vector<Edge> find_outer_component_as_edges();


private:

    void dfsTwoConnected(int u, int parent);

    void markTwoConnectedComponent(int u, int v);

    void printTwoConnectedComponents() const;


};

#endif