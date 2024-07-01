
#ifndef A_STAR_FOR_FLIPDISTANCE_STATICTRIANGULATION_HPP
#define A_STAR_FOR_FLIPDISTANCE_STATICTRIANGULATION_HPP

#include <utility>

#include "vector"
#include "AStarFlipDistance/GeometricObjects.hpp"

class StaticTriangulation {
public:
    StaticTriangulation()=default;
    StaticTriangulation(GeometricObjects objects, std::vector<Triangle>& triangles);
    int n;
    GeometricObjects _objects;
    std::vector<Point_2> _vertices;
    std::vector<std::vector<int>> _edge_matrix;
    std::vector<std::vector<std::pair<int,int>>> _current_flips_matrix;
    std::vector<int> _interior_edge_representation;
    std::vector<Triangle> _triangle_representation;
    std::vector<triangle_int> _short_representation;

    std::vector<bool> _edges_index_in_triangulation;
};




#endif //A_STAR_FOR_FLIPDISTANCE_STATICTRIANGULATION_HPP
