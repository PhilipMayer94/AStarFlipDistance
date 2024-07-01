

#ifndef A_STAR_FOR_FLIPDISTANCE_FASTEPPSTEINHEURISTIC_HPP
#define A_STAR_FOR_FLIPDISTANCE_FASTEPPSTEINHEURISTIC_HPP


#include <utility>
#include <vector>
#include <queue>
#include "BasicDataStructures.hpp"
#include "GeometricObjects.hpp"
#include "TriangulationHandler.hpp"
#include "StaticTriangulation.hpp"

class FastEppsteinHeuristic {
public:

    FastEppsteinHeuristic(std::vector<Point_2> &vertices,  const std::vector<Triangle>& start_triangulation, std::vector<Triangle>& target_triangulation,GeometricObjects & obj);

    FastEppsteinHeuristic(std::vector<Point_2> &vertices,  const std::vector<Triangle>& start_triangulation, std::vector<Triangle>& target_triangulation );

    void run();

    int _heuristic_value=-1;
    std::vector<Point_2> _vertices;
    TriangulationHandler _start_triangulation;
    StaticTriangulation _target_triangulation;
    GeometricObjects _objects;

};


#endif //A_STAR_FOR_FLIPDISTANCE_FASTEPPSTEINHEURISTIC_HPP
