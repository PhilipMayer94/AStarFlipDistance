
#ifndef A_STAR_FOR_FLIPDISTANCE_HEURISTICDISTANCECALCULATOR_HPP
#define A_STAR_FOR_FLIPDISTANCE_HEURISTICDISTANCECALCULATOR_HPP


#include <vector>
#include "BasicDataStructures.hpp"
#include "StaticTriangulation.hpp"
#include "TriangulationHandler.hpp"

class HeuristicDistanceCalculator {
public:
    HeuristicDistanceCalculator(const std::vector<Point_2>&  vertices , const std::vector<Triangle>& triangulation_1, std::vector<Triangle> triangulation_2);
    HeuristicDistanceCalculator(std::vector<Point_2>  vertices , std::vector<Triangle> triangulation_1, std::vector<Triangle> triangulation_2, GeometricObjects objects);
    void run();

    int run_and_return_only_hanke();

    void set_new_target_triangulation(std::vector<Triangle> triangles){
        StaticTriangulation tmp(_objects,triangles);
        _target=tmp;
    }

    void set_new_start_triangulation(std::vector<Triangle> triangles){
        _start.init_triangulation(triangles);
    }


    int simple_heuristic=-1;
    int eppstein_heuristic=-1;
    int hanke_heuristic=-1;

private:
    GeometricObjects _objects;
    std::vector<Point_2> _vertices;
    StaticTriangulation _target;
    TriangulationHandler _start;


};


#endif //A_STAR_FOR_FLIPDISTANCE_HEURISTICDISTANCECALCULATOR_HPP
