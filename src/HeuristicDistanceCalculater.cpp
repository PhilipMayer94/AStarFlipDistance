
#include <utility>

#include "AStarFlipDistance/HeuristicDistanceCalculator.hpp"
#include "AStarFlipDistance/FastHankeHeuristic.hpp"
#include "AStarFlipDistance/FastEppsteinHeuristic.hpp"

HeuristicDistanceCalculator::HeuristicDistanceCalculator(const std::vector<Point_2>& vertices,
                                                         const std::vector<Triangle>& triangulation_1,
                                                         std::vector<Triangle> triangulation_2) :_objects(GeometricObjects(vertices,false)), _start(vertices,_objects),
                                                                                                  _target(_objects,triangulation_2){
    _start.init_triangulation(triangulation_1);
    _vertices=vertices;
}

void HeuristicDistanceCalculator::run() {
    auto copy_triang=_start.get_triangle_representation();


    FastHankeHeuristic hanke(_objects,_target,_start.get_edge_representation(),_start._current_flips_matrix);
    int hanke_value=hanke.compute_heuristic();

    FastEppsteinHeuristic eppstein(_vertices,  _start.get_triangle_representation(), _target._triangle_representation,_objects);
    eppstein.run();
    int eppstein_value=eppstein._heuristic_value;

    simple_heuristic=_start.get_simple_heuristic(_target);
    eppstein_heuristic=eppstein._heuristic_value;
    hanke_heuristic=hanke_value;

    //std::cout<<std::endl<<"hanke: "<<hanke_heuristic<<" eppstein: "<<eppstein_heuristic<<" simple: "<<simple_heuristic<<std::endl;
    _start.init_triangulation(copy_triang);
}


HeuristicDistanceCalculator::HeuristicDistanceCalculator(std::vector<Point_2> vertices,
                                                         std::vector<Triangle> triangulation_1,
                                                         std::vector<Triangle> triangulation_2,
                                                         GeometricObjects objects):_objects(std::move(objects)), _start(vertices,_objects),
                                                                                   _target(_objects,triangulation_2){
    _start.init_triangulation(triangulation_1);
    _vertices=vertices;

}

int HeuristicDistanceCalculator::run_and_return_only_hanke() {
    FastHankeHeuristic hanke(_objects,_target,_start.get_edge_representation(),_start._current_flips_matrix);
    int hanke_value=hanke.compute_heuristic();

    return hanke_value;
}
