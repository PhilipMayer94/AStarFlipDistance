

#ifndef A_STAR_FOR_FLIPDISTANCE_FASTHANKEHEURISTIC_HPP
#define A_STAR_FOR_FLIPDISTANCE_FASTHANKEHEURISTIC_HPP

#include <vector>
#include "StaticTriangulation.hpp"
#include "SimpleMaxHeap.hpp"
#include "AStarFlipDistance/TriangulationHandler.hpp"

class FastHankeHeuristic {
public:

    FastHankeHeuristic(std::vector<Point_2> &vertices,  const std::vector<Triangle>& start_triangulation, std::vector<Triangle>& target_triangulation );

    FastHankeHeuristic(GeometricObjects & objects,StaticTriangulation & target, std::vector<int>  interior_diagonals , std::vector<std::vector<std::pair<int,int>>> flip_matrix);

    FastHankeHeuristic( std::vector<Point_2> &vertices, const std::vector<Triangle> &start_triangulation,
                       std::vector<Triangle> &target_triangulation,GeometricObjects& objects);


    void do_flip(int diagonal_index);

    int compute_heuristic();

    void set_new_start_triangulation( const std::vector<int>&  interior_diagonals , std::vector<std::vector<std::pair<int,int>>> &flip_matrix);

    int compute_heuristic_for_triangulation(std::vector<int> interior_diagonals , std::vector<std::vector<std::pair<int,int>>> &flip_matrix);

    void initalize_datastructures();





    int different_edge_counter;
    int _nr_of_diagonals;
    std::vector<int> _current_diagonals;
    std::vector<int> _current_diagonal_weights;

    std::vector<int> _number_of_intersections_with_target;

    std::vector<std::vector<int>> _where_in_diagonal_list;
    std::vector<std::vector<std::pair<int,int>>> _flip_matrix;

    int _flip_distance=0;

    MaxHeap _heap;


    std::vector<Point_2> _vertices;
    GeometricObjects _objects;
    StaticTriangulation _target;
};


#endif //A_STAR_FOR_FLIPDISTANCE_FASTHANKEHEURISTIC_HPP
