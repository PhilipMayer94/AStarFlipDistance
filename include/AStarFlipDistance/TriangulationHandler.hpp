

#ifndef A_STAR_FOR_FLIPDISTANCE_TRIANGULATIONHANDLER_HPP
#define A_STAR_FOR_FLIPDISTANCE_TRIANGULATIONHANDLER_HPP

#include <utility>
#include "vector"
#include "AStarFlipDistance/BasicDataStructures.hpp"
#include "AStarFlipDistance/Geometry.hpp"
#include "AStarFlipDistance/GeometricObjects.hpp"
#include "StaticTriangulation.hpp"
#include "AStarFlipDistance/HungarianPrimDual.hpp"


class TriangulationHandler {
public:
    TriangulationHandler()=default;
    TriangulationHandler(const std::vector<Point_2>& vertices, GeometricObjects objects);

    void init_triangulation(std::vector<triangle_int>& triangles);
    void init_triangulation(const std::vector<Triangle>& triangles);
    void reset_triangulation();
    void reset_relative_triangulation();


    int get__simple_heuristic_offset(Edge e, Edge f,StaticTriangulation  & target);

    std::vector<Triangle> do_n_random_flips(int nn, int seed);


    Edge do_random_flip();



    Edge do_flip_for_real(Edge e);
    Edge do_flip_implicit(Edge e);
    std::pair<std::vector<triangle_int>,int> do_flip_implicit_for_simple_heuristic(Edge e,const StaticTriangulation& target);

    int get_simple_heuristic(const StaticTriangulation & target_triangulation);
    bool equals_triangulation(const StaticTriangulation & target_triangulation);

    int get_eppstein_heuristic_naive(const StaticTriangulation & target_triangulation);
    std::pair<std::vector<int>,std::vector<int>> get_relevant_edges_for_eppstein(const StaticTriangulation & target_triangulation);
    std::vector<std::vector<triangle_int>> get_distance_matrix_for_eppstein(const StaticTriangulation & target_triangulation);

    int get_eppstein_heuristic_parent(const StaticTriangulation &target_triangulation);
    int get_eppstein_heuristic_child(const StaticTriangulation &target_triangulation,
                                 int flipped_parent_edge, int child_edge);


    std::vector<Triangle> do_n_valid_flips(int n, int seed);
    int get_heuristic_upper_bound(const StaticTriangulation & target);
    std::vector<int> get_intersection_difference_count(const StaticTriangulation & target);

    std::vector<triangle_int> get_representation();
    std::vector<int> get_edge_representation();
    std::vector<Triangle> get_triangle_representation();


    void print_stuff();
    GeometricObjects _objects;
    std::vector<std::vector<std::pair<int,int>>> _current_flips_matrix;

private:

    HungarianPrimDual _dynamic_hungarian_matrix;
    int n;
    std::vector<Point_2> _vertices;
    std::vector<std::vector<int>> _edge_matrix;
    std::vector<int> _interior_edge_representation;
    std::vector<Triangle> _triangle_representation;
    std::vector<triangle_int> _short_representation;

    int binarySearch(const std::vector<triangle_int>& vec, int val);
    void bubble_to_correct_spot(std::vector<triangle_int>& vec,int pos);

};


#endif //A_STAR_FOR_FLIPDISTANCE_TRIANGULATIONHANDLER_HPP
