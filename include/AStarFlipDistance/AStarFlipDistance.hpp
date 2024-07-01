
#ifndef A_STAR_FOR_FLIPDISTANCE_ASTARFLIPDISTANCE_HPP
#define A_STAR_FOR_FLIPDISTANCE_ASTARFLIPDISTANCE_HPP

#include "AStarFlipDistance/GeometricObjects.hpp"
#include "AStarFlipDistance/StaticTriangulation.hpp"
#include "AStarFlipDistance/TriangulationHandler.hpp"
#include "ListHandler.hpp"
#include <chrono>
class AStarFlipDistance {
public:
    AStarFlipDistance(std::vector<Point_2> vertices,  std::vector<Triangle> start_triangulation,  std::vector<Triangle> target_triangulation );
    void run_combined();
    void run_root_combined();
    void run_epptein();
    void run_simple();
    int get_flipdistance();
    int get_runtime();
    int get_initial_heuristic() const;

    void constrain_the_instance_to_a_polygon(const std::vector<int>& polygon){
        _triangulation._objects.restrict_flips_to_polygon(polygon);
        _target_triangulation._objects.restrict_flips_to_polygon(polygon);
        _triangulation._objects.compute_distance_matrix();
    }


    std::string get_data_as_string();

    ListHandler _list;

    int _timeout_in_minutes=5;

    int counter_simple=0;
    int counter_eppstein=0;

    TriangulationHandler _triangulation;

    int get_nr_of_closed_nodes(){
        return _list.get_nr_of_extended_nodes();
    }

    int get_nr_of_opened_nodes(){
        return _list.get_nr_of_opened_nodes();
    }
private:



    static bool are_triangulation_representations_equal(const std::vector<triangle_int>& t_1, const std::vector<triangle_int>&t_2) ;


    StaticTriangulation _target_triangulation;
    int _flip_distance=-224;
    bool _has_run=false;
    int _runtime=-224;
    int _initial_heuristic=-22;

};


#endif //A_STAR_FOR_FLIPDISTANCE_ASTARFLIPDISTANCE_HPP
