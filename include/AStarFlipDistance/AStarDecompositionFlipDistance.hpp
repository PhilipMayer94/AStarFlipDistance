

#ifndef A_STAR_FOR_FLIPDISTANCE_ASTARDECOMPOSITIONFLIPDISTANCE_HPP
#define A_STAR_FOR_FLIPDISTANCE_ASTARDECOMPOSITIONFLIPDISTANCE_HPP


#include "BasicDataStructures.hpp"
#include "Decomposer.hpp"
#include "AStarFlipDistance.hpp"
#include "OuterBiConnectedComponentFinder.hpp"
#include "FastHankeHeuristic.hpp"

class AStarDecompositionFlipDistance {
private:
    std::vector<Polygon> _polygons;
    std::vector<Triangle> _start_triangulation;
    std::vector<Triangle> _target_triangulation;
    std::vector<Point_2> _vertices;

    std::vector<Edge> _fixed_edges;
    std::vector<Edge> _outer_component;
    bool _have_edges_already_been_fixed= false;


    void compute_fixed_edges();
    void compute_decomposition();

    bool _used_hanke=false;
    int _runtime=0;
    int _flip_distance=-1;


public:

    std::string get_data_as_string() {
        std::string result=std::to_string(get_flipdistance())+" "+std::to_string((double)get_runtime()/1000)+" "+ std::to_string(_used_hanke);

        return result;
    }


    std::vector<AStarFlipDistance> _sub_problems;

    AStarDecompositionFlipDistance(std::vector<Point_2> &vertices,  std::vector<Triangle> &start_triangulation,  std::vector<Triangle> &target_triangulation ):
    _start_triangulation(start_triangulation),_target_triangulation(target_triangulation),_vertices(vertices){}

    void set_fixed_edges(std::vector<Edge> &edges);

    void run_with_combined();

    int get_flipdistance(){

        return _flip_distance;
    }
    int get_runtime(){
        return _runtime;

        int result=0;
        for(auto astar:_sub_problems){
            if(astar.get_runtime()>0){
                result+=astar.get_runtime();
            }
        }
        return result;
    }

    int get_initial_heuristic(){
        int result=0;
        for(auto astar:_sub_problems){
            if(astar.get_initial_heuristic()>0){
                result+=astar.get_initial_heuristic();
            }

        }
        return result;
    }

    int get_extended(){
        int result=0;
        for(auto astar:_sub_problems){
            result+=astar._list.get_nr_of_extended_nodes();
        }
        return result;
    }

    int get_open(){
        int result=0;
        for(auto astar:_sub_problems){
            result+=astar._list.get_nr_of_opened_nodes();
        }
        return result;
    }








};


#endif //A_STAR_FOR_FLIPDISTANCE_ASTARDECOMPOSITIONFLIPDISTANCE_HPP
