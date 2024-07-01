
#ifndef A_STAR_FOR_FLIPDISTANCE_BIDIRECTIONALBFS_HPP
#define A_STAR_FOR_FLIPDISTANCE_BIDIRECTIONALBFS_HPP

#include <queue>
#include "AStarFlipDistance/ListHandler.hpp"
#include "AStarFlipDistance/TriangulationHandler.hpp"
#include "AStarFlipDistance/StaticTriangulation.hpp"

class BiDirectionalBfs {
private:
    std::unordered_map<std::vector<triangle_int>,int,Hasher,EqualFn> _closed_from_start;
    std::unordered_map<std::vector<triangle_int>,int,Hasher,EqualFn> _closed_from_target;
    std::queue<std::pair<std::vector<triangle_int>,int>> _queue_from_start;
    std::queue<std::pair<std::vector<triangle_int>,int>> _queue_from_target;

    TriangulationHandler _triangulation_from_start;
    TriangulationHandler _triangulation_from_target;
    StaticTriangulation _start;
    StaticTriangulation _target;
    GeometricObjects _objects;


public:
    BiDirectionalBfs(std::vector<Point_2> points, std::vector<Triangle> start_triangulation, std::vector<Triangle> target_triangulation);
    void run();
    int _flipdistance;
    int _runtime;
    int _closed;
    std::string get_data_as_string();
    int _timeout_in_minutes=5;
};


#endif //A_STAR_FOR_FLIPDISTANCE_BIDIRECTIONALBFS_HPP
