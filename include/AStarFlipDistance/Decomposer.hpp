
#ifndef A_STAR_FOR_FLIPDISTANCE_DECOMPOSER_HPP
#define A_STAR_FOR_FLIPDISTANCE_DECOMPOSER_HPP


#include "BasicDataStructures.hpp"
#include "Geometry.hpp"

struct Polygon{
    std::vector<int> polygon_walk{};
    std::vector<Point_2> all_points{};
    std::vector<Triangle> start_triangulation{};
    std::vector<Triangle> target_triangulation{};
};

struct Node {
    Point_2 _point;
    std::vector<int> neighbors;


    explicit Node(Point_2 p) : _point(p) {}
};

class Decomposer {
    std::vector<Node> nodes;
    std::vector<Edge> _edges;

    std::vector<std::vector<int>> polygons;
    std::vector<std::vector<int>> interior_points;
    std::vector<std::vector<Triangle>> start_triangulations;
    std::vector<std::vector<Triangle>> target_triangulations;
public:
    Decomposer(const std::vector<Point_2>& points,const std::vector<Edge>& edges);



    void compute_polygons();

    void compute_interior_points();

    void distribute_triangulations(const std::vector<Triangle>& start_triangulation, const std::vector <Triangle>& target_triangulation);

    std::vector<Polygon> compute_polygon_decomposition(const std::vector<Triangle>& start_triangulation, const std::vector <Triangle>& target_triangulation);


private:

    double shoelaceArea(const std::vector<int>& polygon);

    int get_next_neighbor(const Node& pred, Node current);

    static bool isLeft(const Point_2& a, const Point_2& b, const Point_2& c);

    static bool rayIntersectsEdge(const Point_2& point, const Point_2& edgeStart, const Point_2& edgeEnd);

    static bool nodeInPolygon(const Point_2& point, const std::vector<Node>& nodes, const std::vector<int>& polygon);

    static bool edge_in_node_polygon(int e_1,int e_2, const std::vector<Node>& nodes, const std::vector<int>& polygon);

    static bool nodeInPolygon(int index, const std::vector<Node>& nodes, const std::vector<int>& polygon);
};











#endif //A_STAR_FOR_FLIPDISTANCE_DECOMPOSER_HPP
