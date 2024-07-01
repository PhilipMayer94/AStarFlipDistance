
#ifndef A_STAR_FOR_FLIPDISTANCE_GEOMETRICOBJECTS_HPP
#define A_STAR_FOR_FLIPDISTANCE_GEOMETRICOBJECTS_HPP

#include <utility>

#include "vector"
#include "AStarFlipDistance/BasicDataStructures.hpp"
#include "AStarFlipDistance/Geometry.hpp"
#include "Global.hpp"


class GeometricObjects {

public:
     GeometricObjects()=default;
     GeometricObjects(std::vector<Point_2> vertices,bool compute_distance=true);

    void print_stuff();

    Edge get_edge(int e_index);
    int get_edge(int e_0,int e_1);
    int get_number_of_vertices();
    int get_number_of_triangles();
    Triangle get_triangle(int i);
    int get_triangle(Triangle t);
    int get_triangle(int i,int j, int k);
    bool is_flip_valid(int e_0, int e_1,int f_0,int f_1);
    bool is_flip_valid(int edge_0, int edge_1);

    int get_number_of_edges();
    void compute_distance_matrix();

    int get_distance(int edge_ind_1,int edge_ind_2);

    std::vector<Point_2> get_vertices();
    void restrict_flips_to_order(int order);



    static bool isLeft(const Point_2& a, const Point_2& b, const Point_2& c){
        return ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) > 0;
    }

    static bool rayIntersectsEdge(const Point_2& point, const Point_2& edgeStart, const Point_2& edgeEnd){
        return (edgeStart.y > point.y) != (edgeEnd.y > point.y) &&
               point.x < (edgeEnd.x - edgeStart.x) * (point.y - edgeStart.y) / (edgeEnd.y - edgeStart.y) + edgeStart.x;
    }

    bool point_in_polygon(const Point_2& point, const std::vector<int>& polygon){
        bool inside = false;
        int n__ = polygon.size();
        for (int i = 0, j = n__ - 1; i < n__; j = i++) {
            Point_2 vi = _vertices[polygon[i]];
            Point_2 vj = _vertices[polygon[j]];

            if ((vi.y > point.y) != (vj.y > point.y) &&
                point.x < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y) + vi.x) {
                inside = !inside;
            }
        }

        return inside;
    }

    bool point_in_polygon(int index, const std::vector<int> &polygon) {
        if(std::find(polygon.begin(),polygon.end(),index)!=polygon.end()){
            return false;
        }
        return point_in_polygon(_vertices[index],polygon);
    }

    void restrict_flips_to_polygon(const std::vector<int>& polygon){
        std::vector<bool> is_inside_polygon(_edge_list.size(),true);
        for(int n=0;n<_edge_list.size();n++){
            auto e =_edge_list[n];
            int e_1=e.src; int e_2=e.dst;
            int index_1=-3; int index_2=-6;

            for(int i=0;i<polygon.size();i++){
                if(polygon[i]==e_1){
                    index_1=i;
                }
                if(polygon[i]==e_2){
                    index_2=i;
                }
            }
            if((index_1==0 &&index_2==polygon.size()-1)||(index_2==0 &&index_1==polygon.size()-1)){
                is_inside_polygon[n]=(true);
                continue;
            }
            if(abs(index_1-index_2)<=1){
                is_inside_polygon[n]=(true);
                continue;
            }

            bool intersects=false;
            for(int p=0;p<polygon.size();p++){
                int q=(p+1)%polygon.size();
                Point_2 point_1=_vertices[polygon[p]];
                Point_2 point_2=_vertices[polygon[q]];
                if(do_arcs_intersect(point_1,point_2,_vertices[e.src],_vertices[e.dst])){
                    intersects=true;
                }
            }
            if(intersects){
                //std::cout<<e.src<<"----------"<<e.dst<<std::endl;
                is_inside_polygon[n]=(false);
                continue;
            }

            Point_2 x(0.5 * (_vertices[e_1].x + _vertices[e_2].x), 0.5 * (_vertices[e_1].y + _vertices[e_2].y));
            if(point_in_polygon(x,polygon)){
                is_inside_polygon[n]=(true);
                continue;
            }
            else {
                //std::cout<<e.src<<"----------"<<e.dst<<std::endl;
                is_inside_polygon[n]=(false);
                continue;
            }
        }

        for(int i=0;i<_edge_list.size();i++){
            for(int j=0;j<_edge_list.size();j++){
                if(!is_inside_polygon[i]){
                    _allowed_flips[i][j]=0;
                }

            }
        }
    }


    std::vector<std::vector<bool>> _allowed_flips;
    std::vector<std::vector<triangle_int>> _edge_distance_matrix;
private:
    int n=0;
    int m=0;
    std::vector<Point_2> _vertices;
    std::vector<Edge> _edge_list;
    std::vector<std::vector<int>> _edges_matrix;
    std::vector<Triangle> _triangle_list;
    std::vector<bool> _triangle_is_allowed_list;
    std::vector<std::vector<std::vector<int>>> _triangle_matrix;
    std::vector<triangle_int> _edge_distance_flattened;
    int _edge_matrix_size;


};


#endif //A_STAR_FOR_FLIPDISTANCE_GEOMETRICOBJECTS_HPP
