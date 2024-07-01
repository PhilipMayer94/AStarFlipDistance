
#ifndef HIGHERORDERDELAUNAY_COMBINATORIALDELAUNAYTRIANGULATION_HPP
#define HIGHERORDERDELAUNAY_COMBINATORIALDELAUNAYTRIANGULATION_HPP


#include <vector>
#include <AStarFlipDistance/BasicDataStructures.hpp>
#include<iostream>
#include<algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>


template<class T>
class CombinatorialDelaunayTriangulation {
public:
    CombinatorialDelaunayTriangulation() = default;

    explicit CombinatorialDelaunayTriangulation(std::vector<T> _vertices) : vertices(_vertices) {
    }

    void triangulate();


    std::pair<int, Edge> edge_triangle_intersection(int e_1, int e_2, int tr, int forbidden_triangle);


    void print_triangles();

    bool has_boundary = true;
    int boundary_size;

    std::vector<T> vertices;
    std::vector<Triangle> triangles;
    std::vector<std::vector<int>> triangles_adjacent_to_vertices;
    std::vector<std::vector<int>> neighbors_of_triangles;
    std::vector<std::vector<Edge>> shared_edge_of_triangle_with_neighbor;
    std::vector<std::vector<int>> adjacency_matrix_of_edges;


private:

    void generate_adjacency_matrix();

    Edge get_shared_edge(int tr1, int tr2);

    void triangulate(std::vector<T> const &_vertices);


};


template<class T>
std::pair<int, Edge>
CombinatorialDelaunayTriangulation<T>::edge_triangle_intersection(int e_1, int e_2, int tr, int forbidden_triangle) {
    std::pair<int, Edge> result(-1, Edge(-1, -1));
    T e1 = vertices[e_1];
    T e2 = vertices[e_2];
    Triangle t = triangles[tr];
    for (int i = 0; i < neighbors_of_triangles[tr].size(); i++) {
        int tr_current = neighbors_of_triangles[tr][i];
        if (tr_current != forbidden_triangle) {
            Edge e = shared_edge_of_triangle_with_neighbor[tr][i];
            if (do_arcs_intersect(vertices[e.src], vertices[e.dst], vertices[e_1], vertices[e_2])) {
                result.first = tr_current, result.second = e;
                return result;
            }
        }
    }


    return result;
}

template<class T>
Edge CombinatorialDelaunayTriangulation<T>::get_shared_edge(int tr1, int tr2) {
    std::array<int, 2> e = {-1, -1};
    Triangle t1 = triangles[tr1];
    Triangle t2 = triangles[tr2];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (t1.vertices[i] == t2.vertices[j]) {
                if (e[0] < 0) {
                    e[0] = t1.vertices[i];
                } else {
                    e[1] = t1.vertices[i];
                }
            }
        }
    }
    if (e[1] < 0) {
        return Edge(-1, -1);
    }


    return Edge(e[0], e[1]);
}


template<class T>
void CombinatorialDelaunayTriangulation<T>::generate_adjacency_matrix() {
    adjacency_matrix_of_edges.clear();
    adjacency_matrix_of_edges.resize(vertices.size());
    for (int i = 0; i < adjacency_matrix_of_edges.size(); i++) {
        adjacency_matrix_of_edges[i].resize(vertices.size());
        std::fill(adjacency_matrix_of_edges[i].begin(), adjacency_matrix_of_edges[i].end(), -1);
    }

    for (auto t: triangles) {
        int e_1 = t.vertices[0];
        int e_2 = t.vertices[1];
        int e_3 = t.vertices[2];

        int minimum = std::min(e_1, e_2);
        int maximum = std::max(e_1, e_2);
        adjacency_matrix_of_edges[minimum][maximum] = 1;

        minimum = std::min(e_1, e_3);
        maximum = std::max(e_1, e_3);
        adjacency_matrix_of_edges[minimum][maximum] = 1;

        minimum = std::min(e_3, e_2);
        maximum = std::max(e_3, e_2);
        adjacency_matrix_of_edges[minimum][maximum] = 1;

    }
    int counter = 0;
    for (int i = 0; i < adjacency_matrix_of_edges.size(); i++) {
        for (auto j: adjacency_matrix_of_edges[i]) {
            if (j == 1) {
                counter++;
            }
        }
    }
}

template<class T>
void CombinatorialDelaunayTriangulation<T>::print_triangles() {
    std::cout << "Delauanay triangles: (Nr.: " << triangles.size() << ")";
    std::cout << std::endl;
    for (auto &triangle: triangles) {
        triangle.print();
        std::cout << std::endl;
    }
}

template<class T>
void CombinatorialDelaunayTriangulation<T>::triangulate(std::vector<T> const &_vertices) {
    std::cout << "wrong datatype in Triangulation";
}

template<class T>
void CombinatorialDelaunayTriangulation<T>::triangulate() {
    triangulate(vertices);
}


template<>
void CombinatorialDelaunayTriangulation<Point_2>::triangulate(std::vector<Point_2> const &_vertices);

template<>
void CombinatorialDelaunayTriangulation<Point_s>::triangulate(std::vector<Point_s> const &_vertices);


#endif //HIGHERORDERDELAUNAY_COMBINATORIALDELAUNAYTRIANGULATION_HPP
