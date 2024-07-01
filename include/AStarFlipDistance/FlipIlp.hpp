
#ifndef HIGHERORDERDELAUNAY_FLIPILP_HPP
#define HIGHERORDERDELAUNAY_FLIPILP_HPP

#include <gurobi_c++.h>
#include "BasicDataStructures.hpp"
#include "HigherOrderDelaunayObjects.hpp"
#include "AStarFlipDistance/GeometricTriangulation.hpp"


template<class T>
class FlipIlp {
public:




protected:

    FlipIlp(GeometricTriangulation<T> triangulation_1,
            GeometricTriangulation<T> triangulation_2,
            HigherOrderDelaunayObjects<T> _objects);
    ~FlipIlp(){
        delete model;
    }

    virtual void solve() = 0;

    virtual double get_number_of_flips() = 0;

    virtual void print_solution() = 0;

    virtual void create_model() = 0;

    std::vector<Triangle> triangles;
    std::vector<Edge> edges;
    std::vector<std::vector<int>> left_triangles_adjacent_to_edge;
    std::vector<std::vector<int>> right_triangles_adjacent_to_edge;
    std::vector<double> triangle_weights;
    std::vector<T> vertices;
    int m;
    int layers;

    bool has_boundary;
    double runtime;


    GeometricTriangulation<T> start_triangulation;
    GeometricTriangulation<T> target_triangulation;


    GRBModel *model;
    double optimization_result = -1;
    double primal_bound = -1;

    int get_precise_layer_bound();

    int get_layer_bound();


};

template<class T>
FlipIlp<T>::FlipIlp(GeometricTriangulation<T> triangulation_1,
                    GeometricTriangulation<T> triangulation_2,
                    HigherOrderDelaunayObjects<T> _objects) {
    start_triangulation = triangulation_1;
    target_triangulation = triangulation_2;
    triangles = _objects.get_useful_triangles();
    edges = _objects.get_useful_edges();
    left_triangles_adjacent_to_edge = _objects.get_left_and_right_adjacent_triangles().first;
    right_triangles_adjacent_to_edge = _objects.get_left_and_right_adjacent_triangles().second;
    vertices = _objects.get_vertices();
    m = (2 * vertices.size() - 4);
    has_boundary = _objects.manifold_is_bounded();
    layers = -1;

}

template<class T>
int FlipIlp<T>::get_layer_bound() {
    int counter = 0;
    for (int i = 0; i < start_triangulation.edges.size(); i++) {
        auto e_1 = start_triangulation.edges[i];
        for (int j = 0; j < target_triangulation.edges.size(); j++) {
            auto e_2 = target_triangulation.edges[j];
            if (e_1.src == e_2.src || e_1.src == e_2.dst || e_1.dst == e_2.src || e_1.dst == e_2.dst) {
                continue;
            } else {
                if (do_arcs_intersect(
                        vertices[e_1.src], vertices[e_1.dst], vertices[e_2.src], vertices[e_2.dst])) {
                    counter++;
                }
            }
        }
    }

    return counter;
}

template<class T>
int FlipIlp<T>::get_precise_layer_bound() {
    int layer_counter = 0;

    std::vector<std::vector<int>> adjacency_matrix;
    adjacency_matrix.resize(vertices.size());
    for (auto &i: adjacency_matrix) {
        i.resize(vertices.size());
        std::fill(i.begin(), i.end(), -1);
    }
    for (int i = 0; i < start_triangulation.edges.size(); i++) {
        int minimum = std::min(start_triangulation.edges[i].src, start_triangulation.edges[i].dst);
        int maximum = std::max(start_triangulation.edges[i].src, start_triangulation.edges[i].dst);
        adjacency_matrix[minimum][maximum] = i;
    }


    GeometricTriangulation<T> tmp_triangulation = start_triangulation;
    bool not_finished = true;
    while (not_finished) {
        int index = -1;
        int counter_max = 0;
        for (int i = 0; i < tmp_triangulation.edges.size(); i++) {
            int counter = 0;
            for (int j = 0; j < tmp_triangulation.edges.size(); j++) {
                auto e_1 = tmp_triangulation.edges[i];
                auto e_2 = target_triangulation.edges[j];
                if (e_1.src == e_2.src || e_1.src == e_2.dst || e_1.dst == e_2.src || e_1.dst == e_2.dst) {
                    continue;
                } else {
                    if (do_arcs_intersect(
                            vertices[e_1.src], vertices[e_1.dst], vertices[e_2.src], vertices[e_2.dst])) {
                        counter++;
                    }
                }
            }
            if (counter > counter_max) {
                index = i;
                counter_max = counter;
            }
        }

        Edge e = tmp_triangulation.edges[index];
        std::vector<Triangle> empty_triangles;

        for (int i = 0; i < tmp_triangulation.vertices.size(); i++) {
            if (i != e.src && i != e.dst) {
                if (adjacency_matrix[std::min(e.src, i)][std::max(e.src, i)] >= 0 &&
                    adjacency_matrix[std::min(e.dst, i)][std::max(e.dst, i)] >= 0) {
                    bool is_empty = true;
                    for (int j = 0; j < vertices.size(); j++) {
                        if (j != i && j != e.src && j != e.dst &&
                            is_point_in_triangle(vertices[j], vertices[i], vertices[e.src], vertices[e.dst])) {
                            is_empty = false;
                            break;
                        }
                    }
                    if (is_empty) {
                        empty_triangles.emplace_back(e.src, e.dst, i);
                    }
                }
            }
        }


        int new_src;
        int new_dst;


        if (index < 0) {
            not_finished = false;
        } else {
            layer_counter++;
            for (auto t: empty_triangles) {
                //t.print();
                /*for(int i=0;i<3;i++){
                    std::cout<<"i: "<<i<<" ";
                    vertices[t.vertices[i]].print();
                }*/
            }

            //std::cout << std::endl;

            assert(empty_triangles.size() == 2);
            for (int i = 0; i < 3; i++) {
                if (empty_triangles[0].vertices[i] != e.src && empty_triangles[0].vertices[i] != e.dst) {
                    new_src = empty_triangles[0].vertices[i];
                }
                if (empty_triangles[1].vertices[i] != e.src && empty_triangles[1].vertices[i] != e.dst) {
                    new_dst = empty_triangles[1].vertices[i];
                }
            }
            tmp_triangulation.edges[index] = Edge(new_dst, new_src);
            adjacency_matrix[std::min(new_src, new_dst)][std::max(new_src, new_dst)] = index;
            adjacency_matrix[std::min(e.src, e.dst)][std::max(e.src, e.dst)] = -1;
        }
    }
    return layer_counter;
}

#endif //HIGHERORDERDELAUNAY_FLIPILPEDGEBASED_HPP