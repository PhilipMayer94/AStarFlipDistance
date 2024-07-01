
#ifndef HIGHERORDERDELAUNAY_HIGHERORDERDELAUNAYOBJECTS_HPP
#define HIGHERORDERDELAUNAY_HIGHERORDERDELAUNAYOBJECTS_HPP

#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>

#include <AStarFlipDistance/CombinatorialDelaunayTriangulation.hpp>
#include <AStarFlipDistance/Geometry.hpp>
#include "AStarFlipDistance/BasicDataStructures.hpp"

template<class T>
class HigherOrderDelaunayObjects {


public:

    HigherOrderDelaunayObjects() = default;

    explicit HigherOrderDelaunayObjects(std::vector<T> _data_points, bool with_pre_computation = true);

    void compute_useful_order_k_edges_and_triangles(int order) ;

    std::vector<int>  compute_order_k_objects_already_pre_computed(int order);

    int get_useful_order_of_triangle(Triangle t);

    std::tuple<std::vector<Edge>, std::vector<Triangle>, std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>>> get_polygon_objects();

    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> get_left_and_right_adjacent_triangles();

    std::vector<Edge> get_useful_edges();

    std::vector<Triangle> get_useful_triangles();

    std::vector<Triangle> get_delaunay_triangles();

    std::vector<T> get_vertices();

    int get_size_of_convex_hull();

    std::vector<Edge> get_fixed_edges();


    std::vector<int> get_delauany_triangulation_as_indices();

    std::vector<bool>  get_which_edges_are_inside_simple_polygon();

    bool manifold_is_bounded();




    void print_numbers();

    void print_edges_with_adjacent_triangles();

    void print_triangles();

    void print_triangle_with_order();

    void print_edges();




private:
    bool defining_circle_has_order_k(int v, int e_1, int e_2, int order, std::vector<int> &candidates);

    bool defining_circle_has_order_k(int v, int e_1, int e_2, int order,
                                     std::vector<int> &candidates, int &useful_order);

    void compute_order_k_useful_edges(int order);


    void compute_order_k_triangles(int order);

    void compute_useful_triangles_one_edge(int e, int order,
                                           std::unordered_set<std::pair<Triangle, int>, triangle_pair_hash> &result);

    bool is_edge_useful(int e_1, int e_2, int order, std::vector<int> &candidates, int &useful_order);

    bool hull_kod_test(int e_1, int e_2, int order, std::vector<int> &all_hull_vertices, std::vector<int> &candidates,
                       int &useful_order);

    std::vector<int> collect_delaunay_vertices(int e_1, int e_2, int order);

    int get_defining_point_on_one_side(int e_1, int e_2, std::vector<int> &candidates);

    std::pair<int, Edge> get_next_delaunay_triangle(int e1, int e2, int triangle, int previous_triangle);

    bool triangle_has_order_k_relative_to_set(Triangle t, int order, std::vector<int> &candidates, int &useful_order);

    bool circle_is_defining(int current_vertex, std::vector<double> circle, std::vector<int> &candidates);

    bool is_vector_not_in_triangle(Triangle t, std::vector<int> &vector);

    std::vector<std::vector<int>> get_adjacency_matrix();

    void add_one_useful_triangle(std::pair<Triangle, int> p, std::vector<std::vector<int>> &adjacency_matrix);



    int size_of_convex_hull = -1;

    int last_calculated_order = -1;
    bool pre_computation = true;
    std::vector<T> vertices;

    std::vector<Triangle> result_triangles;
    std::vector<Edge> result_edges;
    std::vector<std::vector<int>> triangles_adjacent_to_useful_edges;


    std::vector<std::vector<int>> candidate_vertices;

    std::vector<int> edges_useful_order;
    std::vector<int> triangles_useful_order;

    CombinatorialDelaunayTriangulation<T> delaunay;



};

template<class T>
std::vector<Edge> HigherOrderDelaunayObjects<T>::get_fixed_edges() {
    std::vector<Edge> result;
    for(auto e_1 : result_edges){
        int counter=0;
        for(auto e_2 : result_edges){
            if(e_1.src==e_2.src ||
                    e_1.src==e_2.dst ||
                    e_1.dst==e_2.src ||
                    e_1.dst==e_2.dst){
                continue;
            }
            if(do_arcs_intersect(vertices[e_1.src],vertices[e_1.dst],vertices[e_2.src],vertices[e_2.dst])){
                counter++;
            }
        }
        if(counter==0){
            result.emplace_back(e_1);
        }
    }
    return result;
}


template<class T>
void HigherOrderDelaunayObjects<T>::compute_order_k_useful_edges(int order) {

    for (int e_1 = 0; e_1 < vertices.size(); e_1++) {
        for (int e_2 = e_1 + 1; e_2 < vertices.size(); e_2++) {
            std::vector<int> current_candidates = {};
            int useful_order;
            if (is_edge_useful(e_1, e_2, order, current_candidates, useful_order)) {
                result_edges.emplace_back(e_1, e_2);
                edges_useful_order.emplace_back(useful_order);
                candidate_vertices.emplace_back(current_candidates);
                triangles_adjacent_to_useful_edges.emplace_back(std::vector<int>());
            }
        }
    }

}

template<class T>
bool HigherOrderDelaunayObjects<T>::is_edge_useful(int e_1, int e_2, int order, std::vector<int> &candidates,
                                                   int &useful_order) {
    std::vector<int> intersected_delaunay_endpoints = collect_delaunay_vertices(e_1, e_2, order);

    if (intersected_delaunay_endpoints.empty()) {
        useful_order = 0;
        return true;
    }
    if (intersected_delaunay_endpoints.at(0) == -INT_MAX) {
        return false;
    }

    bool result = hull_kod_test(e_1, e_2, order, intersected_delaunay_endpoints, candidates, useful_order);

    return result;
}

template<class T>
std::vector<int> HigherOrderDelaunayObjects<T>::collect_delaunay_vertices(int e_1, int e_2, int order) {
    std::vector<int> result;
    result.reserve(order * 2 + 3);
    int previous_triangle = -1;
    std::pair<int, Edge> intersection(-1, Edge(-1, -1));

    for (int i = 0; i < delaunay.triangles_adjacent_to_vertices[e_1].size(); i++) {
        int tr = delaunay.triangles_adjacent_to_vertices[e_1][i];

        for (int k = 0; k < 2; k++) {
            if (delaunay.triangles[tr].vertices[k] == e_2) {
                return {};
            }
        }


        intersection = delaunay.edge_triangle_intersection(e_1, e_2, tr, -10);
        if (intersection.first >= 0) {
            previous_triangle = tr;
            break;
        }
    }

    if (intersection.first < 0) {
        return {};
    }
    result.push_back(intersection.second.src);
    result.push_back(intersection.second.dst);

    int current_triangle = intersection.first;
    int i = 0;
    while (true) {
        if (delaunay.triangles[current_triangle].vertices[0] == e_2
            || delaunay.triangles[current_triangle].vertices[1] == e_2
            || delaunay.triangles[current_triangle].vertices[2] == e_2) {
            return result;
        }

        std::pair<int, Edge> next_triangle = get_next_delaunay_triangle(e_1, e_2, current_triangle, previous_triangle);


        int tmp;

        if (next_triangle.first < 0) {
            std::cout << e_1<<" "<<e_2<< " edge is colinear with Delaunay edge" << std::endl;

            result.clear();
            result.emplace_back(-INT_MAX);
            return result;

        }

        int v_1 = next_triangle.second.src;
        int v_2 = next_triangle.second.dst;
        if (v_1 == result[result.size() - 2]) {
            tmp = result[result.size() - 2];
            result[result.size() - 2] = result[result.size() - 1];
            result[result.size() - 1] = tmp;
            result.emplace_back(v_2);
        } else if (v_2 == result[result.size() - 2]) {
            tmp = result[result.size() - 2];
            result[result.size() - 2] = result[result.size() - 1];
            result[result.size() - 1] = tmp;
            result.emplace_back(v_1);
        } else if (v_1 == result[result.size() - 1]) {
            result.emplace_back(v_2);
        } else if (v_2 == result[result.size() - 1]) {
            result.emplace_back(v_1);
        } else {
            result.emplace_back(v_1);
            result.emplace_back(v_2);
        }

        previous_triangle = current_triangle;
        current_triangle = next_triangle.first;

        i++;
        if (i > 2 * order + 2) {
            return result;
        }
    }

    return result;
}

template<class T>
std::pair<int, Edge>
HigherOrderDelaunayObjects<T>::get_next_delaunay_triangle(int e1, int e2, int triangle, int previous_triangle) {
    auto intersection = delaunay.edge_triangle_intersection(e1, e2, triangle, previous_triangle);
    std::pair<int, Edge> tmp = intersection;

    return tmp;
}

template<class T>
bool HigherOrderDelaunayObjects<T>::hull_kod_test(int e_1, int e_2, int order, std::vector<int> &all_hull_vertices,
                                                  std::vector<int> &candidates, int &useful_order) {
    std::vector<int> left_hull;
    left_hull.reserve(2 * order);
    std::vector<int> right_hull;
    right_hull.reserve(2 * order);
    for (auto v: all_hull_vertices) {
        if (is_left_oriented(vertices[v], vertices[e_1], vertices[e_2])) {
            left_hull.emplace_back(v);
        } else {
            right_hull.emplace_back(v);
        }
    }

    if (left_hull.size() > order || right_hull.size() > order) {
        return false;
    }

    int l_defining = get_defining_point_on_one_side(e_1, e_2, left_hull);
    int r_defining = get_defining_point_on_one_side(e_1, e_2, right_hull);

    int l_order = INT_MAX;
    int r_order = INT_MAX;


    bool l_kod = defining_circle_has_order_k(l_defining, e_1, e_2, order, candidates, l_order);
    bool r_kod = defining_circle_has_order_k(r_defining, e_1, e_2, order, candidates, r_order);

    useful_order = std::max(l_order, r_order);

    return l_kod && r_kod;
}

template<class T>
int HigherOrderDelaunayObjects<T>::get_defining_point_on_one_side(int e_1, int e_2, std::vector<int> &candidates) {
    for (auto v: candidates) {
        std::vector<double> circle = circle_through_three_points(vertices[v], vertices[e_1], vertices[e_2]);
        if (circle_is_defining(v, circle, candidates)) {
            return v;
        }
    }
    std::cout << "hullerror";
    return -1;
}

template<class T>
bool HigherOrderDelaunayObjects<T>::circle_is_defining(int current_vertex, std::vector<double> circle,
                                                       std::vector<int> &candidates) {
    for (auto v: candidates) {
        if (v == current_vertex) {
            continue;
        }
        //risky if a candidate lies on the circle but in the context its fine
        if (is_point_in_circle(vertices[v], circle)) {
            return false;
        }
    }
    return true;
}


template<class T>
bool HigherOrderDelaunayObjects<T>::defining_circle_has_order_k(int v, int e_1, int e_2, int order,
                                                                std::vector<int> &candidates) {
    int counter = 0;

    std::vector<double> circle = circle_through_three_points(vertices[v], vertices[e_1], vertices[e_2]);
    for (int i = 0; i < vertices.size(); i++) {
        if (is_point_in_circle(vertices[i], circle) && i != e_1 && i != e_2 && i != v) {
            counter++;
            candidates.emplace_back(i);
            if (counter > order) {
                return false;
            }
        }
    }

    if (counter > order) {
        return false;
    } else {
        return true;
    }

}


template<class T>
bool HigherOrderDelaunayObjects<T>::defining_circle_has_order_k(int v, int e_1, int e_2, int order,
                                                                std::vector<int> &candidates, int &useful_order) {
    int counter = 0;

    std::vector<double> circle = circle_through_three_points(vertices[v], vertices[e_1], vertices[e_2]);
    for (int i = 0; i < vertices.size(); i++) {
        if (is_point_in_circle(vertices[i], circle) && i != e_1 && i != e_2 && i != v) {
            counter++;
            candidates.emplace_back(i);
            if (counter > order) {

                return false;
            }
        }
    }

    if (counter > order) {
        return false;
    } else {
        useful_order = counter;
        return true;
    }

}

template<class T>
void HigherOrderDelaunayObjects<T>::compute_order_k_triangles(int order) {
    std::vector<std::vector<int>> adjacency_matrix = get_adjacency_matrix();

    for (auto t: delaunay.triangles) {
        add_one_useful_triangle(std::pair<Triangle, int>(t, 0), adjacency_matrix);

    }

    std::unordered_set<std::pair<Triangle, int>, triangle_pair_hash> result;

    for (int e = 0; e < result_edges.size(); e++) {
        compute_useful_triangles_one_edge(e, order, result);
    }
    for (auto t: result) {
        add_one_useful_triangle(t, adjacency_matrix);
    }

}

template<class T>
void HigherOrderDelaunayObjects<T>::compute_useful_triangles_one_edge(int e, int order,
                                                                      std::unordered_set<std::pair<Triangle, int>, triangle_pair_hash> &result) {
    if (candidate_vertices[e].empty()) {
        return;
    }
    for (auto c: candidate_vertices[e]) {
        Triangle current_triangle(c, result_edges[e].src, result_edges[e].dst);
        int useful_order = INT_MAX;

        if (is_vector_not_in_triangle(current_triangle, candidate_vertices[e]) &&
            triangle_has_order_k_relative_to_set(current_triangle, order, candidate_vertices[e], useful_order)) {
            result.insert(std::pair<Triangle, int>(current_triangle, useful_order));
        }
    }
}

template<class T>
bool HigherOrderDelaunayObjects<T>::triangle_has_order_k_relative_to_set(Triangle t, int order,
                                                                         std::vector<int> &candidates,
                                                                         int &useful_order) {
    std::vector<double> circle = circle_through_three_points(vertices[t.vertices[0]], vertices[t.vertices[1]],
                                                             vertices[t.vertices[2]]);
    if (t == Triangle(7, 9, 10)) {

    }
    int i = 0;
    for (auto c: candidates) {
        if (c == t.vertices[0] || c == t.vertices[1] || c == t.vertices[2]) {
            continue;
        }
        if (is_point_in_circle(vertices[c], circle) &&
            !(c == t.vertices[0] || c == t.vertices[1] || c == t.vertices[2])) {
            i++;
        }
        if (i > order) {
            return false;
        }
    }
   // std::cout<<std::endl;
   // t.print();std::cout <<"| "<<i<<std::endl;
    useful_order = i;
    return true;
}

template<class T>
bool HigherOrderDelaunayObjects<T>::is_vector_not_in_triangle(Triangle t, std::vector<int> &vector) {
    Triangle test(0, 4, 6);
    bool result = true;
    for (auto v: vector) {
        if (v == t.vertices[0] || v == t.vertices[1] || v == t.vertices[2])
            continue;

        if (is_point_in_triangle(vertices[v], vertices[t.vertices[0]], vertices[t.vertices[1]],
                                 vertices[t.vertices[2]])) {
            result = false;
        }
    }
    return result;
}

template<class T>
std::vector<std::vector<int>> HigherOrderDelaunayObjects<T>::get_adjacency_matrix() {
    std::vector<std::vector<int>> adjacency_matrix;
    adjacency_matrix.resize(vertices.size());
    for (int i = 0; i < adjacency_matrix.size(); i++) {
        adjacency_matrix[i].resize(vertices.size());
        std::fill(adjacency_matrix[i].begin(), adjacency_matrix[i].end(), -1);
    }

    for (int i = 0; i < result_edges.size(); i++) {
        int minimum = std::min(result_edges[i].src, result_edges[i].dst);
        int maximum = std::max(result_edges[i].src, result_edges[i].dst);
        adjacency_matrix[minimum][maximum] = i;
    }


    return adjacency_matrix;
}

template<class T>
std::vector<int> HigherOrderDelaunayObjects<T>::compute_order_k_objects_already_pre_computed(int order) {
    std::vector<Edge> old_useful_edges(result_edges);
    std::vector<int> old_edge_useful_order(edges_useful_order);
    std::vector<Triangle> old_triangles(result_triangles);
    std::vector<int> old_triangles_useful_order(triangles_useful_order);
    std::vector<int> indices_in_old_triangle_list;

    result_edges.clear();
    edges_useful_order.clear();

    result_triangles.clear();
    triangles_useful_order.clear();
    triangles_adjacent_to_useful_edges.clear();


    for (int i = 0; i < old_useful_edges.size(); i++) {
        if (old_edge_useful_order[i] <= order) {
            result_edges.emplace_back(old_useful_edges[i]);
            edges_useful_order.emplace_back(old_edge_useful_order[i]);
            triangles_adjacent_to_useful_edges.emplace_back(std::vector<int>());
        }
    }


    std::vector<std::vector<int>> adjacency_matrix = get_adjacency_matrix();
    for (int i = 0; i < old_triangles.size(); i++) {
        if (old_triangles_useful_order[i] <= order) {
            indices_in_old_triangle_list.emplace_back(i);
            add_one_useful_triangle(std::pair<Triangle, int>(old_triangles[i], old_triangles_useful_order[i]),
                                    adjacency_matrix);
        }
    }
    return indices_in_old_triangle_list;

}

template<class T>
void HigherOrderDelaunayObjects<T>::add_one_useful_triangle(std::pair<Triangle, int> p,
                                                            std::vector<std::vector<int>> &adjacency_matrix) {
    Triangle t = p.first;
    int e_1 = -1;
    int e_2 = -1;
    int e_3 = -1;

    int minimum = std::min(t.vertices[0], t.vertices[1]);
    int maximum = std::max(t.vertices[0], t.vertices[1]);
    if (adjacency_matrix[minimum][maximum] >= 0) {
        e_1 = adjacency_matrix[minimum][maximum];
    }
    minimum = std::min(t.vertices[0], t.vertices[2]);
    maximum = std::max(t.vertices[0], t.vertices[2]);
    if (adjacency_matrix[minimum][maximum] >= 0) {
        e_2 = adjacency_matrix[minimum][maximum];
    }
    minimum = std::min(t.vertices[1], t.vertices[2]);
    maximum = std::max(t.vertices[1], t.vertices[2]);
    if (adjacency_matrix[minimum][maximum] >= 0) {
        e_3 = adjacency_matrix[minimum][maximum];
    }
    if (e_1 >= 0 && e_2 >= 0 && e_3 >= 0) {
        triangles_adjacent_to_useful_edges[e_1].emplace_back(result_triangles.size());
        triangles_adjacent_to_useful_edges[e_2].emplace_back(result_triangles.size());
        triangles_adjacent_to_useful_edges[e_3].emplace_back(result_triangles.size());
        result_triangles.emplace_back(t);
        triangles_useful_order.emplace_back(p.second);

    }
}

template<class T>
std::vector<Triangle> HigherOrderDelaunayObjects<T>::get_delaunay_triangles(){
    return delaunay.triangles;

}

template<class T>
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
HigherOrderDelaunayObjects<T>::get_left_and_right_adjacent_triangles() {
    std::vector<std::vector<int>> left;
    std::vector<std::vector<int>> right;
    for (int i = 0; i < triangles_adjacent_to_useful_edges.size(); i++) {
        left.emplace_back(std::vector<int>());
        right.emplace_back(std::vector<int>());
        Edge e = result_edges[i];
        for (int &j: triangles_adjacent_to_useful_edges[i]) {
            Triangle t = result_triangles[j];
            int v = t.third_point(e.src, e.dst);
            if (is_left_oriented(vertices[v], vertices[e.src], vertices[e.dst])) {
                left[i].emplace_back(j);
            } else {
                right[i].emplace_back(j);
            }
        }
    }

    return std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>{left, right};
}

template<class T>
std::vector<Edge> HigherOrderDelaunayObjects<T>::get_useful_edges() {
    return result_edges;
}

template<class T>
std::vector<Triangle> HigherOrderDelaunayObjects<T>::get_useful_triangles() {
    return result_triangles;
}

template<class T>
std::vector<T> HigherOrderDelaunayObjects<T>::get_vertices() {
    return vertices;
}

template<class T>
bool HigherOrderDelaunayObjects<T>::manifold_is_bounded() {

    return delaunay.has_boundary;
}

template<class T>
HigherOrderDelaunayObjects<T>::HigherOrderDelaunayObjects(std::vector<T> _data_points, bool with_pre_computation)
        : vertices(
        _data_points) {
    delaunay = CombinatorialDelaunayTriangulation(vertices);
    delaunay.triangulate();
    pre_computation = with_pre_computation;
    size_of_convex_hull = delaunay.boundary_size;

}

template<class T>
void HigherOrderDelaunayObjects<T>::print_numbers() {
    std::cout << "n:                              " << vertices.size() << std::endl;
    std::cout << "order:                          " << last_calculated_order << std::endl;
    std::cout << "Number of useful Edges:         " << result_edges.size() << std::endl;
    std::cout << "Number of useful Triangles:     " << result_triangles.size() << std::endl;
    std::cout << "Number of Delaunay Triangles:   " << delaunay.triangles.size() << std::endl;
    std::cout << std::endl;
}

template<class T>
void HigherOrderDelaunayObjects<T>::print_edges_with_adjacent_triangles() {
    std::cout << "Edges with adjacent triangles: " << std::endl;
    for (int e = 0; e < result_edges.size(); e++) {
        result_edges[e].print();
        std::cout << ": " << "Nr = " << triangles_adjacent_to_useful_edges[e].size() << "; ";
        for (auto t: triangles_adjacent_to_useful_edges[e]) {
            result_triangles[t].print();
        }
        std::cout << std::endl;
    }
}

template<class T>
void HigherOrderDelaunayObjects<T>::print_triangles() {
    std::cout << "The useful triangles: " << std::endl;
    std::vector<Triangle> inOrder(result_triangles);
    std::sort(inOrder.begin(), inOrder.end());
    for (auto t: inOrder) {
        t.print();
        std::cout << std::endl;
    }
}

template<class T>
void HigherOrderDelaunayObjects<T>::print_triangle_with_order() {
    std::cout << "The useful triangles: " << std::endl;

    for (int i = 0; i < result_triangles.size(); i++) {
        Triangle t = result_triangles[i];
        t.print();
        std::cout << ", " << triangles_useful_order[i] << std::endl;
    }
}

template<class T>
void HigherOrderDelaunayObjects<T>::print_edges() {

    std::cout << "The useful edges: " << std::endl;
    std::vector<Edge> inOrder(result_edges);
    std::sort(inOrder.begin(), inOrder.end());
    for (int i = 0; i < inOrder.size(); i++) {
        Edge e = inOrder[i];
        e.print();
        std::cout << ", " << edges_useful_order[i];
        std::cout << std::endl;
    }
}

template<class T>
int HigherOrderDelaunayObjects<T>::get_useful_order_of_triangle(Triangle t) {
    for (int i = 0; i < result_triangles.size(); i++) {
        if (t == result_triangles[i]) {
            return triangles_useful_order[i];
        }
    }
    return -1;
}

template<class T>
int HigherOrderDelaunayObjects<T>::get_size_of_convex_hull() {
    return size_of_convex_hull;
}

template<class T>
std::vector<int> HigherOrderDelaunayObjects<T>::get_delauany_triangulation_as_indices() {
    std::vector<int> result;
    for(auto t: delaunay.triangles){
        for(int i=0;i<result_triangles.size();i++){
            if(t==result_triangles[i]){
                result.emplace_back(i);
            }
        }
    }
    return result;
}


template<class T>
std::vector<bool>  HigherOrderDelaunayObjects<T>::get_which_edges_are_inside_simple_polygon() {
    std::vector<bool>  edges_inside_simple_polygon;
    for(auto e:result_edges){
        if(edge_in_polygon(e.src,e.dst,vertices)){
            edges_inside_simple_polygon.emplace_back(true);
            //std::cout<<std::endl;
            //e.print();

        }
        else{
            edges_inside_simple_polygon.emplace_back(false);
        }
    }
    return edges_inside_simple_polygon;
}

template<class T>
std::tuple<std::vector<Edge>, std::vector<Triangle>, std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>>>
HigherOrderDelaunayObjects<T>::get_polygon_objects() {
    std::vector<bool>  edges_inside_simple_polygon= get_which_edges_are_inside_simple_polygon();
    std::vector<Edge> edges;
    std::vector<Triangle> triangles;
    std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>> left_right_triangles;

    for(int i =0;i<result_edges.size();i++){
        if(edges_inside_simple_polygon[i]){
            edges.emplace_back(result_edges[i]);
            left_right_triangles.first.emplace_back(std::vector<int>());
            left_right_triangles.second.emplace_back(std::vector<int>());
        }
    }

    std::vector<std::vector<int>> adjacency_matrix;
    adjacency_matrix.resize(vertices.size());
    for (auto & i : adjacency_matrix) {
        i.resize(vertices.size());
        std::fill(i.begin(), i.end(), -1);
    }
    for (int i = 0; i < edges.size(); i++) {
        int minimum = std::min(edges[i].src, edges[i].dst);
        int maximum = std::max(edges[i].src, edges[i].dst);
        adjacency_matrix[minimum][maximum] = i;
    }

    int i=0;
    for(int k=0;k<result_triangles.size();k++){
        auto t =result_triangles[k];
        int e_1,e_2,e_3;
        int v_1,v_2,v_3;
        if(adjacency_matrix[std::min(t.vertices[0],t.vertices[1])][std::max(t.vertices[0],t.vertices[1])]>=0){
            e_1=adjacency_matrix[std::min(t.vertices[0],t.vertices[1])][std::max(t.vertices[0],t.vertices[1])];
            v_1= t.vertices[2];
        }
        else{continue;}
        if(adjacency_matrix[std::min(t.vertices[0],t.vertices[2])][std::max(t.vertices[0],t.vertices[2])]>=0){
            e_2=adjacency_matrix[std::min(t.vertices[0],t.vertices[2])][std::max(t.vertices[0],t.vertices[2])];
            v_2=t.vertices[1];
        }
        else{continue;}
        if(adjacency_matrix[std::min(t.vertices[1],t.vertices[2])][std::max(t.vertices[2],t.vertices[1])]>=0){
            e_3=adjacency_matrix[std::min(t.vertices[1],t.vertices[2])][std::max(t.vertices[2],t.vertices[1])];
            v_3=t.vertices[0];
        }
        else{continue;}


        if(is_left_oriented(vertices[v_1],vertices[edges[e_1].src],vertices[edges[e_1].dst])){
            left_right_triangles.first[e_1].emplace_back(i);
        }
        else{
            left_right_triangles.second[e_1].emplace_back(i);
        }

        if(is_left_oriented(vertices[v_2],vertices[edges[e_2].src],vertices[edges[e_2].dst])){
            left_right_triangles.first[e_2].emplace_back(i);
        }
        else{
            left_right_triangles.second[e_2].emplace_back(i);
        }

        if(is_left_oriented(vertices[v_3],vertices[edges[e_3].src],vertices[edges[e_3].dst])){
            left_right_triangles.first[e_3].emplace_back(i);
        }
        else{
            left_right_triangles.second[e_3].emplace_back(i);
        }
        //t.print();
        //std::cout<< std::endl;
        triangles.emplace_back(t);
        i++;
    }

    return {edges,triangles,left_right_triangles};
}

template<class T>
void HigherOrderDelaunayObjects<T>::compute_useful_order_k_edges_and_triangles(int order) {
    if ( order < last_calculated_order && pre_computation) {
        compute_order_k_objects_already_pre_computed(order);
    } else {
        result_triangles.clear();
        result_edges.clear();
        triangles_adjacent_to_useful_edges.clear();
        candidate_vertices.clear();
        edges_useful_order.clear();
        triangles_useful_order.clear();

        compute_order_k_useful_edges(order);
        compute_order_k_triangles(order);
    }

    last_calculated_order = order;
}


#endif //HIGHERORDERDELAUNAY_HIGHERORDERDELAUNAYOBJECTS_HPP
