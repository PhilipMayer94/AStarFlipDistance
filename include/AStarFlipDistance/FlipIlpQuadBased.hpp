
#ifndef HIGHERORDERDELAUNAY_FLIPILPQUADBASED_HPP
#define HIGHERORDERDELAUNAY_FLIPILPQUADBASED_HPP

#include <gurobi_c++.h>
#include "BasicDataStructures.hpp"
#include "HigherOrderDelaunayObjects.hpp"
#include "AStarFlipDistance/FlipIlp.hpp"
#include "AStarFlipDistance/GeometricTriangulation.hpp"


template<class T>
class FlipIlpQuadBased : FlipIlp<T> {
    using FlipIlp<T>::layers;
    using FlipIlp<T>::triangles;
    using FlipIlp<T>::edges;
    using FlipIlp<T>::left_triangles_adjacent_to_edge;
    using FlipIlp<T>::right_triangles_adjacent_to_edge;
    using FlipIlp<T>::triangle_weights;
    using FlipIlp<T>::vertices;
    using FlipIlp<T>::m;

    using FlipIlp<T>::has_boundary;

    using FlipIlp<T>::start_triangulation;
    using FlipIlp<T>::target_triangulation;

    using FlipIlp<T>::model;
    using FlipIlp<T>::optimization_result;

public:
    FlipIlpQuadBased(GeometricTriangulation<T> triangulation_1, GeometricTriangulation<T> triangulation_2,
                     HigherOrderDelaunayObjects<T> _objects,bool only_useful_flips=true);

    void solve();

    double get_number_of_flips();

    void print_solution();

    //double result = -1;

    std::vector<Edge> get_e_vars_one_layer(int layer);

    std::vector<Edge> get_d_vars_one_layer(int layer);

    std::pair<std::vector<std::vector<Edge>>, std::vector<std::vector<Edge>>> get_solution();


private:

    int number_off_edges_per_triangulation{};

    void create_model();

    void calculate_all_intersecting_edges();

    void calculate_all_valid_flips();

    std::vector<std::vector<int>> get_adjacency_matrix_with_indices();

    bool is_flip_valid(int a, int b, int c, int d);

    std::pair<std::vector<int>,std::vector<int>> get_relevant_flips(int edge);

    std::vector<int> get_triangulation_edge_indices(GeometricTriangulation<T> triangulation);

    std::vector<std::vector<int>> edge_intersections;
    std::vector<std::pair<int, int>> flips;

    std::vector<std::vector<GRBVar>> e_vars;
    std::vector<std::vector<GRBVar>> a_vars;
    std::vector<std::vector<GRBVar>> f_vars;


    void calculate_all_useful_valid_flips();

    bool triangle_in_vector_of_triangles(Triangle t, std::vector<int> vec);
};


template<class T>
void FlipIlpQuadBased<T>::solve() {
    try {
        model->optimize();
        optimization_result = model->get(GRB_DoubleAttr_ObjVal);
        this->runtime=model->get(GRB_DoubleAttr_Runtime);
        std::cout << "resulting number of flips: " << get_number_of_flips() <<" rt: "<<this->runtime <<std::endl;
    }
    catch (GRBException &e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch (...) {
        std::cout << "Error during optimization" << std::endl;
    }


}

template<class T>
double FlipIlpQuadBased<T>::get_number_of_flips() {

    return optimization_result;
}

template<class T>
void FlipIlpQuadBased<T>::print_solution() {
    std::pair<std::vector<std::vector<Edge>>, std::vector<std::vector<Edge>>> solution = get_solution();
    std::cout << "Maximum number of flips Layers: " << layers << std::endl;
    std::cout << "Number of Flips:  " << get_number_of_flips() << std::endl;
    for (int i = 0; i < solution.first.size(); i++) {

        if (i > 0) {
            std::cout << "D_" << i << ": ";
            for (auto t: solution.second[i]) {
                t.print();
            }
        }
        std::cout << std::endl;
        std::cout << "E_" << i << ": ";
        for (auto t: solution.first[i]) {
            t.print();
        }
        std::cout << std::endl;


    }
}

template<class T>
void FlipIlpQuadBased<T>::create_model() {
    GRBEnv env;
    env.set(GRB_IntParam_OutputFlag, 0);
    model = new GRBModel(env);
    model->set(GRB_DoubleParam_TimeLimit, 180.0);
    for (int i = 0; i < this->layers + 1; i++) {
        e_vars.emplace_back(std::vector<GRBVar>{});
        a_vars.emplace_back(std::vector<GRBVar>{});
        f_vars.emplace_back(std::vector<GRBVar>{});
        for (int j = 0; j < edges.size(); j++) {
            std::string str = "E_";
            str += std::to_string(i) + ". " + std::to_string(j);

            e_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, str));

            str = "A_";
            str += std::to_string(i) + ". " + std::to_string(j);
            a_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY, str));
        }
        for(int j=0;j<flips.size();j++){
            std::string str = "F_";
            str += std::to_string(i) + ". " + std::to_string(j);
            f_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY, str));
        }
    }

    for (auto e_vars_on_layer: e_vars) {
        int layer = -1;
        layer++;
        GRBLinExpr expr_sum;
        for (int i = 0; i < edges.size(); i++) {
            expr_sum += e_vars_on_layer[i];

            for (int j = 0; j < edge_intersections[i].size(); j++) {
                GRBLinExpr expr;
                expr = e_vars_on_layer[i] + e_vars_on_layer[edge_intersections[i][j]];
                std::string str = "C_" + std::to_string(layer) + "__" + std::to_string(i) + "_" +
                                  std::to_string(edge_intersections[i][j]);
                model->addConstr(expr, GRB_LESS_EQUAL, 1, str);
            }

        }
        std::string str = "C_" + std::to_string(layer) + "__" + "sum";

        model->addConstr(expr_sum, GRB_EQUAL, number_off_edges_per_triangulation, str);
    }


    //setting start and target triangulation values
    std::vector<int> start_indices = get_triangulation_edge_indices(start_triangulation);
    std::vector<int> target_indices = get_triangulation_edge_indices(target_triangulation);
    std::vector<int> start_values(edges.size(), 0);
    std::vector<int> target_values(edges.size(), 0);
    for (int i = 0; i < start_indices.size(); i++) {
        start_values[start_indices[i]] = 1;
        target_values[target_indices[i]] = 1;
    }

    for (int i = 0; i < edges.size(); i++) {
        GRBLinExpr expr_start;
        GRBLinExpr expr_target;
        expr_start += e_vars[0][i];
        expr_target += e_vars[e_vars.size() - 1][i];
        std::string str = "C_start_" + std::to_string(i);
        model->addConstr(expr_start, GRB_EQUAL, start_values[i], str);
        str = "C_target_" + std::to_string(i);
        model->addConstr(expr_target, GRB_EQUAL, target_values[i], str);
    }

   //d values
    GRBLinExpr expr_objective;
    for (int i = 0; i < e_vars.size() - 1; i++) {
        GRBLinExpr expr_all;


        for(int j=0;j<flips.size();j++){
            expr_objective += f_vars[i][j];
            expr_all += f_vars[i][j];
        }

        for (int j = 0; j < edges.size(); j++) {
            std::pair<std::vector<int>,std::vector<int>>  relevant_flips= get_relevant_flips(j);

            GRBLinExpr expr_in;
            GRBLinExpr expr_out;
            expr_in+=e_vars[i][j]-e_vars[i+1][j];
            expr_out+=e_vars[i+1][j]-e_vars[i][j];
            for(int f : relevant_flips.first){
                expr_in-=f_vars[i][f];
            }

            for(int f : relevant_flips.second){
                expr_out-=f_vars[i][f];
            }

            std::string str = "C_e_out_" + std::to_string(i) + "__" + std::to_string(j);
            model->addConstr(expr_out, GRB_LESS_EQUAL, 0, str);

            str = "C_e_in_" + std::to_string(i) + "__" + std::to_string(j);
            //model->addConstr(expr_in, GRB_LESS_EQUAL, 0, str);
        }
        std::string str = "C_d_sum_" + std::to_string(i);
        model->addConstr(expr_all, GRB_LESS_EQUAL, 1, str);
    }
    model->setObjective(expr_objective, GRB_MINIMIZE);
}

template<class T>
void FlipIlpQuadBased<T>::calculate_all_intersecting_edges() {
    for (int i = 0; i < edges.size(); i++) {
        edge_intersections.emplace_back(std::vector<int>());
    }
    for (int i = 0; i < edges.size(); i++) {
        for (int j = i + 1; j < edges.size(); j++) {
            if (do_arcs_intersect(vertices[edges[i].src], vertices[edges[i].dst], vertices[edges[j].src],
                                  vertices[edges[j].dst])) {
                edge_intersections[i].emplace_back(j);
                edge_intersections[j].emplace_back(i);
            }
        }
    }

}


template<class T>
std::vector<int> FlipIlpQuadBased<T>::get_triangulation_edge_indices(GeometricTriangulation<T> triangulation) {
    std::vector<int> result_indices;
    for (auto e: triangulation.edges) {
        for (int i = 0; i < edges.size(); i++) {
            if (e == edges[i]) {
                result_indices.emplace_back(i);
            }
        }
    }
    return result_indices;
}

template<class T>
FlipIlpQuadBased<T>::FlipIlpQuadBased(GeometricTriangulation<T> triangulation_1,
                                      GeometricTriangulation<T> triangulation_2,
                                      HigherOrderDelaunayObjects<T> _objects,bool only_useful_flips) : FlipIlp<T>(triangulation_1,
                                                                                           triangulation_2,
                                                                                           _objects) {
    number_off_edges_per_triangulation = 3 * _objects.get_vertices().size() - 3 - _objects.get_size_of_convex_hull();
    if (!has_boundary) {
        number_off_edges_per_triangulation = number_off_edges_per_triangulation - 3;
    }

    calculate_all_intersecting_edges();
    if(only_useful_flips){
        calculate_all_useful_valid_flips();
    }
    else{
        calculate_all_valid_flips();
    }

    create_model();

}

template<class T>
std::vector<Edge> FlipIlpQuadBased<T>::get_e_vars_one_layer(int layer) {
    std::vector<Edge> result_triangulation;
    for (int i = 0; i < edges.size(); i++) {
        if (e_vars[layer][i].get(GRB_DoubleAttr_X) > 0.9) {
            result_triangulation.emplace_back(edges[i]);
        }
    }

    return result_triangulation;
}


template<class T>
std::vector<Edge> FlipIlpQuadBased<T>::get_d_vars_one_layer(int layer) {
    std::vector<Edge> result_triangulation;
    for (int i = 0; i < edges.size(); i++) {
        if (f_vars[layer][i].get(GRB_DoubleAttr_X) > 0.9) {
            result_triangulation.emplace_back(edges[i]);
        }
    }

    return result_triangulation;
}

template<class T>
std::pair<std::vector<std::vector<Edge>>, std::vector<std::vector<Edge>>>
FlipIlpQuadBased<T>::get_solution() {
    std::pair<std::vector<std::vector<Edge>>, std::vector<std::vector<Edge>>> result_pair;


    for (int i = 0; i < layers; i++) {
        result_pair.second.emplace_back(get_d_vars_one_layer(i));
        result_pair.first.emplace_back(get_e_vars_one_layer(i));
    }

    return result_pair;
}





template<class T>
void FlipIlpQuadBased<T>::calculate_all_valid_flips() {
    get_adjacency_matrix_with_indices() ;
    flips.clear();
    for (int i = 0; i < edge_intersections.size(); i++) {
        for (int j = 0; j < edge_intersections[i].size(); j++) {
            int a = edges[edge_intersections[i][j]].src;
            int b= edges[edge_intersections[i][j]].dst;
            int c= edges[i].src;
            int d = edges[i].dst;
            if (is_flip_valid(a,b,c,d)) {
                flips.emplace_back(i,edge_intersections[i][j]);
            }
        }
    }
}

int get_third_point(Edge e,Triangle t){
    for(int i=0;i<3;i++){
        if(t.vertices[i]!=e.src&&t.vertices[i]!=e.dst){
            return t.vertices[i];
        }
    }
    return -1;
}

template<class T>
bool FlipIlpQuadBased<T>::triangle_in_vector_of_triangles(Triangle t,std::vector<int> vec){
    for(int & i : vec){
        if(t==triangles[i]){
            return true;
        }
    }
    return false;
}

template<class T>
void FlipIlpQuadBased<T>::calculate_all_useful_valid_flips() {
    std::vector<std::vector<int>> matrix= get_adjacency_matrix_with_indices();
    flips.clear();
    for(int i=0;i<edges.size();i++){
        int e_1=i;
        for(auto t_1: left_triangles_adjacent_to_edge[e_1]){
            for(auto t_2: right_triangles_adjacent_to_edge[e_1]){
                int v_1= get_third_point(edges[e_1],triangles[t_1]);
                int v_2= get_third_point(edges[e_1],triangles[t_2]);
                int e_2=-1;
                if(matrix[v_1][v_2]<0){
                    continue;
                }
                else{
                    e_2=matrix[v_1][v_2];
                }
                Triangle quad_t_1(v_1,v_2,edges[e_1].src);
                Triangle quad_t_2(v_1,v_2,edges[e_1].dst);
                if( (triangle_in_vector_of_triangles(quad_t_1,left_triangles_adjacent_to_edge[e_2]) &&
                             triangle_in_vector_of_triangles(quad_t_2,right_triangles_adjacent_to_edge[e_2]) )
                             || (triangle_in_vector_of_triangles(quad_t_2,left_triangles_adjacent_to_edge[e_2]) &&
                                 triangle_in_vector_of_triangles(quad_t_1,right_triangles_adjacent_to_edge[e_2])
                                                     )){
                    int a = edges[e_1].src;
                    int b= edges[e_1].dst;
                    int c= edges[e_2].src;
                    int d = edges[e_2].dst;

                    if (is_flip_valid(a,b,c,d)) {
                        flips.emplace_back(e_1,e_2);
                    }



                }
            }
        }
    }


}



template<class T>
std::pair<std::vector<int>,std::vector<int>> FlipIlpQuadBased<T>::get_relevant_flips(int edge) {
    std::pair<std::vector<int>,std::vector<int>> result_flips;
    for(int i=0;i<flips.size();i++){
        if(flips[i].first==edge){
            result_flips.first.emplace_back(i);
        }
        if(flips[i].second==edge){
            result_flips.second.emplace_back(i);
        }
    }


    return result_flips;
}

template<class T>
bool FlipIlpQuadBased<T>::is_flip_valid(int a, int b, int c, int d){
    if(!is_quadrailateral_convex(vertices[a],vertices[b],vertices[c],vertices[d])){
        return false;
    }
    for(int i=0;i<vertices.size();i++){
        if(i==a||i==b||i==c||i==d){
            continue;
        }
        else{
            if(is_point_in_triangle(vertices[i],vertices[b],vertices[c],vertices[a])){
                return false;
            }
            if(is_point_in_triangle(vertices[i],vertices[b],vertices[d],vertices[a])){
                return false;
            }
        }
    }
    return true;
}


template<class T>
std::vector<std::vector<int>> FlipIlpQuadBased<T>::get_adjacency_matrix_with_indices() {
    std::vector<std::vector<int>> result_matrix;
    for(int i=0;i<vertices.size();i++){
        result_matrix.emplace_back();
        result_matrix[i].resize(vertices.size());
        std::fill(result_matrix[i].begin(), result_matrix[i].end(),-1);
    }
    for(int i=0;i<edges.size();i++){
        result_matrix[edges[i].src][edges[i].dst]=i;
        result_matrix[edges[i].dst][edges[i].src]=i;
    }
    return result_matrix;
}



#endif //HIGHERORDERDELAUNAY_FLIPILPQUADBASED_HPP
