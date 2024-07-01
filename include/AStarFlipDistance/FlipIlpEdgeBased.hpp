#include <gurobi_c++.h>
#include "BasicDataStructures.hpp"
#include "HigherOrderDelaunayObjects.hpp"
#include "AStarFlipDistance/FlipIlp.hpp"
#include "AStarFlipDistance/GeometricTriangulation.hpp"


#ifndef HIGHERORDERDELAUNAY_FLIPILPEDGEBASED_HPP
#define HIGHERORDERDELAUNAY_FLIPILPEDGEBASED_HPP




template<class T>
class FlipIlpEdgeBased : FlipIlp<T> {
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
    using FlipIlp<T>::primal_bound;

public:
    FlipIlpEdgeBased(GeometricTriangulation<T> triangulation_1, GeometricTriangulation<T> triangulation_2,
            HigherOrderDelaunayObjects<T> _objects);

    FlipIlpEdgeBased(GeometricTriangulation<T> triangulation_1, GeometricTriangulation<T> triangulation_2,
                     HigherOrderDelaunayObjects<T> _objects, int layer_number);

    void solve();
    double get_number_of_flips();
    double get_primal();
    double get_runtime();
    void print_solution();
    double result = -1;
    double runtime=-1;


    int _timeout_in_minutes=5   *60;

    std::vector<Edge> get_e_vars_one_layer(int layer);
    std::vector<Edge> get_d_vars_one_layer(int layer);
    std::pair<std::vector<std::vector<Edge>>, std::vector<std::vector<Edge>>> get_solution() ;


private:

    int number_off_edges_per_triangulation{};
    void create_model();
    void calculate_all_intersecting_edges();


    std::vector<int> get_triangulation_edge_indices(GeometricTriangulation<T> triangulation);

    std::vector<std::vector<int>> edge_intersections;
    std::vector<std::vector<int>> edge_intersections_both_ways;

    std::vector<std::vector<GRBVar>> e_vars;
    std::vector<std::vector<GRBVar>> a_vars;
    std::vector<std::vector<GRBVar>> d_vars;




};

template<class T>
double FlipIlpEdgeBased<T>::get_primal() {
    return primal_bound;
}

template<class T>
double FlipIlpEdgeBased<T>::get_runtime() {
    return this->runtime;
}


template<class T>
void FlipIlpEdgeBased<T>::solve() {
    try {
        model->optimize();
        this->runtime=model->get(GRB_DoubleAttr_Runtime);
        if(this->runtime>=_timeout_in_minutes){
            optimization_result=-model->get(GRB_DoubleAttr_ObjBound);
            primal_bound       =-model->get(GRB_DoubleAttr_ObjVal);
        }
        else{
            optimization_result = model->get(GRB_DoubleAttr_ObjVal);
            primal_bound=model->get(GRB_DoubleAttr_ObjVal);
        }

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
double FlipIlpEdgeBased<T>::get_number_of_flips() {

    return optimization_result;
}

template<class T>
void FlipIlpEdgeBased<T>::print_solution() {
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
void FlipIlpEdgeBased<T>::create_model() {
    try {
        std::cout<<"number of layers: " <<layers<<"-----------------------------" <<std::endl;
        GRBEnv env;
        env.set(GRB_IntParam_OutputFlag, 1);
        model = new GRBModel(env);
        model->set(GRB_DoubleParam_TimeLimit, _timeout_in_minutes);
        model->set(GRB_IntParam_Threads,1);
        model->set(GRB_StringParam_LogFile, "log_gurobi");
        for (int i = 0; i < this->layers + 1; i++) {
            e_vars.emplace_back(std::vector<GRBVar>{});
            a_vars.emplace_back(std::vector<GRBVar>{});
            d_vars.emplace_back(std::vector<GRBVar>{});
            for (int j = 0; j < edges.size(); j++) {
                std::string str = "E_";
                str += std::to_string(i) + ". " + std::to_string(j);

                e_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, str));

                str = "A_";
                str += std::to_string(i) + ". " + std::to_string(j);
                a_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY, str));
                str = "D_";
                str += std::to_string(i) + ". " + std::to_string(j);
                d_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY, str));
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

//        for (auto e_vars_on_layer: e_vars) {
//            int layer = -1;
//            layer++;
//            for (int i = 0; i < edges.size(); i++) {
//                GRBLinExpr expr;
//                expr+=e_vars_on_layer[i];
//                for (int j = 0; j < edge_intersections_both_ways[i].size(); j++) {
//                    expr +=  e_vars_on_layer[edge_intersections_both_ways[i][j]];
//                }
//                std::string str = "C_" + std::to_string(layer) + "__" + std::to_string(i) + "_D";
//                model->addConstr(expr, GRB_GREATER_EQUAL, 1, str);
//            }
//
//        }


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
            for (int j = 0; j < edges.size(); j++) {
                GRBLinExpr expr;
                expr_objective += d_vars[i][j];
                expr_all += d_vars[i][j];
                expr = expr + e_vars[i][j] - e_vars[i + 1][j] - d_vars[i][j];
                std::string str = "C_d_e_" + std::to_string(i) + "__" + std::to_string(j);
                model->addConstr(expr, GRB_LESS_EQUAL, 0, str);
            }
            std::string str = "C_d_sum_" + std::to_string(i);
            model->addConstr(expr_all, GRB_LESS_EQUAL, 1, str);
        }
        model->setObjective(expr_objective, GRB_MINIMIZE);
    }
    catch (GRBException &e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch (...) {
        std::cout << "Error during optimization" << std::endl;
    }
}

template<class T>
void FlipIlpEdgeBased<T>::calculate_all_intersecting_edges() {
    for(int i=0;i<edges.size();i++){
        edge_intersections.emplace_back(std::vector<int>());
        edge_intersections_both_ways.emplace_back(std::vector<int>());
    }
    for(int i=0;i<edges.size();i++){
        for(int j=i+1;j<edges.size();j++){
            if(do_arcs_intersect(vertices[edges[i].src],vertices[edges[i].dst],vertices[edges[j].src],vertices[edges[j].dst])){
                edge_intersections[i].emplace_back(j);
                //edge_intersections[j].emplace_back(i);

                edge_intersections_both_ways[i].emplace_back(j);
                edge_intersections_both_ways[j].emplace_back(i);
            }
        }
    }

}



template<class T>
std::vector<int> FlipIlpEdgeBased<T>::get_triangulation_edge_indices(GeometricTriangulation<T> triangulation) {
    std::vector<int> result_indices;
    for (auto e: triangulation.edges) {
        for (int i = 0; i < edges.size(); i++) {
            if (e== edges[i]) {
                result_indices.emplace_back(i);
            }
        }
    }
    return result_indices;
}

template<class T>
FlipIlpEdgeBased<T>::FlipIlpEdgeBased(GeometricTriangulation <T> triangulation_1,
                                      GeometricTriangulation <T> triangulation_2,
                                      HigherOrderDelaunayObjects<T> _objects) : FlipIlp<T>(triangulation_1, triangulation_2,
                                                                                      _objects)  {
    number_off_edges_per_triangulation=3*_objects.get_vertices().size()-3-_objects.get_size_of_convex_hull();
    layers=this->get_precise_layer_bound();
    if(!has_boundary){
        number_off_edges_per_triangulation=number_off_edges_per_triangulation-3;
    }
    calculate_all_intersecting_edges();
    create_model();

}

template<class T>
FlipIlpEdgeBased<T>::FlipIlpEdgeBased(GeometricTriangulation <T> triangulation_1,
                                      GeometricTriangulation <T> triangulation_2,
                                      HigherOrderDelaunayObjects<T> _objects, int layer_number) : FlipIlp<T>(triangulation_1, triangulation_2,
                                                                                           _objects)  {
    number_off_edges_per_triangulation=3*_objects.get_vertices().size()-3-_objects.get_size_of_convex_hull();
    layers=layer_number;
    if(!has_boundary){
        number_off_edges_per_triangulation=number_off_edges_per_triangulation-3;
    }
    calculate_all_intersecting_edges();
    create_model();

}

template<class T>
std::vector<Edge> FlipIlpEdgeBased<T>::get_e_vars_one_layer(int layer) {
    std::vector<Edge> result_triangulation;
    for (int i = 0; i < edges.size(); i++) {
        if (e_vars[layer][i].get(GRB_DoubleAttr_X) > 0.9) {
            result_triangulation.emplace_back(edges[i]);
        }
    }

    return result_triangulation;
}


template<class T>
std::vector<Edge> FlipIlpEdgeBased<T>::get_d_vars_one_layer(int layer) {
    std::vector<Edge> result_triangulation;
    for (int i = 0; i < edges.size(); i++) {
        if (d_vars[layer][i].get(GRB_DoubleAttr_X) > 0.9) {
            result_triangulation.emplace_back(edges[i]);
        }
    }

    return result_triangulation;
}

template<class T>
std::pair<std::vector<std::vector<Edge>>, std::vector<std::vector<Edge>>>
FlipIlpEdgeBased<T>::get_solution() {
    std::pair<std::vector<std::vector<Edge>>, std::vector<std::vector<Edge>>> result_pair;


    for (int i = 0; i < layers; i++) {
        result_pair.second.emplace_back(get_d_vars_one_layer(i));
        result_pair.first.emplace_back(get_e_vars_one_layer(i));
    }

    return result_pair;
}




#endif //HIGHERORDERDELAUNAY_FLIPILPEDGEBASED_HPP