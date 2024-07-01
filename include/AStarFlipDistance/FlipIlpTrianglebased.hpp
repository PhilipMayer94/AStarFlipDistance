
#ifndef HIGHERORDERDELAUNAY_FLIPILPTRIANGLEBASED_HPP
#define HIGHERORDERDELAUNAY_FLIPILPTRIANGLEBASED_HPP

#include <gurobi_c++.h>
#include "BasicDataStructures.hpp"
#include "HigherOrderDelaunayObjects.hpp"
#include "AStarFlipDistance/FlipIlp.hpp"
#include "AStarFlipDistance/GeometricTriangulation.hpp"



template<class T>
class FlipIlpTriangleBased : FlipIlp<T> {
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
    FlipIlpTriangleBased(GeometricTriangulation<T> triangulation_1, GeometricTriangulation<T> triangulation_2,
                         HigherOrderDelaunayObjects<T> _objects);


    void solve();

    double get_number_of_flips();

    void print_solution();

    std::pair<std::vector<std::vector<Triangle>>, std::vector<std::vector<Triangle>>> get_solution();

    std::vector<Triangle> get_t_vars_one_layer(int layer);

    std::vector<Triangle> get_d_vars_one_layer(int layer);


private:

    void create_model();


    std::vector<int> get_triangulation_indices(GeometricTriangulation<T> triangulation);

    std::vector<std::vector<GRBVar>> t_vars;
    std::vector<std::vector<GRBVar>> a_vars;
    std::vector<std::vector<GRBVar>> d_vars;


};

template<class T>
FlipIlpTriangleBased<T>::FlipIlpTriangleBased(GeometricTriangulation<T> triangulation_1,
                                              GeometricTriangulation<T> triangulation_2,
                                              HigherOrderDelaunayObjects<T> _objects) :FlipIlp<T>(triangulation_1,
                                                                                                  triangulation_2,
                                                                                                  _objects) {
    create_model();
}

template<class T>
void FlipIlpTriangleBased<T>::create_model() {

    GRBEnv env;
    env.set(GRB_IntParam_OutputFlag, 0);
    model = new GRBModel(env);
    model->set(GRB_DoubleParam_TimeLimit, 180.0);
    for (int i = 0; i < this->layers + 1; i++) {
        t_vars.emplace_back(std::vector<GRBVar>{});
        a_vars.emplace_back(std::vector<GRBVar>{});
        d_vars.emplace_back(std::vector<GRBVar>{});
        for (int j = 0; j < triangles.size(); j++) {
            std::string str = "T_";
            str += std::to_string(i) + ". " + std::to_string(j);

            t_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY, str));

            str = "A_";
            str += std::to_string(i) + ". " + std::to_string(j);
            a_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY, str));
            str = "D_";
            str += std::to_string(i) + ". " + std::to_string(j);
            d_vars[i].emplace_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY, str));
        }
    }
    for (auto t_vars_on_layer: t_vars) {
        int layer = -1;
        layer++;
        for (int i = 0; i < edges.size(); i++) {
            GRBLinExpr expr;
            if (!left_triangles_adjacent_to_edge[i].empty() && !right_triangles_adjacent_to_edge[i].empty()) {
                for (auto t: left_triangles_adjacent_to_edge[i]) {
                    expr = expr + t_vars_on_layer[t];
                }
                for (auto t: right_triangles_adjacent_to_edge[i]) {
                    expr = expr + (-1) * t_vars_on_layer[t];
                }
                std::string str = "C_" + std::to_string(layer) + "__" + std::to_string(edges[i].src) + "_" +
                                  std::to_string(edges[i].dst);
                model->addConstr(expr, GRB_EQUAL, 0, str);
            } else if (!left_triangles_adjacent_to_edge[i].empty()) {
                for (auto t: left_triangles_adjacent_to_edge[i]) {
                    expr += t_vars_on_layer[t];
                }
                std::string str = "C_" + std::to_string(layer) + "__" + std::to_string(edges[i].src) + "_" +
                                  std::to_string(edges[i].dst);
                model->addConstr(expr, GRB_EQUAL, 1, str);
            } else if (!right_triangles_adjacent_to_edge[i].empty()) {
                for (auto t: right_triangles_adjacent_to_edge[i]) {
                    expr += t_vars_on_layer[t];
                }
                std::string str = "C_" + std::to_string(layer) + "__" + std::to_string(edges[i].src) + "_" +
                                  std::to_string(edges[i].dst);
                model->addConstr(expr, GRB_EQUAL, 1, str);
            }
        }

        if (!has_boundary) {
            GRBLinExpr expr;
            for (int i = 0; i < triangles.size(); i++) {
                expr += t_vars_on_layer[i];
            }
            std::string str = "C_+ std::to_string(layer)+__size";
            model->addConstr(expr, GRB_EQUAL, m, str);
        }
    }
    if (!has_boundary) {
        for (int l = 0; l < layers; l++) {
            GRBLinExpr expr;
            for (int i = 0; i < triangles.size(); i++) {
                expr += t_vars[l][i];
            }
            std::string str = "C_size";
            model->addConstr(expr, GRB_EQUAL, m, str);
        }
    }

    //setting start and target triangulation values
    std::vector<int> start_indices = get_triangulation_indices(start_triangulation);
    std::vector<int> target_indices = get_triangulation_indices(target_triangulation);
    std::vector<int> start_values(triangles.size(), 0);
    std::vector<int> target_values(triangles.size(), 0);
    for (int i = 0; i < start_indices.size(); i++) {
        start_values[start_indices[i]] = 1;
        target_values[target_indices[i]] = 1;
    }

    for (int i = 0; i < triangles.size(); i++) {
        GRBLinExpr expr_start;
        GRBLinExpr expr_target;
        expr_start += t_vars[0][i];
        expr_target += t_vars[t_vars.size() - 1][i];
        std::string str = "C_start_" + std::to_string(i);
        model->addConstr(expr_start, GRB_EQUAL, start_values[i], str);
        str = "C_target_" + std::to_string(i);
        model->addConstr(expr_target, GRB_EQUAL, target_values[i], str);
    }

    //d values
    GRBLinExpr expr_objective;
    for (int i = 1; i < t_vars.size(); i++) {
        GRBLinExpr expr_all;
        for (int j = 0; j < triangles.size(); j++) {
            GRBLinExpr expr;
            expr_objective += d_vars[i][j];
            expr_all += d_vars[i][j];
            expr = expr + t_vars[i][j] - t_vars[i - 1][j] - d_vars[i][j];
            std::string str = "C_d_e_" + std::to_string(i) + "__" + std::to_string(j);
            model->addConstr(expr, GRB_LESS_EQUAL, 0, str);
        }
        std::string str = "C_d_sum_" + std::to_string(i);
        model->addConstr(expr_all, GRB_LESS_EQUAL, 2, str);
    }
    model->setObjective(expr_objective, GRB_MINIMIZE);
}


template<class T>
void FlipIlpTriangleBased<T>::solve() {
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
std::vector<int> FlipIlpTriangleBased<T>::get_triangulation_indices(GeometricTriangulation<T> triangulation) {
    std::vector<int> result_indices;
    for (auto t: triangulation.triangles) {
        for (int i = 0; i < triangles.size(); i++) {
            if (t == triangles[i]) {
                result_indices.template emplace_back(i);
            }
        }
    }
    return result_indices;
}

template<class T>
double FlipIlpTriangleBased<T>::get_number_of_flips() {

    return optimization_result / 2;
}

template<class T>
void FlipIlpTriangleBased<T>::print_solution() {
    std::pair<std::vector<std::vector<Triangle>>, std::vector<std::vector<Triangle>>> solution = get_solution();
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
        std::cout << "T_" << i << ": ";
        for (auto t: solution.first[i]) {
            t.print();
        }
        std::cout << std::endl;


    }

}

template<class T>
std::vector<Triangle> FlipIlpTriangleBased<T>::get_t_vars_one_layer(int layer) {
    std::vector<Triangle> result_triangulation;
    for (int i = 0; i < triangles.size(); i++) {
        if (t_vars[layer][i].get(GRB_DoubleAttr_X) > 0.9) {
            result_triangulation.emplace_back(triangles[i]);
        }
    }

    return result_triangulation;
}


template<class T>
std::vector<Triangle> FlipIlpTriangleBased<T>::get_d_vars_one_layer(int layer) {
    std::vector<Triangle> result_triangulation;
    for (int i = 0; i < triangles.size(); i++) {
        if (d_vars[layer][i].get(GRB_DoubleAttr_X) > 0.9) {
            result_triangulation.emplace_back(triangles[i]);
        }
    }

    return result_triangulation;
}

template<class T>
std::pair<std::vector<std::vector<Triangle>>, std::vector<std::vector<Triangle>>>
FlipIlpTriangleBased<T>::get_solution() {
    std::pair<std::vector<std::vector<Triangle>>, std::vector<std::vector<Triangle>>> result_pair;


    for (int i = 0; i < layers; i++) {
        result_pair.second.emplace_back(get_d_vars_one_layer(i));
        result_pair.first.emplace_back(get_t_vars_one_layer(i));
    }

    return result_pair;
}


#endif //HIGHERORDERDELAUNAY_FLIPILPTRIANGLEBASED_HPP