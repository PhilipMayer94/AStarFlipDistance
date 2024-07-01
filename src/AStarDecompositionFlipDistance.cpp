
#include "AStarFlipDistance/AStarDecompositionFlipDistance.hpp"

void AStarDecompositionFlipDistance::compute_fixed_edges() {
    if(_have_edges_already_been_fixed){
        return;
    }
    _fixed_edges.clear();

    int n=_vertices.size();
    std::vector<std::vector<bool>> matrix_start (n, std::vector<bool>(n, false));
    std::vector<std::vector<bool>> matrix_target(n, std::vector<bool>(n, false));

    for(auto t:_start_triangulation){
        int v0=t.vertices[0];
        int v1=t.vertices[1];
        int v2=t.vertices[2];

        Edge e1(v0,v1);
        Edge e2(v0,v2);
        Edge e3(v2,v1);
        matrix_start[e1.src][e1.dst]=true;
        matrix_start[e2.src][e2.dst]=true;
        matrix_start[e3.src][e3.dst]=true;
    }
    for(auto t:_target_triangulation){
        int v0=t.vertices[0];
        int v1=t.vertices[1];
        int v2=t.vertices[2];

        Edge e1(v0,v1);
        Edge e2(v0,v2);
        Edge e3(v2,v1);
        matrix_target[e1.src][e1.dst]=true;
        matrix_target[e2.src][e2.dst]=true;
        matrix_target[e3.src][e3.dst]=true;
    }

    for (int i=0;i<n;i++) {
        for (int j = i; j < n; j++) {
            if (matrix_start[i][j] && matrix_target[i][j]) {
                _fixed_edges.emplace_back(i, j);
            }
        }
    }
}

void AStarDecompositionFlipDistance::set_fixed_edges(std::vector<Edge> &edges) {
    _fixed_edges=edges;
    _have_edges_already_been_fixed=true;
}

void AStarDecompositionFlipDistance::compute_decomposition() {
    compute_fixed_edges();
    OuterBiConnectedComponentFinder outer_component(_vertices,_fixed_edges);
    outer_component.findTwoConnectedComponents();

    _outer_component=outer_component.find_outer_component_as_edges();
    Decomposer decomposer(_vertices,_outer_component);
    _polygons=decomposer.compute_polygon_decomposition(_start_triangulation,_target_triangulation);
}

void AStarDecompositionFlipDistance::run_with_combined() {
    auto start_time = std::chrono::high_resolution_clock::now();
    compute_decomposition();
    int flipdist=0;

    _sub_problems.reserve(_polygons.size());
    for(auto& p:_polygons) {
        if (p.all_points.size() > 18) {
            FastHankeHeuristic tmp(p.all_points, p.start_triangulation, p.target_triangulation);
            auto x=tmp.compute_heuristic();
            flipdist+=x;
            _used_hanke=true;
            std::cout<<"used hanke heuristic instead of optimal solving, to make sure runtime is small"<<std::endl;
        }
        else{
            AStarFlipDistance tmp(p.all_points, p.start_triangulation, p.target_triangulation);
            _sub_problems.push_back(tmp);
            _sub_problems[_sub_problems.size() - 1].constrain_the_instance_to_a_polygon(p.polygon_walk);
        }
    }

    for(auto & astar:_sub_problems){
        astar.run_combined();
        flipdist+=astar.get_flipdistance();

    }
    _flip_distance=flipdist;
    auto end_time = std::chrono::high_resolution_clock::now();
    _runtime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
}
