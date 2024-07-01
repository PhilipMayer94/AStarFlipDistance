

#include "AStarFlipDistance/FastEppsteinHeuristic.hpp"

FastEppsteinHeuristic::FastEppsteinHeuristic(std::vector<Point_2> &vertices,
                                             const std::vector<Triangle> &start_triangulation,
                                             std::vector<Triangle> &target_triangulation, GeometricObjects &obj) {
    _vertices=vertices;
    _objects=obj;
    _start_triangulation=TriangulationHandler(vertices,_objects);
    _start_triangulation.init_triangulation(start_triangulation);
    _target_triangulation= StaticTriangulation(_objects,target_triangulation);
}

FastEppsteinHeuristic::FastEppsteinHeuristic(std::vector<Point_2> &vertices,
                                             const std::vector<Triangle> &start_triangulation,
                                             std::vector<Triangle> &target_triangulation) {
    _vertices=vertices;
    _objects=GeometricObjects(vertices,false);
    _start_triangulation=TriangulationHandler(vertices,_objects);
    _start_triangulation.init_triangulation(start_triangulation);
    _target_triangulation= StaticTriangulation(_objects,target_triangulation);
}

void FastEppsteinHeuristic::run() {
    int n=_objects.get_number_of_edges();
    std::vector<triangle_int> _current_distances;

    _current_distances.resize(n,-1);
    std::vector<std::vector<triangle_int>> cost_matrix(_start_triangulation.get_edge_representation().size(),
                                                       std::vector<triangle_int>(_start_triangulation.get_edge_representation().size()));



    auto edges_start=_start_triangulation.get_edge_representation();
    auto edges_target=_target_triangulation._interior_edge_representation;

    std::vector<std::vector<int>> neighbors(n);

    for(int i=0;i<_objects._allowed_flips.size();i++){
        for(int j=i+1;j<_objects._allowed_flips[i].size();j++){
            if(_objects._allowed_flips[i][j]){
                neighbors[i].emplace_back(j);
                neighbors[j].emplace_back(i);
            }
        }
    }



    std::vector<int> visited(n,-1);
    std::queue<std::pair<int,int>> queue=std::queue<std::pair<int,int>>();;
    for(int i=0;i<edges_start.size();i++){
        auto e =edges_start[i];
        queue.emplace(e,0);
        visited[e]=i;
        while(!queue.empty()){
            auto current=queue.front();
            queue.pop();
            _current_distances[current.first]=static_cast<triangle_int>(current.second);
            for(auto next:neighbors[current.first]){
                if(visited[next]!=i){
                    queue.emplace(next,current.second+1);
                    visited[next]=i;
                }
            }
        }

        for(int k=0;k<edges_target.size();k++){
            auto f =edges_target[k];
            if(visited[f]==i){
                cost_matrix[i][k]=_current_distances[f];
            }
            else{
                cost_matrix[i][k]=10000;
            }
        }
    }
    HungarianPrimDual hungarian(cost_matrix);
    _heuristic_value= hungarian.parent_execute(edges_target,_objects.get_number_of_edges());
    std::cout<<"heuristic: " <<_heuristic_value<<std::endl;
}
