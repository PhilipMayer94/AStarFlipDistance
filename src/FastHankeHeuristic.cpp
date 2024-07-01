
#include <utility>

#include "AStarFlipDistance/FastHankeHeuristic.hpp"


FastHankeHeuristic::FastHankeHeuristic(std::vector<Point_2> &vertices, const std::vector<Triangle> &start_triangulation,
                                       std::vector<Triangle> &target_triangulation)  {


    auto obj=GeometricObjects(vertices,false);
    TriangulationHandler handler(vertices,obj);
    handler.init_triangulation(start_triangulation);

    auto interior_diagonals=handler.get_edge_representation();
    auto flip_matrix=handler._current_flips_matrix;
    _target= StaticTriangulation(obj,target_triangulation);;
    _heap=MaxHeap( interior_diagonals.size());


    for(auto i : interior_diagonals ){
        _current_diagonals.emplace_back(i);
    }


    _objects=obj;
    _nr_of_diagonals=interior_diagonals.size();
    _flip_matrix=std::vector<std::vector<std::pair<int,int>>> (flip_matrix.size(), std::vector<std::pair<int, int>>(flip_matrix.size()));
    _where_in_diagonal_list=std::vector<std::vector<int>> (flip_matrix.size(),std::vector<int>(flip_matrix.size(),-1));
    for(int i=0;i<flip_matrix.size();i++){
        for(int j=0;j<flip_matrix.size();j++){
            _flip_matrix[i][j]=std::pair<int,int>(flip_matrix[i][j].first,flip_matrix[i][j].second);
        }
    }

    initalize_datastructures();
}


FastHankeHeuristic::FastHankeHeuristic(GeometricObjects &objects, StaticTriangulation &target,
                                       std::vector<int> interior_diagonals,
                                       std::vector<std::vector<std::pair<int, int>>> flip_matrix) : _heap(interior_diagonals.size()){
   _heap=MaxHeap( interior_diagonals.size());
    _objects=objects;
    _target=target;
    for(auto i : interior_diagonals ){
        _current_diagonals.emplace_back(i);
    }
    _nr_of_diagonals=interior_diagonals.size();
    _flip_matrix=std::vector<std::vector<std::pair<int,int>>> (flip_matrix.size(), std::vector<std::pair<int, int>>(flip_matrix.size()));
    _where_in_diagonal_list=std::vector<std::vector<int>> (flip_matrix.size(),std::vector<int>(flip_matrix.size(),-1));

    for(int i=0;i<flip_matrix.size();i++){
        for(int j=0;j<flip_matrix.size();j++){
            _flip_matrix[i][j]=std::pair<int,int>(flip_matrix[i][j].first,flip_matrix[i][j].second);
        }
    }

    initalize_datastructures();


}

void FastHankeHeuristic::do_flip(int diagonal_index){
    int e_ind=_current_diagonals[diagonal_index];
    Edge e=_objects.get_edge(e_ind);
    int e_0=e.src;
    int e_1=e.dst;

    Edge f(_flip_matrix[e_0][e_1].first,_flip_matrix[e_0][e_1].second);
    int f_0=f.src;
    int f_1=f.dst;

    if(f_0<0 || f_1 <0 ){
        std::cout<< "flip_was_not_valid_1  ";
        return;
    }

    int f_ind=_objects.get_edge(f_0,f_1);

    if(e==Edge(0, 6)){
        int debug=123;
    }

    if(!_objects.is_flip_valid(e_0,e_1,f_0,f_1)){
        std::cout<< "flip_was_not_valid_2zz  ";
        return;
    }

    if(_target._edges_index_in_triangulation[f_ind]){
        different_edge_counter--;
    }
    if(different_edge_counter==0){
        _flip_distance++;
        return;
    }

    _where_in_diagonal_list[e_0][e_1]=-5;
    _where_in_diagonal_list[f_0][f_1]=diagonal_index;

    _current_diagonals[diagonal_index]=f_ind;
    _flip_matrix[e_0][e_1]={-1,-1};
    _flip_matrix[f_0][f_1]={e_0,e_1};



    std::vector<Edge> additional_edges_that_need_to_be_changed;

    //0 0
    if (_flip_matrix[std::min(e_0, f_0)][std::max(e_0, f_0)].first == e_1) {
        _flip_matrix[std::min(e_0, f_0)][std::max(e_0, f_0)].first = f_1;
        Edge tmp = Edge(e_0, f_0);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }

    } else {
        assert(_flip_matrix[std::min(e_0, f_0)][std::max(e_0, f_0)].second == e_1);
        _flip_matrix[std::min(e_0, f_0)][std::max(e_0, f_0)].second = f_1;
        Edge tmp = Edge(e_0, f_0);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }
    }

    //1 0
    if (_flip_matrix[std::min(e_1, f_0)][std::max(e_1, f_0)].first == e_0) {
        _flip_matrix[std::min(e_1, f_0)][std::max(e_1, f_0)].first = f_1;
        Edge tmp = Edge(e_1, f_0);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }
    } else {
        assert(_flip_matrix[std::min(e_1, f_0)][std::max(e_1, f_0)].second == e_0);
        _flip_matrix[std::min(e_1, f_0)][std::max(e_1, f_0)].second = f_1;
        Edge tmp = Edge(e_1, f_0);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }

    }

    //0 1
    if(_flip_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].first==e_1){
        _flip_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].first=f_0;
        Edge tmp=Edge(e_0,f_1);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }
    }
    else{
        assert(_flip_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].second==e_1);
        _flip_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].second=f_0;
        Edge tmp=Edge(e_0,f_1);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }
    }

    //1 1
    if (_flip_matrix[std::min(e_1, f_1)][std::max(e_1, f_1)].first == e_0) {
        _flip_matrix[std::min(e_1, f_1)][std::max(e_1, f_1)].first = f_0;
        Edge tmp = Edge(e_1, f_1);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }
    } else {
        assert(_flip_matrix[std::min(e_1, f_1)][std::max(e_1, f_1)].second == e_0);
        _flip_matrix[std::min(e_1, f_1)][std::max(e_1, f_1)].second = f_0;
        Edge tmp = Edge(e_1, f_1);
        if(tmp.src>=0){
            additional_edges_that_need_to_be_changed.emplace_back(tmp);
        }
    }

    for(auto g :additional_edges_that_need_to_be_changed){



        int h_1 = std::min(_flip_matrix[g.src][g.dst].first, _flip_matrix[g.src][g.dst].second);
        int h_2 = std::max(_flip_matrix[g.src][g.dst].first, _flip_matrix[g.src][g.dst].second);




        int position = _where_in_diagonal_list[g.src][g.dst];
        int tmp_value=-1;

        if(h_1<0||position<0){
            continue;
        }
        if(_objects.is_flip_valid(g.src,g.dst,h_1,h_2)){
            tmp_value = _number_of_intersections_with_target[_objects.get_edge(g.src,g.dst)] -
                        _number_of_intersections_with_target[_objects.get_edge(h_1, h_2)];
        }else{
            tmp_value=-1000;
        }

        if(tmp_value!=_current_diagonal_weights[position]){
            _current_diagonal_weights[position]=tmp_value;
            _heap.changeKey(position,tmp_value);
        }

    }

    int new_value=_number_of_intersections_with_target[f_ind]-_number_of_intersections_with_target[_objects.get_edge(_flip_matrix[f_0][f_1].first,_flip_matrix[f_0][f_1].second)];
    assert(new_value==(_number_of_intersections_with_target[f_ind]-_number_of_intersections_with_target[e_ind]));
    if(_objects.is_flip_valid(f_0,f_1,e_0,e_1)){
        _current_diagonal_weights[diagonal_index]=new_value;
    }else{
        _current_diagonal_weights[diagonal_index]=-1000;
    }

    _heap.insert(diagonal_index,new_value);
    _flip_distance++;

}



int FastHankeHeuristic::compute_heuristic() {
    while(different_edge_counter>0){
        auto edge= _heap.extractMax();
        do_flip(edge.second);
    }
    return _flip_distance;
}

void FastHankeHeuristic::set_new_start_triangulation(const std::vector<int>& interior_diagonals,
                                                     std::vector<std::vector<std::pair<int, int>>> &flip_matrix) {
    _flip_distance=0;
    _heap=MaxHeap(interior_diagonals.size());


    for(auto i : _current_diagonals){
        auto x=_objects.get_edge(i);
        _where_in_diagonal_list[x.src][x.dst]=-1;
    }

    _current_diagonals.clear();
    _current_diagonal_weights.clear();
    for(auto i : interior_diagonals ){
        _current_diagonals.emplace_back(i);
    }

    for(int i=0;i<flip_matrix.size();i++){
        for(int j=0;j<flip_matrix.size();j++){
            _flip_matrix[i][j]=std::pair<int,int>(flip_matrix[i][j].first,flip_matrix[i][j].second);
        }
    }
    int pos=0;
    for(auto i : _current_diagonals){
        auto x=_objects.get_edge(i);
        _where_in_diagonal_list[x.src][x.dst]=pos;
        pos++;
    }


    int i=0;
    for(int _current_diagonal : _current_diagonals){
        auto e=_objects.get_edge(_current_diagonal);
        int current_value;
        if(_objects.is_flip_valid(e.src,e.dst,_flip_matrix[e.src][e.dst].first,_flip_matrix[e.src][e.dst].second)){
            current_value=_number_of_intersections_with_target[_current_diagonal]-_number_of_intersections_with_target[_objects.get_edge(_flip_matrix[e.src][e.dst].first,_flip_matrix[e.src][e.dst].second)];
        }else{
            current_value=-1000;
        }
        _current_diagonal_weights.emplace_back(current_value);
        _heap.insert(i,_current_diagonal_weights[i]);
        i++;
    }

    different_edge_counter=0;
    for(auto d :_current_diagonals){
        if(!_target._edges_index_in_triangulation[d]){
            different_edge_counter++;
        }
    }

}

int FastHankeHeuristic::compute_heuristic_for_triangulation(std::vector<int> interior_diagonals,
                                                            std::vector<std::vector<std::pair<int, int>>> &flip_matrix) {
    set_new_start_triangulation(interior_diagonals,flip_matrix);
    compute_heuristic();
    return _flip_distance;
}

FastHankeHeuristic::FastHankeHeuristic(std::vector<Point_2> &vertices, const std::vector<Triangle> &start_triangulation,
                                       std::vector<Triangle> &target_triangulation, GeometricObjects &objects) {

    _objects=objects;

    TriangulationHandler handler(vertices,_objects);
    handler.init_triangulation(start_triangulation);

    auto interior_diagonals=handler.get_edge_representation();
    auto flip_matrix=handler._current_flips_matrix;
    _target= StaticTriangulation(_objects,target_triangulation);;
    _heap=MaxHeap( interior_diagonals.size());


    for(auto i : interior_diagonals ){
        _current_diagonals.emplace_back(i);
    }


    _nr_of_diagonals=interior_diagonals.size();
    _flip_matrix=std::vector<std::vector<std::pair<int,int>>> (flip_matrix.size(), std::vector<std::pair<int, int>>(flip_matrix.size()));
    _where_in_diagonal_list=std::vector<std::vector<int>> (flip_matrix.size(),std::vector<int>(flip_matrix.size(),-1));
    for(int i=0;i<flip_matrix.size();i++){
        for(int j=0;j<flip_matrix.size();j++){
            _flip_matrix[i][j]=std::pair<int,int>(flip_matrix[i][j].first,flip_matrix[i][j].second);
        }
    }

    initalize_datastructures();
}

void FastHankeHeuristic::initalize_datastructures() {
    int pos=0;
    for(auto i : _current_diagonals){
        auto x=_objects.get_edge(i);
        _where_in_diagonal_list[x.src][x.dst]=pos;
        pos++;
    }


    _vertices=_objects.get_vertices();
    for(int i=0;i<_objects.get_number_of_edges();i++){
        int counter=0;
        for(auto j :_target._interior_edge_representation){
            if(i==j){
                continue;
            }
            auto e1=_objects.get_edge(i);
            auto e2=_objects.get_edge(j);
            if(e1.src==e2.src||e1.dst==e2.dst||e1.src==e2.dst||e1.dst==e2.src){
                continue;
            }
            if(do_arcs_intersect(_vertices[e1.src],_vertices[e1.dst],_vertices[e2.src],_vertices[e2.dst])){
                counter++;
            }

        }
        _number_of_intersections_with_target.emplace_back(counter);
    }
    int i=0;
    for(int _current_diagonal : _current_diagonals){
        auto e=_objects.get_edge(_current_diagonal);
        int current_value;
        if(_objects.is_flip_valid(e.src,e.dst,_flip_matrix[e.src][e.dst].first,_flip_matrix[e.src][e.dst].second)){
            current_value=_number_of_intersections_with_target[_current_diagonal]-_number_of_intersections_with_target[_objects.get_edge(_flip_matrix[e.src][e.dst].first,_flip_matrix[e.src][e.dst].second)];
        }else{
            current_value=-1000;
        }
        _current_diagonal_weights.emplace_back(current_value);
        _heap.insert(i,_current_diagonal_weights[i]);
        i++;
    }

    different_edge_counter=0;
    for(auto d :_current_diagonals){
        if(!_target._edges_index_in_triangulation[d]){
            different_edge_counter++;
        }
    }

}

