
#include <utility>

#include "AStarFlipDistance/TriangulationHandler.hpp"


TriangulationHandler::TriangulationHandler(const std::vector<Point_2> &vertices, GeometricObjects objects)
            :_vertices(vertices),_objects(std::move(objects)){
    n=(int)_vertices.size();
    std::vector<std::vector<int>> edge_matrix(n, std::vector<int>(n, -1));
    _edge_matrix=edge_matrix;
    //assert(_objects.get_number_of_triangles()<32767 ); //maxsize int16
    std::vector<std::vector<std::pair<int, int>>> flip_matrix(n, std::vector<std::pair<int, int>>(n, {-1, -1}));
    _current_flips_matrix=flip_matrix;
}


void TriangulationHandler::init_triangulation(std::vector<triangle_int> &triangles) {
    std::vector<Edge> edge_representation;
    if(!std::is_sorted(triangles.begin(),triangles.end())){
        std::sort(triangles.begin(),triangles.end());
    }
    if(!_triangle_representation.empty()){
        reset_triangulation();
    }

    for(auto ind:triangles){
        _short_representation.emplace_back((triangle_int)ind);
        _triangle_representation.emplace_back(_objects.get_triangle(ind));
    }

    for(auto t:_triangle_representation){
        int v1=t.vertices[0];
        int v2=t.vertices[1];
        int v3=t.vertices[2];

        _edge_matrix[v1][v2]=1;
        if(_current_flips_matrix[v1][v2].first<0){
            edge_representation.emplace_back(v1,v2);
            _current_flips_matrix[v1][v2].first=v3;
        }
        else{
            if(_current_flips_matrix[v1][v2].second>0){
                std::cerr <<"there where three flipable edges";
            }
            else{
                //_edge_representation.emplace_back(v1,v2);
                _current_flips_matrix[v1][v2].second=v3;
            }
        }

        _edge_matrix[v1][v3]=1;
        if(_current_flips_matrix[v1][v3].first<0){
            edge_representation.emplace_back(v1,v3);
            _current_flips_matrix[v1][v3].first=v2;
        }
        else{
            if(_current_flips_matrix[v1][v3].second>0){
                std::cerr <<"there where three flipable edges";
            }
            else{
                //_edge_representation.emplace_back(v1,v3);
                _current_flips_matrix[v1][v3].second=v2;
            }
        }

        _edge_matrix[v2][v3]=1;
        if(_current_flips_matrix[v2][v3].first<0){
            edge_representation.emplace_back(v2,v3);
            _current_flips_matrix[v2][v3].first=v1;
        }
        else{
            if(_current_flips_matrix[v2][v3].second>0){
                std::cerr <<"there where three flipable edges";
            }
            else{
                //_edge_representation.emplace_back(v2,v3);
                _current_flips_matrix[v2][v3].second=v1;
            }
        }
    }
    std::vector<int> int_edge_rep;
    for(auto e:edge_representation){
        if(_current_flips_matrix[e.src][e.dst].second!=-1){
            int_edge_rep.emplace_back(_objects.get_edge(e.src,e.dst));
        }
    }
    std::sort(int_edge_rep.begin(),int_edge_rep.end());
    _interior_edge_representation=int_edge_rep;

}

void TriangulationHandler::init_triangulation(const std::vector<Triangle>& triangles) {
    std::vector<triangle_int> triangles_ind;
    for(auto t:triangles){
        triangles_ind.emplace_back(_objects.get_triangle(t));
    }
    init_triangulation(triangles_ind);
}

void TriangulationHandler::reset_triangulation() {
    for (int i = 0; i < _current_flips_matrix.size(); ++i) {
        for (int j = 0; j < _current_flips_matrix[i].size(); ++j) {
            _current_flips_matrix[i][j] = {-1, -1};
        }
    }
    for (int i = 0; i < _edge_matrix.size(); ++i) {
        for (int j = 0; j < _edge_matrix[i].size(); ++j) {
            _edge_matrix[i][j] = -1;
        }
    }
    _interior_edge_representation.clear();
    _triangle_representation.clear();
    _short_representation.clear();

}

void TriangulationHandler::print_stuff() {
    for(auto e:_interior_edge_representation){
        _objects.get_edge(e).print();
    }
    std::cout<<std::endl;
    for(auto e:_triangle_representation){
        e.print();
    }
    std::cout<<std::endl;
    for(auto e:_short_representation){
        std::cout<<e<<" ";
    }
    std::cout<<std::endl;

    for(int i=0;i<_edge_matrix.size();i++){
        for(int j=0;j<_edge_matrix.size();j++){
            std::cout<<_edge_matrix[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    for(int i=0;i<_current_flips_matrix.size();i++){
        for(int j=0;j<_current_flips_matrix.size();j++){
            std::cout<<"["<<_current_flips_matrix[i][j].first<<", "<<_current_flips_matrix[i][j].second<<"]";
        }
        std::cout<<std::endl;
    }

}

void TriangulationHandler::reset_relative_triangulation() {
    for(auto i :_interior_edge_representation){
        auto e=_objects.get_edge(i);
        _edge_matrix[e.src][e.dst]=-1;
        _current_flips_matrix[e.src][e.dst]={-1,-1};
    }
    _interior_edge_representation.clear();
    _triangle_representation.clear();
    _short_representation.clear();

}

Edge TriangulationHandler::do_flip_for_real(Edge e) {
    int e_0=e.src;
    int e_1=e.dst;


    int f_0=std::min(_current_flips_matrix[e_0][e_1].first,_current_flips_matrix[e_0][e_1].second);
    int f_1=std::max(_current_flips_matrix[e_0][e_1].first,_current_flips_matrix[e_0][e_1].second);
    if(f_0<0 || f_1 <0 ){
        return {-1,-1};
    }
    if(!_objects.is_flip_valid(e_0,e_1,f_0,f_1)){
        return {-1,-1};
    }

    _edge_matrix[e_0][e_1]=-1;
    _edge_matrix[f_0][f_1]=1;

    for(auto & i : _interior_edge_representation){
        if(i==_objects.get_edge(e.src,e.dst)){
            i=_objects.get_edge(f_0,f_1);
        }
    }


    _current_flips_matrix[e_0][e_1]={-1,-1};
    _current_flips_matrix[f_0][f_1]={e_0,e_1};





    //change the flips for adjacent edges
    if(_current_flips_matrix[std::min(e_0,f_0)][std::max(e_0,f_0)].first==e_1){
        _current_flips_matrix[std::min(e_0,f_0)][std::max(e_0,f_0)].first=f_1;
    }
    else{
        assert(_current_flips_matrix[std::min(e_0,f_0)][std::max(e_0,f_0)].second==e_1);
        _current_flips_matrix[std::min(e_0,f_0)][std::max(e_0,f_0)].second=f_1;
    }

    if(_current_flips_matrix[std::min(e_1,f_0)][std::max(e_1,f_0)].first==e_0){
        _current_flips_matrix[std::min(e_1,f_0)][std::max(e_1,f_0)].first=f_1;
    }
    else{
        assert(_current_flips_matrix[std::min(e_1,f_0)][std::max(e_1,f_0)].second==e_0);
        _current_flips_matrix[std::min(e_1,f_0)][std::max(e_1,f_0)].second=f_1;
    }

    if(_current_flips_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].first==e_1){
        _current_flips_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].first=f_0;
    }
    else{
        assert(_current_flips_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].second==e_1);
        _current_flips_matrix[std::min(e_0,f_1)][std::max(e_0,f_1)].second=f_0;
    }

    if(_current_flips_matrix[std::min(e_1,f_1)][std::max(e_1,f_1)].first==e_0){
        _current_flips_matrix[std::min(e_1,f_1)][std::max(e_1,f_1)].first=f_0;
    }
    else{
        assert(_current_flips_matrix[std::min(e_1,f_1)][std::max(e_1,f_1)].second==e_0);
        _current_flips_matrix[std::min(e_1,f_1)][std::max(e_1,f_1)].second=f_0;
    }

    //changing triangle representation
    Triangle old_1{e_0,e_1,f_0};
    Triangle old_2{e_0,e_1,f_1};
    int ind_o1=_objects.get_triangle(old_1);
    int ind_o2=_objects.get_triangle(old_2);
    Triangle new_1{f_0,f_1,e_0};
    Triangle new_2{f_0,f_1,e_1};
    int ind_n1=_objects.get_triangle(new_1);
    int ind_n2=_objects.get_triangle(new_2);

    for(int i=0;i<_triangle_representation.size();i++){
        if(_short_representation[i]==ind_o1){
            _short_representation[i]=(triangle_int)ind_n1;
            _triangle_representation[i]=new_1;
        }
        if(_short_representation[i]==ind_o2){
            _short_representation[i]=(triangle_int)ind_n2;
            _triangle_representation[i]=new_2;
        }
    }
    std::sort(_short_representation.begin(),_short_representation.end());
    std::sort(_triangle_representation.begin(),_triangle_representation.end());
    std::sort(_interior_edge_representation.begin(),_interior_edge_representation.end());



    return {f_0,f_1};

}

int TriangulationHandler::get_simple_heuristic(const StaticTriangulation &target_triangulation) {
    assert(_edge_matrix[0].size()==target_triangulation._edge_matrix[0].size());
    assert(_edge_matrix[0].size()==target_triangulation._edge_matrix[0].size());
    int counter=0;
    for (int i=0;i<_edge_matrix.size();i++){
        for(int j=0;j<target_triangulation._edge_matrix.size();j++){
            if(_edge_matrix[i][j]>target_triangulation._edge_matrix[i][j]){
                counter++;
            }
        }
    }

    return counter;
}

std::vector<triangle_int> TriangulationHandler::get_representation() {
    std::vector<triangle_int> result(_short_representation);
    return result;
}

bool TriangulationHandler::equals_triangulation(const StaticTriangulation &target_triangulation) {
    return std::equal(_short_representation.begin(), _short_representation.end(),
                      target_triangulation._short_representation.begin());

}

std::vector<int> TriangulationHandler::get_edge_representation() {
    return _interior_edge_representation;
}

std::vector<Triangle> TriangulationHandler::get_triangle_representation() {
    return _triangle_representation;
}

Edge TriangulationHandler::do_flip_implicit(Edge e) {
    int e_0=e.src;
    int e_1=e.dst;


    int f_0=std::min(_current_flips_matrix[e_0][e_1].first,_current_flips_matrix[e_0][e_1].second);
    int f_1=std::max(_current_flips_matrix[e_0][e_1].first,_current_flips_matrix[e_0][e_1].second);
    if(f_0<0 || f_1 <0 ){
        return {-1,-1};
    }
    if(!_objects.is_flip_valid(e_0,e_1,f_0,f_1)){
        return {-1,-1};
    }

    _edge_matrix[e_0][e_1]=-1;
    _edge_matrix[f_0][f_1]=1;


    _current_flips_matrix[e_0][e_1]={-1,-1};
    _current_flips_matrix[f_0][f_1]={e_0,e_1};





    //changing triangle representation
//    Triangle old_1{e_0,e_1,f_0};
//    Triangle old_2{e_0,e_1,f_1};
    int ind_o1=_objects.get_triangle(e_0,e_1,f_0);
    int ind_o2=_objects.get_triangle(e_0,e_1,f_1);
//    Triangle new_1{f_0,f_1,e_0};
//    Triangle new_2{f_0,f_1,e_1};
    int ind_n1=_objects.get_triangle(f_0,f_1,e_0);
    int ind_n2=_objects.get_triangle(f_0,f_1,e_1);


    int pos= binarySearch(_short_representation,ind_o1);
    _short_representation[pos]=(triangle_int)ind_n1;
    bubble_to_correct_spot(_short_representation,pos);

    pos= binarySearch(_short_representation,ind_o2);
    _short_representation[pos]=(triangle_int)ind_n2;
    bubble_to_correct_spot(_short_representation,pos);


    return {f_0,f_1};
}

int TriangulationHandler::binarySearch(const std::vector<triangle_int> &vec, int val) {
    int low = 0;
    int high = vec.size() - 1;

    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (vec[mid] == val) {
            return mid;  // Found the element
        } else if (vec[mid] < val) {
            low = mid + 1;  // Search in the right half
        } else {
            high = mid - 1;  // Search in the left half
        }
    }

    // Element not found, return the position where it should be inserted
    return low;
}

void TriangulationHandler::bubble_to_correct_spot(std::vector<triangle_int> &vec, int pos) {
    int nn = vec.size();
    bool isSwapped = true;

    while (pos > 0 && isSwapped) {
        isSwapped = false;

        if (vec[pos] < vec[pos - 1]) {
            std::swap(vec[pos], vec[pos - 1]);
            --pos;
            isSwapped = true;
        }
    }

    isSwapped = true;
    while (pos < nn - 1 && isSwapped) {
        isSwapped = false;

        if (vec[pos] > vec[pos + 1]) {
            std::swap(vec[pos], vec[pos + 1]);
            ++pos;
            isSwapped = true;
        }
    }
}

std::pair<std::vector<triangle_int>,int> TriangulationHandler::do_flip_implicit_for_simple_heuristic(Edge e,const StaticTriangulation& target){
    int offset=0;

    int e_0=e.src;
    int e_1=e.dst;

/**/    int f_0=std::min(_current_flips_matrix[e_0][e_1].first,_current_flips_matrix[e_0][e_1].second);
    int f_1=std::max(_current_flips_matrix[e_0][e_1].first,_current_flips_matrix[e_0][e_1].second);
    if(f_0<0 || f_1 <0 ){
        return {{},-1};
    }
    if(!_objects.is_flip_valid(e_0,e_1,f_0,f_1)){
        return {{},-1};
    }


    //changing triangle representation
    int ind_o1=_objects.get_triangle(e_0,e_1,f_0);
    int ind_o2=_objects.get_triangle(e_0,e_1,f_1);
    int ind_n1=_objects.get_triangle(f_0,f_1,e_0);
    int ind_n2=_objects.get_triangle(f_0,f_1,e_1);

    std::vector<triangle_int> result(_short_representation);


    int pos= binarySearch(result,ind_o1);
    result[pos]=(triangle_int)ind_n1;
    bubble_to_correct_spot(result,pos);

    pos= binarySearch(result,ind_o2);
    result[pos]=(triangle_int)ind_n2;
    bubble_to_correct_spot(result,pos);

    int e_ind=_objects.get_edge(e_0,e_1);
    int f_ind=_objects.get_edge(f_0,f_1);

    if(target._edges_index_in_triangulation[e_ind]&&target._edges_index_in_triangulation[f_ind]){
        offset=0;
    }
    if(!target._edges_index_in_triangulation[e_ind]&&!target._edges_index_in_triangulation[f_ind]){
        offset=0;
    }
    if(target._edges_index_in_triangulation[e_ind]&&!target._edges_index_in_triangulation[f_ind]){
        offset=1;
    }
    if(!target._edges_index_in_triangulation[e_ind]&&target._edges_index_in_triangulation[f_ind]){
        offset=-1;
    }


    return {result,offset};
}


/*int TriangulationHandler::get_eppstein_heuristic_naive(const StaticTriangulation &target_triangulation) {



    auto matrix=get_distance_matrix_for_eppstein(target_triangulation);

    HungarianPrimDual hungi(matrix);
    auto result=hungi.execute();
    int val=0;
    for(int i=0;i<result.size();i++){
        val+=matrix[i][result[i]];
    }
    //std::cout<<val<<std::endl;
    return val;


}*/

int TriangulationHandler::get_eppstein_heuristic_parent(const StaticTriangulation &target_triangulation) {
    auto matrix=get_distance_matrix_for_eppstein(target_triangulation);
    auto edges= get_relevant_edges_for_eppstein(target_triangulation);
    _dynamic_hungarian_matrix=HungarianPrimDual(matrix);
    auto result=_dynamic_hungarian_matrix.parent_execute(edges.first,_objects.get_number_of_edges());
    return result;
}

int TriangulationHandler::get_eppstein_heuristic_child(const StaticTriangulation &target_triangulation,
                                                       int flipped_parent_edge, int child_edge) {


    auto edges_target= get_relevant_edges_for_eppstein(target_triangulation).second;
    std::vector<triangle_int> new_row_costs;
    new_row_costs.reserve(edges_target.size());
    for(int i=0;i<edges_target.size();i++){
        new_row_costs.emplace_back(_objects.get_distance(child_edge,edges_target[i]));
    }
    int result=_dynamic_hungarian_matrix.child_execute(flipped_parent_edge,new_row_costs);
    return result;
}


std::pair<std::vector<int>, std::vector<int>>
TriangulationHandler::get_relevant_edges_for_eppstein(const StaticTriangulation &target_triangulation) {
    return {_interior_edge_representation,target_triangulation._interior_edge_representation};
}


std::vector<std::vector<triangle_int>>
TriangulationHandler::get_distance_matrix_for_eppstein(const StaticTriangulation &target_triangulation) {
    auto x=get_relevant_edges_for_eppstein(target_triangulation);

    int imimax=_objects.get_number_of_edges()*_objects.get_number_of_edges();
    std::vector<std::vector<triangle_int>> matrix;
    matrix.resize(x.first.size(), std::vector<triangle_int>(x.second.size()));
    for(int i=0;i<x.first.size();i++){
        for(int j=0;j<x.first.size();j++){
            matrix[i][j]=_objects.get_distance(x.first[i],x.second[j]);
        }
    }

    return matrix;
}



int getRandomInteger(int min, int max) {
    return min + (std::rand() % (max - min + 1));
}

std::vector<Triangle> TriangulationHandler::do_n_valid_flips(int n, int seed) {
    srand(seed);
    std::vector<int> already_considered;
    int counter=0;
    while(counter<n){
        int ra= getRandomInteger(0,_interior_edge_representation.size()-1);

        auto e=_objects.get_edge(_interior_edge_representation[ra]);
        auto x= do_flip_for_real(e);
        if(x.src==-1){
            continue;
        }
        else{
            already_considered.emplace_back(_interior_edge_representation[ra]);
            counter++;
        }

    }
    return _triangle_representation;
}

std::vector<int> TriangulationHandler::get_intersection_difference_count(const StaticTriangulation &target) {
    std::vector<int> result;
    for (int i=0;i<_interior_edge_representation.size();i++){
        int counter=0;
        auto e = _objects.get_edge(_interior_edge_representation[i]);
        for (int j=0;j<target._interior_edge_representation.size();j++){

            auto f= _objects.get_edge(target._interior_edge_representation[j]);
            if(e.src==f.src||e.dst==f.src|| e.src==f.dst || e.dst==f.dst){
                continue;
            }
            if(do_arcs_intersect(_vertices[e.src],_vertices[e.dst],_vertices[f.src],_vertices[f.dst])){
                counter++;
            }
        }
        e=do_flip_for_real(e);
        if(e.src==-1){
            counter=-10000000;
        }
        for (int j=0;j<target._interior_edge_representation.size();j++){

            auto f= _objects.get_edge(target._interior_edge_representation[j]);
            if(e.src==f.src||e.dst==f.src|| e.src==f.dst || e.dst==f.dst){
                continue;
            }
            if(do_arcs_intersect(_vertices[e.src],_vertices[e.dst],_vertices[f.src],_vertices[f.dst])){
                counter--;
            }
        }
        if(counter<-1000){
            result.emplace_back(counter);
        }
        else{
            do_flip_for_real(e);
            result.emplace_back(counter);
        }


    }
    return result;
}

int TriangulationHandler::get_heuristic_upper_bound(const StaticTriangulation & target){
    int flip=0;
    while(true){
        auto x= get_intersection_difference_count(target);
        int ind=-1;
        int max=-1;
        for(int i=0;i<x.size();i++){
            if(x[i]>max){
                max=x[i];
                ind=i;
            }
        }
        do_flip_for_real(_objects.get_edge(_interior_edge_representation[ind]));
        flip++;
        if(equals_triangulation(target)){
            return flip;
        }
    }



    return 0;
}

Edge TriangulationHandler::do_random_flip() {
    int nr_edges=(int)_interior_edge_representation.size();
    int edge_int=std::rand() % (nr_edges );
    Edge e=_objects.get_edge(_interior_edge_representation[edge_int]);
    auto x=do_flip_for_real(e);
    return x;
}

int TriangulationHandler::get__simple_heuristic_offset(Edge e, Edge f, StaticTriangulation &target) {
    int e_ind=_objects.get_edge(e.src,e.dst);
    int f_ind=_objects.get_edge(f.src,f.dst);
    int offset=-50;
    if(target._edges_index_in_triangulation[e_ind]&&target._edges_index_in_triangulation[f_ind]){
        offset=0;
    }
    if(!target._edges_index_in_triangulation[e_ind]&&!target._edges_index_in_triangulation[f_ind]){
        offset=0;
    }
    if(target._edges_index_in_triangulation[e_ind]&&!target._edges_index_in_triangulation[f_ind]){
        offset=1;
    }
    if(!target._edges_index_in_triangulation[e_ind]&&target._edges_index_in_triangulation[f_ind]){
        offset=-1;
    }
    return offset;
}

std::vector<Triangle> TriangulationHandler::do_n_random_flips(int nn, int seed) {
    std::srand(seed);
    for(int i=0;i<nn;i++){
        do_random_flip();
    }

    return get_triangle_representation();
}







