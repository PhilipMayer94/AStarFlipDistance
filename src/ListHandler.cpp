
#include "AStarFlipDistance/ListHandler.hpp"




int ListHandler::Parent(int i) { return (i - 1) / 2; }

int ListHandler::LeftChild(int i) { return 2 * i + 1; }

int ListHandler::RightChild(int i) { return 2 * i + 2; }


void ListHandler::insert(const std::vector<triangle_int>& triangulation, int p_new_estimated_Distance ,int shortest_path_distance, int simple_heuristic){
    _instance_counter++;
    int last_item=_the_items.size()-1;
    int i = _heap.size();

    _the_position_in_the_heap.emplace_back(i);
    _is_in_open_list.emplace_back(true);
    _heap.emplace_back(last_item, NodeInformation{p_new_estimated_Distance,shortest_path_distance,_instance_counter,simple_heuristic,last_item});

    _open[last_item] = 1;

    while (i > 0 and std::get<1>(_heap[Parent(i)]) > std::get<1>(_heap[i])){
        std::swap(_heap[i], _heap[Parent(i)]);
        _the_position_in_the_heap[std::get<1>(_heap[i])._position_in_vector] = i;
        _the_position_in_the_heap[std::get<1>(_heap[Parent(i)])._position_in_vector] = Parent(i);
        i = Parent(i);
    }
}




void ListHandler::decrease_Key(const std::vector<triangle_int>&  triangulation, int p_new_estimated_Distance ,int shortest_path_distance, int simple_heuristic , const std::unordered_map<int,int,HasherInt,EqualFnInt> ::const_iterator& found_element ){
    auto x= found_element;
    int triang_ind=x->first;
    int i=_the_position_in_the_heap[triang_ind];

    _instance_counter++;
    if(_heap[i].second._shortest_path_distance<=shortest_path_distance){
        return;
    }
    _heap[i] =  std::pair(triang_ind, NodeInformation{p_new_estimated_Distance,shortest_path_distance,_instance_counter,simple_heuristic,triang_ind});

    while (i > 0 and std::get<1>(_heap[Parent(i)]) > std::get<1>(_heap[i])){
        std::swap(_heap[i], _heap[Parent(i)]);
        _the_position_in_the_heap[std::get<1>(_heap[i])._position_in_vector] = i;
        _the_position_in_the_heap[std::get<1>(_heap[Parent(i)])._position_in_vector] = Parent(i);
        i = Parent(i);
    }
}



bool ListHandler::is_Empty(){
    if(_heap.size() > 0){
        return false;
    }else{
        return true;
    }
}

std::pair<std::vector<triangle_int>, NodeInformation> ListHandler::extract_min(){
    if(_heap.size() > 1) {
        auto min = _heap[0];
        _heap[0] = _heap.back();
        _heap.pop_back();
        _is_in_open_list[min.first]=false;//optimal
        _the_position_in_the_heap[std::get<1>(min)._position_in_vector]=-1;
        _the_position_in_the_heap[std::get<1>(_heap[0])._position_in_vector]=0;
        int i = 0;
        while (LeftChild(i) < _heap.size()) {
            int j = LeftChild(i);

            if (j + 1 < _heap.size() && std::get<1>(_heap[j + 1]) < std::get<1>(_heap[j])) {
                ++j;
            }

            if (std::get<1>(_heap[i]) > std::get<1>(_heap[j])) {
                std::swap(_heap[i], _heap[j]);
                _the_position_in_the_heap[std::get<1>(_heap[i])._position_in_vector] = i;
                _the_position_in_the_heap[std::get<1>(_heap[j])._position_in_vector] = j;
                i = j;
            } else {
                break;
            }
        }
        std::pair<std::vector<triangle_int>, NodeInformation> result={_the_items[min.first],min.second};
        return result;
    }else{
        auto min = _heap[0];
        _heap[0] = _heap.back();
        _heap.pop_back();
        _is_in_open_list[min.first]=false;
        _the_position_in_the_heap[std::get<1>(min)._position_in_vector]=-1;
        std::pair<std::vector<triangle_int>, NodeInformation> result={_the_items[min.first],min.second};
        return result;
    }
}

bool ListHandler::insert_or_decrease_key_if_not_in_closed(const std::vector<triangle_int>& triangulation, int p_new_estimated_Distance ,int shortest_path_distance,int simple_heuristic) {
    int last_item=_the_items.size();
    _the_items.emplace_back(triangulation);


    if(_open.empty()){
        insert(triangulation,p_new_estimated_Distance,shortest_path_distance,simple_heuristic);
        return true;
    }
    auto possible_map_object=_open.find(last_item);

    if(_open.empty()||possible_map_object==_open.end()){
        insert(triangulation,p_new_estimated_Distance,shortest_path_distance,simple_heuristic);
        return true;
    }

    if(_is_in_open_list[possible_map_object->first]==0){
        _the_items.pop_back();
        return false;
    }

    else{
        _the_items.pop_back();
        decrease_Key(triangulation,p_new_estimated_Distance,shortest_path_distance,simple_heuristic,possible_map_object);
    }
    return true;
}


bool ListHandler::insert_or_decrease_key_if_not_in_closed(const std::vector<triangle_int>& triangulation, int p_new_estimated_Distance ,int shortest_path_distance) {
    return insert_or_decrease_key_if_not_in_closed(triangulation, p_new_estimated_Distance ,shortest_path_distance,-42);
}

int ListHandler::get_nr_of_extended_nodes() {
    int counter=0;
    for(auto iter :_is_in_open_list){
        if(!(iter)){
            counter++;
        }
    }
    return counter;
}

int ListHandler::get_nr_of_opened_nodes() {
    int counter=0;
    for(auto iter :_is_in_open_list){
        if((iter)){
            counter++;
        }
    }
    return counter;
}

void ListHandler::set_hasher() {
    HasherInt::set_item_vector(_the_items);
    EqualFnInt::set_item_vector(_the_items);
}

