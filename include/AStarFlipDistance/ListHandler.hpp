
#ifndef A_STAR_FOR_FLIPDISTANCE_LISTHANDLER_HPP
#define A_STAR_FOR_FLIPDISTANCE_LISTHANDLER_HPP
#include <unordered_map>
#include <vector>
#include <iostream>
#include "Global.hpp"

class Hasher
{
public:

    std::size_t operator() (std::vector<triangle_int> const& key) const
        {
            std::size_t seed = key.size();

            // Combine each element and its position to compute the hash
            for (std::size_t i = 0; i < key.size(); ++i) {
                std::hash<triangle_int> hasher;
                seed ^= hasher(key[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2) + i;
            }

            return seed;
        }

};
class EqualFn {
public:
    bool operator()(std::vector<triangle_int> const &t1, std::vector<triangle_int> const &t2) const {
        return std::equal(t1.begin(), t1.end(), t2.begin());
    }
};



class HasherInt
{

public:
    inline static std::vector<std::vector<triangle_int>> * _items= {};
    static void set_item_vector(std::vector<std::vector<triangle_int>> & items){
        _items=&items;
    }

    std::size_t operator() (int const& key_int) const
    {
        auto items =_items;
        const std::vector<triangle_int>& key=(*items)[key_int];

        std::size_t seed = key.size();

        // Combine each element and its position to compute the hash
        for (std::size_t i = 0; i < key.size(); ++i) {
            std::hash<triangle_int> hasher;
            seed ^= hasher(key[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2) + i;
        }

        return seed;
    }


};

class EqualFnInt {
public:
    inline static std::vector<std::vector<triangle_int>> * _items= {};
    static void set_item_vector(std::vector<std::vector<triangle_int>> & items){
        _items=&items;
    }
    bool operator()(int const &t1_ind, int const &t2_ind) const {
        auto items =_items;
        const std::vector<triangle_int>& t1=(*items)[t1_ind];
        const std::vector<triangle_int>& t2=(*items)[t2_ind];
        return std::equal(t1.begin(), t1.end(), t2.begin());
    }

};


struct NodeInformation{
    int32_t _heuristic;
    int32_t _shortest_path_distance;
    int32_t _number;
    int32_t _simple_heuristic=-10000;


    int32_t _position_in_vector=-1;



    NodeInformation(int heuristic, int distance,int number,int simple_heuristic, int position_in_vector){
        _heuristic=heuristic;
        _shortest_path_distance=distance;
        _number=number;
        _simple_heuristic=simple_heuristic;
        _position_in_vector=position_in_vector;
    }

    bool operator<(const NodeInformation& other) const {
        if(_heuristic < other._heuristic){
            return true;
        }
        else if (_heuristic > other._heuristic){
            return false;
        }
        else {
            if(_shortest_path_distance > other._shortest_path_distance)
            {
                return true;
            }
            else if (_shortest_path_distance <other._shortest_path_distance) {
                return false;
            }
            else {
                return _number>other._number;

            }
        }
    }

    bool operator>(const NodeInformation& other) const {
        if (_heuristic > other._heuristic) {
            return true;
        } else if (_heuristic < other._heuristic) {
            return false;
        } else {
            return _shortest_path_distance < other._shortest_path_distance;
        }
    }
};

class ListHandler {
public:
    ListHandler() {
        _the_items={};
        HasherInt::set_item_vector(_the_items);
        EqualFnInt::set_item_vector(_the_items);
    }
    std::pair<std::vector<triangle_int>,NodeInformation> extract_min();
    bool insert_or_decrease_key_if_not_in_closed(const std::vector<triangle_int>& triangulation, int p_new_estimated_Distance ,int shortest_path_distance);
    bool insert_or_decrease_key_if_not_in_closed(const std::vector<triangle_int>& triangulation, int p_new_estimated_Distance ,int shortest_path_distance,int simple_heuristic);
    int get_nr_of_extended_nodes();
    int get_nr_of_opened_nodes();

    void set_hasher();


    std::vector<std::vector<triangle_int>> _the_items;
private:

    std::vector<int> _the_position_in_the_heap;
    std::vector<bool> _is_in_open_list;

    std::unordered_map<int,int,HasherInt,EqualFnInt> _open;
    std::vector<std::pair<int,NodeInformation>> _heap;



    void insert(const std::vector<triangle_int>& triangulation, int p_estimated_Distance,int shortest_path_distance, int simple_heuristic);
    void decrease_Key(const std::vector<triangle_int>&  triangulation, int p_new_estimated_Distance ,int shortest_path_distance, int simple_heuristic , const std::unordered_map<int,int,HasherInt,EqualFnInt>::const_iterator& found_element );

    bool is_Empty();
    static int Parent(int i);
    static int LeftChild(int i);
    static int RightChild(int i);

    int _instance_counter=0;


};




#endif //A_STAR_FOR_FLIPDISTANCE_LISTHANDLER_HPP
