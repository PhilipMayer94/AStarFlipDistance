
#include "AStarFlipDistance/AStarFlipDistance.hpp"
#include "AStarFlipDistance/HeuristicDistanceCalculator.hpp"
#include "AStarFlipDistance/FastHankeHeuristic.hpp"

AStarFlipDistance::AStarFlipDistance( std::vector<Point_2> vertices, std::vector<Triangle> start_triangulation,
                                     std::vector<Triangle> target_triangulation) : _triangulation(vertices,GeometricObjects(vertices)),
                                                                                          _target_triangulation(GeometricObjects(vertices),target_triangulation) {
    _triangulation.init_triangulation(start_triangulation);
}




void AStarFlipDistance::run_combined() {
    _list.set_hasher();

    if(_triangulation.get_representation().size()==1){
        _runtime=0;
        _flip_distance=0;
        _has_run=true;
        return;
    }

    FastHankeHeuristic hanke(_target_triangulation._objects,_target_triangulation,_triangulation.get_edge_representation(),_triangulation._current_flips_matrix);
    int hanke_value=hanke.compute_heuristic();

    auto t1 = std::chrono::high_resolution_clock::now();

    //initial values for start_triangulation
    int heuristic_distance=_triangulation.get_eppstein_heuristic_parent(_target_triangulation);
    int simple_heuristic_root=_triangulation.get_simple_heuristic(_target_triangulation);
    int shortest_path_distance=0;
    int estimated_distance=heuristic_distance;



    _initial_heuristic=heuristic_distance;
    if(heuristic_distance==0){
        _flip_distance=0;
        _has_run=true;
        _runtime=0;
        return;
    }

    if(heuristic_distance==hanke_value){
        _flip_distance=hanke_value;
        _has_run=true;
        _runtime=0;
        return;
    }

    if(heuristic_distance==simple_heuristic_root){
        run_simple();
        _has_run=true;
        return;
    }


    _list.insert_or_decrease_key_if_not_in_closed(_triangulation.get_representation(),estimated_distance,shortest_path_distance,simple_heuristic_root);

    while(true){
        bool use_eppstein=true;

        auto t_timeout = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::high_resolution_clock ::now() - t1;
        if(duration > std::chrono::minutes(_timeout_in_minutes)){
            std::cout<<" timed out->";
            _flip_distance=-_list.extract_min().second._heuristic;
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            _has_run=true;
            return;
        }


        //avoid memory overflow on our servers
        if(_list._the_items.size()>900000000){
            std::cout<<" timed out->";
            _flip_distance=-_list.extract_min().second._heuristic;
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            _has_run=true;
            return;
        }


        auto current=_list.extract_min();
        _triangulation.init_triangulation(current.first);
        shortest_path_distance=current.second._shortest_path_distance;
        auto possible_flips=_triangulation.get_edge_representation();



        if(current.second._heuristic-current.second._shortest_path_distance==current.second._simple_heuristic){
            use_eppstein=false;
            counter_simple++;
        }
        else{
            counter_eppstein++;
        }
        int parent_h=-22;
        if(use_eppstein){
            parent_h=_triangulation.get_eppstein_heuristic_parent(_target_triangulation);
        }
        else{
            parent_h=current.second._simple_heuristic;
        }


        if(_triangulation.equals_triangulation(_target_triangulation)){
            _flip_distance=shortest_path_distance;
            auto t2 = std::chrono::high_resolution_clock::now();
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
            std::cout<<"simple used: "<< counter_simple<< " eppsteinused: " <<counter_eppstein  <<std::endl;
            return;
        }

        for(auto e:possible_flips){

            auto f=_triangulation.do_flip_implicit(_triangulation._objects.get_edge(e));
            if(f.src<0){
                continue;
            }
            int offset_simple=_triangulation.get__simple_heuristic_offset(_triangulation._objects.get_edge(e),f,_target_triangulation);
            int simple_heuristic=current.second._simple_heuristic+offset_simple;

            if(use_eppstein){
                heuristic_distance=_triangulation.get_eppstein_heuristic_child(_target_triangulation,e,_triangulation._objects.get_edge(f.src,f.dst));
            }
            else{
                heuristic_distance=simple_heuristic;
            }


            estimated_distance=heuristic_distance+(shortest_path_distance+1);
            _list.insert_or_decrease_key_if_not_in_closed(_triangulation.get_representation(),estimated_distance,shortest_path_distance+1,simple_heuristic);

            _triangulation.do_flip_implicit(f);
        }
    }
}



void AStarFlipDistance::run_epptein() {
    FastHankeHeuristic hanke(_target_triangulation._objects,_target_triangulation,_triangulation.get_edge_representation(),_triangulation._current_flips_matrix);
    int hanke_value=hanke.compute_heuristic();

    auto t1 = std::chrono::high_resolution_clock::now();

    bool found=false;
    //initial values for start_triangulation
    int heuristic_distance=_triangulation.get_eppstein_heuristic_parent(_target_triangulation);
    _initial_heuristic=heuristic_distance;
    int shortest_path_distance=0;
    int estimated_distance=heuristic_distance;
    if(heuristic_distance==0){
        _flip_distance=0;
        _has_run=true;
        _runtime=0;
        return;
    }
    if(heuristic_distance==hanke_value){
        _flip_distance=hanke_value;
        _has_run=true;
        _runtime=0;
        return;
    }
    _list.insert_or_decrease_key_if_not_in_closed(_triangulation.get_representation(),estimated_distance,shortest_path_distance);

    while(true){


        auto t_timeout = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::high_resolution_clock ::now() - t1;
        if(duration > std::chrono::minutes(_timeout_in_minutes)){
            _flip_distance=-_list.extract_min().second._heuristic;
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            _has_run=true;
            return;
        }

        //avoid memory overflow on our servers
        if(_list._the_items.size()>900000000){
            std::cout<<" timed out->";
            _flip_distance=-_list.extract_min().second._heuristic;
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            _has_run=true;
            return;
        }
        auto current=_list.extract_min();
        _triangulation.init_triangulation(current.first);
        shortest_path_distance=current.second._shortest_path_distance;
        auto possible_flips=_triangulation.get_edge_representation();
        int parent_h=_triangulation.get_eppstein_heuristic_parent(_target_triangulation);


        if(_triangulation.equals_triangulation(_target_triangulation)){
            _flip_distance=shortest_path_distance;
            auto t2 = std::chrono::high_resolution_clock::now();
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
            return;
        }

        for(auto e:possible_flips){

            auto f=_triangulation.do_flip_implicit(_triangulation._objects.get_edge(e));
            if(f.src<0){
                continue;
            }



            heuristic_distance=_triangulation.get_eppstein_heuristic_child(_target_triangulation,e,_triangulation._objects.get_edge(f.src,f.dst));

            int x=(parent_h-heuristic_distance);
            if (abs(x)>1 ) {
                throw std::invalid_argument( "not consistent " +std::to_string(abs(x))+" "+ std::to_string(heuristic_distance)+ " " + std::to_string(parent_h));
            }
            estimated_distance=heuristic_distance+(shortest_path_distance+1);
            _list.insert_or_decrease_key_if_not_in_closed(_triangulation.get_representation(),estimated_distance,shortest_path_distance+1);
            _triangulation.do_flip_implicit(f);
        }
    }
}


void AStarFlipDistance::run_simple() {
    FastHankeHeuristic hanke(_target_triangulation._objects,_target_triangulation,_triangulation.get_edge_representation(),_triangulation._current_flips_matrix);
    int hanke_value=hanke.compute_heuristic();
    auto t1 = std::chrono::high_resolution_clock::now();


    bool found=false;
    //initial values for start_triangulation
    int heuristic_distance=_triangulation.get_simple_heuristic(_target_triangulation);
    _initial_heuristic=heuristic_distance;
    int shortest_path_distance=0;
    int estimated_distance=heuristic_distance;

    if(heuristic_distance==0){
        _flip_distance=0;
        _has_run=true;
        _runtime=0;
        return;
    }
    if(heuristic_distance==hanke_value){
        _flip_distance=hanke_value;
        _has_run=true;
        _runtime=0;
        return;
    }


    _list.insert_or_decrease_key_if_not_in_closed(_triangulation.get_representation(),estimated_distance,shortest_path_distance);
    while(!found){

        auto duration = std::chrono::high_resolution_clock ::now() - t1;
        if(duration > std::chrono::minutes(_timeout_in_minutes)){
            std::cout<<" timed out";
            auto x=_list.extract_min();
            _flip_distance=-x.second._heuristic;
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            _has_run=true;
            return;
        }

        //avoid memory overflow on our servers
        if(_list._the_items.size()>900000000){
            std::cout<<" timed out->";
            _flip_distance=-_list.extract_min().second._heuristic;
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            _has_run=true;
            return;
        }

        auto current=_list.extract_min();
        _triangulation.init_triangulation(current.first);
        shortest_path_distance=current.second._shortest_path_distance;

        int parent_heuristic=current.second._heuristic-shortest_path_distance;
        //int parent_heuristic=current.second.first-shortest_path_distance;

        auto possible_flips=_triangulation.get_edge_representation();

        if(_triangulation.equals_triangulation(_target_triangulation)){
            _flip_distance=shortest_path_distance;
            auto t2 = std::chrono::high_resolution_clock::now();
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
            return;
        }

        for(auto e:possible_flips){
            auto f=_triangulation.do_flip_implicit_for_simple_heuristic(_triangulation._objects.get_edge(e),_target_triangulation);
            if(f.first.empty()){
                continue;
            }



            int delta_heuristic=f.second;
            //int delta_heuristic=f.second;;

            estimated_distance=parent_heuristic+delta_heuristic+(shortest_path_distance+1);
            //estimated_distance=parent_heuristic+delta_heuristic+(shortest_path_distance+1);
            _list.insert_or_decrease_key_if_not_in_closed(f.first,estimated_distance,shortest_path_distance+1);
        }
    }
}


int AStarFlipDistance::get_flipdistance() {
    return _flip_distance;
}

int AStarFlipDistance::get_runtime() {
    return _runtime;
}

bool AStarFlipDistance::are_triangulation_representations_equal(const std::vector<triangle_int>& t_1, const std::vector<triangle_int>&t_2) {
    for (int i = 0; i < t_1.size(); i++) {
        if (t_2[i] != t_1[i]) {
            return false;
        }
    }

    return true;
}

std::string AStarFlipDistance::get_data_as_string() {
    std::string result=std::to_string(get_flipdistance())+" "+std::to_string( _initial_heuristic)+" "+std::to_string(get_runtime())+" "
            + std::to_string(_list.get_nr_of_extended_nodes())+ " "+ std::to_string(_list.get_nr_of_opened_nodes()); //+
            //+ " "+ std::to_string(counter_eppstein)+ " " + std::to_string(counter_simple);


    return result;
}

void AStarFlipDistance::run_root_combined() {
    FastHankeHeuristic hanke(_target_triangulation._objects,_target_triangulation,_triangulation.get_edge_representation(),_triangulation._current_flips_matrix);
    int hanke_value=hanke.compute_heuristic();

    auto t1 = std::chrono::high_resolution_clock::now();

    bool found=false;
    //initial values for start_triangulation
    int eppstein_distance=_triangulation.get_eppstein_heuristic_parent(_target_triangulation);
    int simple_heuristic_root=_triangulation.get_simple_heuristic(_target_triangulation);
    _initial_heuristic=eppstein_distance;
    int shortest_path_distance=0;
    int estimated_distance=eppstein_distance;
    if(eppstein_distance==0){
        _flip_distance=0;
        _has_run=true;
        _runtime=0;
        return;
    }
    if(eppstein_distance==hanke_value){
        _flip_distance=hanke_value;
        _has_run=true;
        _runtime=0;
        return;
    }
    if(eppstein_distance==simple_heuristic_root){
        run_simple();
    }
    else{
        run_epptein();
    }
    _has_run=true;

}

int AStarFlipDistance::get_initial_heuristic() const {
    return _initial_heuristic;
}

