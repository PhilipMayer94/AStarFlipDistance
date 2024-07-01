
#include <utility>
#include "AStarFlipDistance/BiDirectionalBfs.hpp"

BiDirectionalBfs::BiDirectionalBfs(std::vector<Point_2> points, std::vector<Triangle> start_triangulation,
                                   std::vector<Triangle> target_triangulation) : _objects(points,false){
    _start=StaticTriangulation(_objects,start_triangulation);
    _target=StaticTriangulation(_objects,target_triangulation);
    _triangulation_from_start=TriangulationHandler(points,_objects);
    _triangulation_from_target=TriangulationHandler(points,_objects);
    _triangulation_from_start.init_triangulation(start_triangulation);
    _triangulation_from_target.init_triangulation(target_triangulation);
}

void BiDirectionalBfs::run() {
    auto t1 = std::chrono::high_resolution_clock::now();
    _closed=0;
    if(_triangulation_from_start.equals_triangulation(_target)){
        _flipdistance=0;
    }
    _queue_from_start.emplace(_triangulation_from_start.get_representation(),0);
    _queue_from_target.emplace(_triangulation_from_target.get_representation(),0);
    _closed_from_start.emplace(_triangulation_from_start.get_representation(),0);
    _closed_from_target.emplace(_triangulation_from_target.get_representation(),0);
    bool do_start=true;
    int d=0;
    while(true){
        _closed++;

        auto t_timeout = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::high_resolution_clock ::now() - t1;
        if(duration > std::chrono::minutes(_timeout_in_minutes)){
            _flipdistance=-1;
            _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            return;
        }


        if(do_start){


            //get queue element
            auto cur=_queue_from_start.front();
            if(cur.second>d){
                do_start=false;
                continue;
            }
            int parent_distance=cur.second;
            _queue_from_start.pop();
            _triangulation_from_start.init_triangulation(cur.first);

            //check if its in the other front
            auto x=_closed_from_target.find(cur.first);
            if(x!=_closed_from_target.end()){
                _flipdistance= parent_distance+x->second;
                auto t2 = std::chrono::high_resolution_clock::now();
                _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
                return;
            }
            //_closed_from_start.emplace(cur);
            //extend it
            auto possible_flips=_triangulation_from_start.get_edge_representation();
            for(auto e:possible_flips){

                auto f=_triangulation_from_start.do_flip_implicit(_triangulation_from_start._objects.get_edge(e));
                if(f.src==-1){
                    continue;
                }

                int distance=parent_distance+1;
                auto x_2=_closed_from_start.find(_triangulation_from_start.get_representation());
                if((x_2==_closed_from_start.end())){
                    _queue_from_start.push(std::pair<std::vector<triangle_int>,int>(_triangulation_from_start.get_representation(),distance));
                    _closed_from_start.emplace(_triangulation_from_start.get_representation(),distance);
                }


                _triangulation_from_start.do_flip_implicit(f);

            }
        }
        else{
            if(_queue_from_target.empty()){
                _flipdistance=-22;
                return ;
            }


            //get queue element
            auto cur=_queue_from_target.front();
            if(cur.second>d){
                do_start=true;
                d++;
                continue;
            }
            int parent_distance=cur.second;
            _queue_from_target.pop();
            _triangulation_from_target.init_triangulation(cur.first);

            //check if its in the other front
            auto x=_closed_from_start.find(cur.first);
            if(x!=_closed_from_start.end()){
                _flipdistance= parent_distance+x->second;
                auto t2 = std::chrono::high_resolution_clock::now();
                _runtime=(int)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
                return;
            }
            //_closed_from_target.emplace(cur);

            //extend it
            auto possible_flips=_triangulation_from_target.get_edge_representation();
            for(auto e:possible_flips){

                auto f=_triangulation_from_target.do_flip_implicit(_triangulation_from_target._objects.get_edge(e));
                if(f.src==-1){
                    continue;
                }

                int distance=parent_distance+1;
                auto x_2=_closed_from_target.find(_triangulation_from_target.get_representation());
                if(x_2==_closed_from_target.end()){
                    _queue_from_target.push(std::pair<std::vector<triangle_int>,int>(_triangulation_from_target.get_representation(),distance));
                    _closed_from_target.emplace(_triangulation_from_target.get_representation(),distance);
                }
                _triangulation_from_target.do_flip_implicit(f);
            }
        }
    }
}

std::string BiDirectionalBfs::get_data_as_string() {
    std::string result=std::to_string(_flipdistance)+" "+std::to_string(_runtime)+" "+std::to_string(_closed)+" "
                       + std::to_string(_closed_from_target.size()+_closed_from_start.size());
    return result;
}
