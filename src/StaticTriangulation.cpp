
#include "AStarFlipDistance/StaticTriangulation.hpp"



StaticTriangulation::StaticTriangulation(GeometricObjects objects, std::vector<Triangle> &triangles_real) : _objects(std::move(objects)){
    n=_objects.get_number_of_vertices();

    std::vector<Edge> edge_representation;
    std::vector<int> triangles;
    for(auto t:triangles_real){
        triangles.emplace_back(_objects.get_triangle(t));
    }
    if(!std::is_sorted(triangles.begin(),triangles.end())){
        std::sort(triangles.begin(),triangles.end());
    }

    for(auto ind:triangles){
        _short_representation.emplace_back((short)ind);
        _triangle_representation.emplace_back(_objects.get_triangle(ind));
    }

    std::vector<std::vector<int>> edge_matrix(n, std::vector<int>(n, -1));
    std::vector<std::vector<std::pair<int, int>>> flip_matrix(n, std::vector<std::pair<int, int>>(n, {-1, -1}));
    _edge_matrix=edge_matrix;
    _current_flips_matrix=flip_matrix;

    std::vector<bool> boolVector(_objects.get_number_of_edges(), false);
    _edges_index_in_triangulation=boolVector;

    for(auto t:_triangle_representation){
        int v1=t.vertices[0];
        int v2=t.vertices[1];
        int v3=t.vertices[2];

        _edge_matrix[v1][v2]=1;
        if(_current_flips_matrix[v1][v2].first<0){
            _edges_index_in_triangulation[_objects.get_edge(v1,v2)]=true;
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
            _edges_index_in_triangulation[_objects.get_edge(v1,v3)]=true;
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
            _edges_index_in_triangulation[_objects.get_edge(v2,v3)]=true;
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