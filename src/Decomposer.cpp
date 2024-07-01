

#include "AStarFlipDistance/Decomposer.hpp"

Decomposer::Decomposer(const std::vector<Point_2> &points, const std::vector<Edge> &edges) {
    for (auto p : points) {
        nodes.emplace_back(p);
    }
    for (auto e : edges) {
        nodes[e.src].neighbors.emplace_back(e.dst);
        nodes[e.dst].neighbors.emplace_back(e.src);
    }
    _edges=edges;
}

void Decomposer::compute_polygons() {
    bool not_all_processed=true;
    while(not_all_processed) {

        for (int i = 0; i < nodes.size(); i++) {
            auto & start = nodes[i];
            if (start.neighbors.empty()) {
                continue;
            }
            std::vector<int> next_polygon;
            bool go_on = true;
            next_polygon.emplace_back(i);
            next_polygon.emplace_back(start.neighbors.back());
            start.neighbors.pop_back();
            int current = next_polygon[1];
            int pred = next_polygon[0];
            while (go_on) {
                if (nodes[current].neighbors.empty()) {
                    std::cout<<current<<" ";
                    std::cout << "there was no neighbor which should not happen"<<std::endl;
                    break;
                }
                int next_in_list = get_next_neighbor(nodes[pred], nodes[current]);
                int next = nodes[current].neighbors[next_in_list];
                nodes[current].neighbors.erase(nodes[current].neighbors.begin() + next_in_list);

                if (next == i) {
                    go_on = false;
                } else {
                    next_polygon.emplace_back(next);
                    pred = current;
                    current = next;
                }
            }
            polygons.emplace_back(next_polygon);

        }
        not_all_processed=false;
        for(const auto& n:nodes){
            if(n.neighbors.size()>0){
                not_all_processed=true;
            }
        }

    }

    for(int i=polygons.size()-1;i>=0;i--){
        std::reverse(polygons[i].begin(), polygons[i].end());
        double area= shoelaceArea(polygons[i]);
        if(area<0){
            polygons.erase(polygons.begin()+i);
        }
    }


}

void Decomposer::compute_interior_points() {
    for(const auto& polygon : polygons){
        std::vector<int> tmp_interior_points;
        for(int j=0;j<nodes.size();j++){
            auto n=nodes[j];
            if(nodeInPolygon(j,nodes,polygon)){
                tmp_interior_points.emplace_back(j);
            }
        }
        interior_points.emplace_back(tmp_interior_points);
    }
}

void Decomposer::distribute_triangulations(const std::vector<Triangle> &start_triangulation,
                                           const std::vector<Triangle> &target_triangulation) {
    for(auto polygon: polygons){
        std::vector<Triangle> tmp_start;
        std::vector<Triangle> tmp_target;
        for(auto t: start_triangulation){
            int v1=t.vertices[0];
            int v2=t.vertices[1];
            int v3=t.vertices[2];
            if(edge_in_node_polygon(v1,v2,nodes,polygon)&&
               edge_in_node_polygon(v1,v3,nodes,polygon)&&
               edge_in_node_polygon(v2,v3,nodes,polygon)
                    ){
                tmp_start.emplace_back(t);
            }

        }
        for(auto t: target_triangulation){
            int v1=t.vertices[0];
            int v2=t.vertices[1];
            int v3=t.vertices[2];
            if(edge_in_node_polygon(v1,v2,nodes,polygon)&&
               edge_in_node_polygon(v1,v3,nodes,polygon)&&
               edge_in_node_polygon(v2,v3,nodes,polygon)
                    ){
                tmp_target.emplace_back(t);
            }

        }
        start_triangulations.emplace_back(tmp_start);
        target_triangulations.emplace_back(tmp_target);
    }


}

std::vector<Polygon> Decomposer::compute_polygon_decomposition(const std::vector<Triangle> &start_triangulation,
                                                               const std::vector<Triangle> &target_triangulation) {
    std::vector<Polygon> result;


    polygons.clear();
    interior_points.clear();
    start_triangulations.clear();
    target_triangulations.clear();

    compute_polygons();
    compute_interior_points();
    distribute_triangulations(start_triangulation,target_triangulation);

    for(int i=0;i<polygons.size();i++){
        Polygon tmp_polygon;
        std::vector<int> new_ids(nodes.size(), -1);
        int counter=0;
        for(auto p:polygons[i]){
            new_ids[p]=counter;
            tmp_polygon.polygon_walk.emplace_back(counter);
            tmp_polygon.all_points.emplace_back(nodes[p]._point);
            counter++;
        }
        for(auto p:interior_points[i]){
            new_ids[p]=counter;
            tmp_polygon.all_points.emplace_back(nodes[p]._point);
            counter++;
        }
        for(auto t:start_triangulations[i]){
            tmp_polygon.start_triangulation.emplace_back(new_ids[t.vertices[0]],new_ids[t.vertices[1]],new_ids[t.vertices[2]]);
        }
        for(auto t:target_triangulations[i]){
            tmp_polygon.target_triangulation.emplace_back(new_ids[t.vertices[0]],new_ids[t.vertices[1]],new_ids[t.vertices[2]]);
        }

        result.emplace_back(tmp_polygon);
    }




    return result;

}

double Decomposer::shoelaceArea(const std::vector<int> &polygon) {
    int n = polygon.size();
    double area = 0.0;

    for (int i = 0; i < n - 1; ++i) {
        area += nodes[polygon[i]]._point.x * nodes[polygon[i + 1]]._point.y - nodes[polygon[i + 1]]._point.x * nodes[polygon[i]]._point.y;
    }

    area += nodes[polygon[n - 1]]._point.x * nodes[polygon[0]]._point.y - nodes[polygon[0]]._point.x * nodes[polygon[n - 1]]._point.y;

//        area = 0.5 * std::abs(area);
    return area;
}

int Decomposer::get_next_neighbor(const Node &pred, Node current) {
    if(current.neighbors.empty()){
        return -1;
    }
    else{
        double min_angle=1000000;
        int min_ind=-2;
        for(int i=0;i<current.neighbors.size();i++){
            if(calculateAngle(pred._point,current._point,nodes[current.neighbors[i]]._point)<min_angle){
                min_angle=calculateAngle(pred._point,current._point,nodes[current.neighbors[i]]._point);
                min_ind=i;
            }
        }
        return min_ind;
    }

}

bool Decomposer::isLeft(const Point_2 &a, const Point_2 &b, const Point_2 &c) {
    return ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) > 0;
}

bool Decomposer::rayIntersectsEdge(const Point_2 &point, const Point_2 &edgeStart, const Point_2 &edgeEnd) {
    return (edgeStart.y > point.y) != (edgeEnd.y > point.y) &&
           point.x < (edgeEnd.x - edgeStart.x) * (point.y - edgeStart.y) / (edgeEnd.y - edgeStart.y) + edgeStart.x;
}

bool Decomposer::nodeInPolygon(const Point_2 &point, const std::vector<Node> &nodes, const std::vector<int> &polygon) {
    bool inside = false;
    int n = polygon.size();

    for (int i = 0, j = n - 1; i < n; j = i++) {
        Point_2 vi = nodes[polygon[i]]._point;
        Point_2 vj = nodes[polygon[j]]._point;

        if ((vi.y > point.y) != (vj.y > point.y) &&
            point.x < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y) + vi.x) {
            inside = !inside;
        }
    }

    return inside;
}

bool Decomposer::edge_in_node_polygon(int e_1, int e_2, const std::vector<Node> &nodes, const std::vector<int> &polygon) {
    int index_1=-3;
    int index_2=-6;
    for(int i=0;i<polygon.size();i++){
        if(polygon[i]==e_1){
            index_1=i;
        }
        if(polygon[i]==e_2){
            index_2=i;
        }
    }
    if((index_1==0 &&index_2==polygon.size()-1)||(index_2==0 &&index_1==polygon.size()-1)){
        return true;
    }
    if(abs(index_1-index_2)<=1){
        return true;
    }
    auto p1=nodes[e_1];
    auto p2=nodes[e_2];
    Point_2 x(0.5 * (p1._point.x + p2._point.x), 0.5 * (p1._point.y + p2._point.y));
    if(nodeInPolygon(x,nodes,polygon)){
        return true;
    }
    else {
        return false;
    }
}

bool Decomposer::nodeInPolygon(int index, const std::vector<Node> &nodes, const std::vector<int> &polygon) {
    if(std::find(polygon.begin(),polygon.end(),index)!=polygon.end()){
        return false;
    }
    return nodeInPolygon(nodes[index]._point,nodes,polygon);
}
