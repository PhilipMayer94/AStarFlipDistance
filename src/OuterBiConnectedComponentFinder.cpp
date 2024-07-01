
#include "AStarFlipDistance/OuterBiConnectedComponentFinder.hpp"


void OuterBiConnectedComponentFinder::dfsTwoConnected(int u, int parent) {
    nodes[u].dfsNum = nodes[u].low = dfsCounter++;
    dfsStack.push(u);

    for (int v : nodes[u].neighbors) {
        if (v == parent) continue;

        if (nodes[v].dfsNum == -1) {
            dfsTwoConnected(v, u);
            nodes[u].low = std::min(nodes[u].low, nodes[v].low);

            if (nodes[v].low >= nodes[u].dfsNum) {
                markTwoConnectedComponent(u, v);
            }
        } else {
            nodes[u].low = std::min(nodes[u].low, nodes[v].dfsNum);
        }
    }
}

void OuterBiConnectedComponentFinder::markTwoConnectedComponent(int u, int v) {
    std::vector<int> component;
    while (!dfsStack.empty() && dfsStack.top() != v) {
        nodes[dfsStack.top()].isArticulation = true;
        component.push_back(dfsStack.top());
        dfsStack.pop();
    }

    nodes[u].isArticulation = true;
    component.push_back(u);
    component.push_back(v);
    dfsStack.pop();  // Pop v

    twoConnectedComponents.push_back(component);
}

void OuterBiConnectedComponentFinder::printTwoConnectedComponents() const {
    std::cout << "Two-Connected Components:\n";
    for (const auto& component : twoConnectedComponents) {
        for (int node : component) {
            std::cout << node << " ";
        }
        std::cout << "\n";
    }
}

std::vector<int> OuterBiConnectedComponentFinder::find_outer_component() {
    std::vector<int> result;
    Node2D max({-DBL_MAX,-DBL_MAX});
    Node2D min({DBL_MAX,DBL_MAX});
    int min_ind=-1;
    int max_ind=1;
    for (int i=0;i<nodes.size();i++){
        auto p=nodes[i];
        if(p._point.x<min._point.x) {
            min=p;
            min_ind=i;
        }
        if(p._point.x>max._point.x){
            max=p;
            max_ind=i;
        }
    }
    for(const auto& comp: twoConnectedComponents){
        int counter=0;
        for(auto p:comp){
            result.emplace_back(p);
            if(p==max_ind){
                counter++;
            }
            if(p==min_ind){
                counter++;
            }
        }
        if(counter!=2){
            result.clear();
        }
        else{
            return result;
        }
    }
    std::cout<<"no correct component was found"<<std::endl;
    return {};
}

void OuterBiConnectedComponentFinder::findTwoConnectedComponents() {
    for (int i = 0; i < nodes.size(); ++i) {
        if (nodes[i].dfsNum == -1) {
            dfsTwoConnected(i, -1);
        }
    }
}

OuterBiConnectedComponentFinder::OuterBiConnectedComponentFinder(const std::vector<Point_2> &points,
                                                                 const std::vector<Edge> &edges)
        : dfsCounter(0) {
    for (auto p : points) {
        nodes.emplace_back(p);
    }
    for (auto e : edges) {
        nodes[e.src].neighbors.emplace_back(e.dst);
        nodes[e.dst].neighbors.emplace_back(e.src);
    }
    _edges=edges;
}

std::vector<Edge> OuterBiConnectedComponentFinder::find_outer_component_as_edges() {
    auto verts=find_outer_component();
    std::vector<Edge> result;


    for(auto e:_edges){
        if(std::find(verts.begin(), verts.end(), e.src) != verts.end()
        && std::find(verts.begin(), verts.end(), e.dst) != verts.end()){
            result.emplace_back(e);
        }
    }


    return result;
}
