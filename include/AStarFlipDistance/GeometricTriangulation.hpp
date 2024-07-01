#include <utility>
#include <unordered_set>
#include <cfloat>
#include <climits>
#include <filesystem>

#include "BasicDataStructures.hpp"

#ifndef HIGHERORDERDELAUNAY_GEOMETRICTRIANGULATION_HPP
#define HIGHERORDERDELAUNAY_GEOMETRICTRIANGULATION_HPP



template<class T>
class GeometricTriangulation {


public:
    std::vector<T> vertices;
    std::vector<Triangle> triangles;
    std::vector<Edge> edges;
    int boundary_size = INT_MIN;
    double value = -DBL_MAX;

    GeometricTriangulation(std::vector<Triangle> _triangles, std::vector<T> _vertices, double _value) {
        GeometricTriangulation(_triangles, _vertices);
        value = _value;
    }

    GeometricTriangulation() = default;

    GeometricTriangulation(std::vector<Triangle> _triangles, std::vector<T> _vertices) : vertices(_vertices), triangles(
            std::move(_triangles)) {
        std::unordered_set<Edge, edge_hash> tmpEdges;
        for (auto t: triangles) {
            tmpEdges.insert(Edge(t.vertices[0], t.vertices[1]));
            tmpEdges.insert(Edge(t.vertices[0], t.vertices[2]));
            tmpEdges.insert(Edge(t.vertices[1], t.vertices[2]));
        }
        edges = std::vector<Edge>(tmpEdges.begin(), tmpEdges.end());
    }

    void print_triangles() {
        std::cout << "The triangles of the triangulation: " << std::endl;
        for (auto t: triangles) {
            t.print();
        }
        std::cout << std::endl;
    }

    void print_edges() {
        std::cout << "The edges of the triangulation: " << std::endl;
        for (auto t: edges) {
            t.print();
        }
        std::cout << std::endl;
    }

    void generate_edge_csv(std::string file_location){
        std::string name="../data/results/"+file_location+".csv";
        std::ofstream fx(name, std::ofstream::out | std::ofstream::trunc);
        fx.close();
        std::ofstream fw(name, std::ofstream::out );
        fw<< "lon1"<<"," << "lat1"<<","<<  "lon2"<<"," << "lat2"<<" "<<"\n";
        for(int i=0;i<edges.size();i++ ){
            int v1=edges.at(i).src;
            int v2=edges.at(i).dst;
            double lon1=vertices[v1].lon*180/M_PI;
            double lon2=vertices[v2].lon *180/M_PI;

            if(lon1>180){
                lon1=lon1-360;
            }
            if(lon2>180){
                lon2=lon2-360;
            }

            fw<<lon1<<"," << vertices[v1].lat*180/M_PI<<", "<<  lon2<<"," << vertices[v2].lat*180/M_PI<<" "<<"\n";
        }
        fw.close();
    }

    void generate_triangle_csv_sorted(const std::string& foldername, int year, int month, int order){
        std::sort(triangles.begin(), triangles.end(), compareTriangles);



        std::string folderPath="../data/results/flipdistance/"+foldername;
        if (!std::filesystem::exists(folderPath)) {
            if (std::filesystem::create_directory(folderPath)) {
                std::cout << "Folder created successfully.\n";
            } else {
                std::cout << "Failed to create the folder.\n";
            }
        } else {
            std::cout << "Folder already exists.\n";
        }

        std::string name2= "../data/results/flipdistance/vertices.csv";
        std::ofstream outputFile2(name2);

        if (outputFile2.is_open()) {
            for (const auto& vertex : vertices) {
                outputFile2 << vertex.x << ' ' <<vertex.y  << '\n';
            }
            outputFile2.close();
            std::cout << "Sorting and writing to file complete.\n";
        } else {
            std::cout << "Unable to open the file.\n";
        }


        std::string name= "../data/results/flipdistance/"+foldername+"/t_"+std::to_string(year)+"_"+std::to_string(month)+"_"+std::to_string(order)+".csv";
        std::ofstream outputFile(name);

        if (outputFile.is_open()) {
            for (const auto& triangle : triangles) {
                outputFile << triangle.vertices[0] << ' ' << triangle.vertices[1] << ' ' << triangle.vertices[2] << '\n';
            }
            outputFile.close();
            std::cout << "Sorting and writing to file complete.\n";
        } else {
            std::cout << "Unable to open the file.\n";
        }
    }

    static bool compareTriangles(const Triangle& t1, const Triangle& t2) {
        if (t1.vertices[0] != t2.vertices[0])
            return t1.vertices[0] < t2.vertices[0];
        else if (t1.vertices[1] != t2.vertices[1])
            return t1.vertices[1] < t2.vertices[1];
        else
            return t1.vertices[2] < t2.vertices[2];
    }



};


#endif //HIGHERORDERDELAUNAY_GEOMETRICTRIANGULATION_HPP