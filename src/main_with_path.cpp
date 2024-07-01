//
// Created by philip on 4/19/24.
//
#include <iostream>
#include <filesystem>
#include <regex>
#include <random>
#include "AStarFlipDistance/Geometry.hpp"
#include "AStarFlipDistance/BasicDataStructures.hpp"
#include "AStarFlipDistance/GeometricObjects.hpp"
#include "AStarFlipDistance/TriangulationHandler.hpp"
#include "AStarFlipDistance/AStarFlipDistance.hpp"
#include "AStarFlipDistance/HeuristicDistanceCalculator.hpp"
#include "AStarFlipDistance/HigherOrderDelaunayObjects.hpp"
#include "AStarFlipDistance/FlipIlpEdgeBased.hpp"
#include "AStarFlipDistance/BiDirectionalBfs.hpp"
#include "AStarFlipDistance/FastHankeHeuristic.hpp"
#include "AStarFlipDistance/OuterBiConnectedComponentFinder.hpp"
#include "AStarFlipDistance/Decomposer.hpp"
#include "AStarFlipDistance/AStarDecompositionFlipDistance.hpp"
#include "AStarFlipDistance/json.hpp"
#include "AStarFlipDistance/FastEppsteinHeuristic.hpp"
#include "AStarFlipDistance/Global.hpp"


const int SIMPLE=0;
const int EPPSTEIN=1;
const int HEURISTIC=2;
const int BFS=3;
const int ILPEDGE=4;
const int COMBINED=5;
const int ROOTCOMBINED=6;
const int ALLFLIP=7;
const int DECOMPOSITION=8;


std::vector<Point_2> generateRandomPoints(int num_points, double square_size, unsigned int seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> distribution(-square_size / 2.0, square_size / 2.0);

    std::vector<Point_2> points;
    for (int i = 0; i < num_points; ++i) {
        double x = distribution(rng);
        double y = distribution(rng);
        points.emplace_back(x, y);
    }

    return points;
}

void writePointsToCSV(const std::vector<Point_2>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& point : points) {
            file << point.x << ' ' << point.y << '\n';
        }
        file.close();
        std::cout << "CSV file '" << filename << "' written successfully.\n";
    } else {
        std::cerr << "Error opening file: " << filename << "\n";
    }
}

void writeEdgesToCSV(const std::vector<Edge>& edges, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& edge : edges) {
            file << edge.src << ' ' << edge.dst << '\n';
        }
        file.close();
        std::cout << "CSV file '" << filename << "' written successfully.\n";
    } else {
        std::cerr << "Error opening file: " << filename << "\n";
    }
}

std::vector<Triangle> get_triangulation_after_n_random_flips( int n , int seed, std::vector<Point_2>& points,  std::vector<Triangle>& start_triangulation){
    GeometricObjects obj(points,false);
    TriangulationHandler triang(points,obj);
    triang.init_triangulation(start_triangulation);
    std::vector<Triangle> result=triang.do_n_random_flips(n,seed);

    return result;

}

std::vector<std::string> get_file_names_in_folder(std::string foldername){


    std::vector<std::string> fileNames;

    for (const auto& entry : std::filesystem::directory_iterator(foldername)) {
        if (entry.is_regular_file()) {
            fileNames.push_back(entry.path().filename().string());
        }
    }
    std::sort(fileNames.begin(),fileNames.end());
    return fileNames;
}

std::vector<std::pair<std::string,std::string>> get_file_pairs_in_folder(std::string foldername){
    auto names= get_file_names_in_folder(foldername);
    std::vector<std::pair<std::string,std::string>> result;
    for(int i=0;i<names.size();i++){
        for(int j=i+1;j<names.size();j++){
            result.emplace_back(foldername+"/"+names[i],foldername+"/"+names[j]);
        }
    }
    std::sort(result.begin(),result.end());
    return result;
}

std::vector<Point_2> load_vertices(std::string name){
    std::vector<Point_2> points;
    std::ifstream file(name);

    if (!file.is_open()) {
        std::cout << "Failed to open file: " << name << std::endl;
        return points;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x;
        double y;
        if (!(iss >> x >> y)) {
            std::cout << "Invalid line: " << line << std::endl;
            continue;
        }
        points.emplace_back(x,y);
    }

    file.close();
    return points;
}

std::vector<Triangle> load_triangulation(std::string name){
    std::ifstream file(name);
    std::string line;
    std::vector<Triangle> result;

    if (file.is_open()) {
        while (getline(file, line)) {
            std::istringstream iss(line);
            int num1, num2, num3;

            if (iss >> num1 >> num2 >> num3) {
                result.emplace_back(num1,num2,num3);
            } else {
                std::cout << "Invalid line format: " << line << std::endl;
            }
        }
        file.close();
    } else {
        std::cout << "Unable to open file: " << name << std::endl;
    }
    return result;
}


int get_year(const std::string& filename){
    std::regex yearRegex(R"(_(\d{4})_)");

    std::smatch match;
    int year=-1;
    if (std::regex_search(filename, match, yearRegex)) {
        std::string yearString = match.str(1);
        year = std::stoi(yearString);
    }




    return year;
}

std::vector<std::string> do_one_product_of_instances_sea_level(const std::string& folder_path,int type){
    std::string outputfile=folder_path+"/results";
    std::string header;
    std::vector<std::string> result;
    if(type==SIMPLE){
        std::cout<<"doing_simple"<<std::endl;
        outputfile+="_simple.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==EPPSTEIN){
        std::cout<<"doing_eppstein"<<std::endl;
        outputfile+="_eppstein.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }

    else if (type==COMBINED){
        std::cout<<"doing_combined"<<std::endl;
        outputfile+="_combined.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==ROOTCOMBINED){
        std::cout<<"doing_rootcombined"<<std::endl;
        outputfile+="_rootcombined.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==HEURISTIC){
        std::cout<<"doing_heuristic"<<std::endl;
        outputfile+="_heuristic.csv";
        header="file1 file2 simple eppstein hanke\n";
    }

    else if (type==BFS){
        std::cout<<"doing_heuristic"<<std::endl;
        outputfile+="_bfs.csv";
        header="file1 file2 flip runtime closed open\n";
    }
    else if (type==ILPEDGE){
        std::cout<<"doing_ILP"<<std::endl;
        outputfile+="_ilp.csv";
        header="file1 file2 flip runtime primal_bound\n";
    }

    else if (type==DECOMPOSITION){
        std::cout<<"doing_decompositon"<<std::endl;
        outputfile+="_decomposition.csv";
        header="file1 file2 flip runtime used_hanke  \n";
    }

    std::string triangulation_path=folder_path+"/triangulation";
    std::string points_path=folder_path+"/points.csv";
    auto file_pairs= get_file_pairs_in_folder(triangulation_path);
    auto points=load_vertices(points_path);

    std::vector<Triangle> loaded_triangulation_1= load_triangulation(file_pairs[0].first);
    std::vector<Triangle> loaded_triangulation_2= load_triangulation(file_pairs[0].second);
    HeuristicDistanceCalculator dists(points,loaded_triangulation_1,loaded_triangulation_2);

    std::ofstream file(outputfile, std::ios::app);
    std::cout<<"outputfile   "<<outputfile<<std::endl;
    file<< header;

    for(auto p:file_pairs){

        loaded_triangulation_1= load_triangulation(p.first);
        loaded_triangulation_2= load_triangulation(p.second);



        if(type==EPPSTEIN){
            AStarFlipDistance a_star_epp(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_epp.run_epptein();
            std::string to_append=p.first+" "+p.second+" "+a_star_epp.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }
        if(type==COMBINED){
            AStarFlipDistance a_star_comb(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_comb.run_combined();
            std::string to_append=p.first+" "+p.second+" "+a_star_comb.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }
        if(type==ROOTCOMBINED){
            AStarFlipDistance a_star_comb(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_comb.run_root_combined();
            std::string to_append=p.first+" "+p.second+" "+a_star_comb.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }

        else if (type==SIMPLE){
            AStarFlipDistance a_star_simple(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_simple.run_simple();
            std::string to_append=p.first+" "+p.second+" "+a_star_simple.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }
        else if (type==HEURISTIC){
            dists.set_new_start_triangulation(loaded_triangulation_1);
            dists.set_new_target_triangulation(loaded_triangulation_2);
            dists.run();
            std::string to_append=p.first+" "+p.second+" "+std::to_string(dists.simple_heuristic)+" "+std::to_string(dists.eppstein_heuristic)+" "+std::to_string(dists.hanke_heuristic);
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }
        else if (type==BFS){
            dists.set_new_start_triangulation(loaded_triangulation_1);
            dists.set_new_target_triangulation(loaded_triangulation_2);
            dists.run();
            if(dists.eppstein_heuristic==dists.hanke_heuristic){
                std::string to_append=p.first+" "+p.second+" "+std::to_string(dists.eppstein_heuristic)+" "+std::to_string(0)+" "+std::to_string(0)+" "
                                      + std::to_string(0);
                file << to_append<<std::endl;
                result.emplace_back(to_append);
                std::cout<<to_append<<std::endl;
            }else{
                BiDirectionalBfs bfs(points,loaded_triangulation_1,loaded_triangulation_2);
                bfs.run();
                std::string to_append=p.first+" "+p.second+" "+bfs.get_data_as_string();
                file << to_append<<std::endl;
                result.emplace_back(to_append);
                std::cout<<to_append<<std::endl;
            }



        }
        else if (type==ILPEDGE) {

            dists.set_new_start_triangulation(loaded_triangulation_1);
            dists.set_new_target_triangulation(loaded_triangulation_2);
            dists.run();
            if(dists.eppstein_heuristic==dists.hanke_heuristic){
                std::string to_append =
                        p.first + " " + p.second + " " + std::to_string((int) dists.eppstein_heuristic) + " " +
                        std::to_string((int) (0)) + " " + std::to_string(dists.eppstein_heuristic);
                file << to_append<<std::endl;
                result.emplace_back(to_append);
                std::cout<<to_append<<std::endl;
                continue;
            }



            AStarFlipDistance a_star_comb(points, loaded_triangulation_1, loaded_triangulation_2);
            a_star_comb.run_root_combined();

            int real_flip_distance = a_star_comb.get_flipdistance();
            if (real_flip_distance > 0) {
                HigherOrderDelaunayObjects obj(points);
                obj.compute_useful_order_k_edges_and_triangles(points.size());
                FlipIlpEdgeBased flip(GeometricTriangulation(loaded_triangulation_1, points),
                                      GeometricTriangulation(loaded_triangulation_2, points), obj, real_flip_distance);
                flip.solve();

                std::string to_append =
                        p.first + " " + p.second + " " + std::to_string((int) flip.get_number_of_flips()) + " " +
                        std::to_string((int) (flip.runtime * 1000)) + " " + std::to_string((int) flip.get_primal());
                file << to_append << std::endl;
                result.emplace_back(to_append);
                std::cout << to_append << std::endl;
            } else {
                HigherOrderDelaunayObjects obj(points);
                obj.compute_useful_order_k_edges_and_triangles(points.size());

                FlipIlpEdgeBased flip(GeometricTriangulation(loaded_triangulation_1, points),
                                      GeometricTriangulation(loaded_triangulation_2, points), obj,dists.hanke_heuristic);
                flip.solve();

                std::string to_append =
                        p.first + " " + p.second + " " + std::to_string((int) flip.get_number_of_flips()) + " " +
                        std::to_string((int) (flip.runtime * 1000)) + " " + std::to_string((int) flip.get_primal());
                file << to_append << std::endl;
                result.emplace_back(to_append);
                std::cout << to_append << std::endl;
            }
        }

        else if (type==DECOMPOSITION){
            AStarDecompositionFlipDistance a_star_decomposition(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_decomposition.run_with_combined();
            std::string to_append=p.first+" "+p.second+" "+a_star_decomposition.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }

    }
    return result;

}


std::vector<std::string> do_one_product_of_instances(const std::string& folder_path,int type){
    std::string outputfile=folder_path+"/results";
    std::string header;
    std::vector<std::string> result;
    if(type==SIMPLE){
        std::cout<<"doing_simple"<<std::endl;
        outputfile+="_simple.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==EPPSTEIN){
        std::cout<<"doing_eppstein"<<std::endl;
        outputfile+="_eppstein.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==COMBINED){
        std::cout<<"doing_combined"<<std::endl;
        outputfile+="_combined.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==ROOTCOMBINED){
        std::cout<<"doing_rootcombined"<<std::endl;
        outputfile+="_rootcombined.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==HEURISTIC){
        std::cout<<"doing_heuristic"<<std::endl;
        outputfile+="_heuristic.csv";
        header="file1 file2 simple eppstein hanke \n";
    }
    else if (type==BFS){
        std::cout<<"doing_BFS"<<std::endl;
        outputfile+="_bfs.csv";
        header="file1 file2 flip runtime closed open \n";
    }
    else if (type==ILPEDGE){
        std::cout<<"doing_ILP"<<std::endl;
        outputfile+="_ilp.csv";
        header="file1 file2 flip runtime\n";
    }
    else if (type==DECOMPOSITION){
        std::cout<<"doing_decompositon"<<std::endl;
        outputfile+="_decomposition.csv";
        header="file1 file2 flip runtime used_hanke\n";
    }

    std::string triangulation_path=folder_path+"/triangulation";
    std::string points_path=folder_path+"/points";
    auto file_pairs= get_file_pairs_in_folder(triangulation_path);
    auto points=load_vertices(points_path);

    std::vector<Triangle> loaded_triangulation_1= load_triangulation(file_pairs[0].first);
    std::vector<Triangle> loaded_triangulation_2= load_triangulation(file_pairs[0].second);
    HeuristicDistanceCalculator dists(points,loaded_triangulation_1,loaded_triangulation_2);

    std::ofstream file(outputfile, std::ios::app);
    std::cout<<"outputfile   "<<outputfile<<std::endl;
    file<< header;
    //int counter=1;
    for(auto p:file_pairs){
        loaded_triangulation_1= load_triangulation(p.first);
        loaded_triangulation_2= load_triangulation(p.second);

        if(type==EPPSTEIN){
            AStarFlipDistance a_star_epp(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_epp.run_epptein();
            std::string to_append=p.first+" "+p.second+" "+a_star_epp.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
            //counter++;
        }
        else if(type==COMBINED){
            AStarFlipDistance a_star_comb(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_comb.run_combined();
            std::string to_append=p.first+" "+p.second+" "+a_star_comb.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
            //counter++;
        }
        else if(type==ROOTCOMBINED){
            AStarFlipDistance a_star_comb(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_comb.run_root_combined();
            std::string to_append=p.first+" "+p.second+" "+a_star_comb.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
            //counter++;
        }

        else if (type==SIMPLE){
            AStarFlipDistance a_star_simple(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_simple.run_simple();
            std::string to_append=p.first+" "+p.second+" "+a_star_simple.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }
        else if (type==HEURISTIC){

            dists.set_new_start_triangulation(loaded_triangulation_1);
            dists.set_new_target_triangulation(loaded_triangulation_2);
            dists.run();
            std::string to_append=p.first+" "+p.second+" "+std::to_string(dists.simple_heuristic)+" "+std::to_string(dists.eppstein_heuristic)+" "+std::to_string(dists.hanke_heuristic);


            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }
        else if (type==BFS){
            BiDirectionalBfs bfs(points,loaded_triangulation_1,loaded_triangulation_2);
            bfs.run();
            dists.set_new_start_triangulation(loaded_triangulation_1);
            dists.set_new_target_triangulation(loaded_triangulation_2);
            dists.run();
            std::string to_append=p.first+" "+p.second+" "+bfs.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            std::cout<<to_append<<std::endl;
        }
        else if (type==ILPEDGE){

            AStarFlipDistance a_star_comb(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_comb.run_root_combined();

            int real_flip_distance=a_star_comb.get_flipdistance();
            if(real_flip_distance>0){
                HigherOrderDelaunayObjects obj(points);
                obj.compute_useful_order_k_edges_and_triangles(points.size());
                FlipIlpEdgeBased flip(GeometricTriangulation(loaded_triangulation_1,points),GeometricTriangulation(loaded_triangulation_2,points),obj,real_flip_distance);
                flip.solve();

                std::string to_append=p.first+" "+p.second+" "+std::to_string((int)flip.get_number_of_flips())+" "+ std::to_string((int)(flip.runtime*1000))+" "+std::to_string((int)flip.get_primal()) ;
                file << to_append<<std::endl;
                result.emplace_back(to_append);
                std::cout<<to_append<<std::endl;
            }
            else{
                HigherOrderDelaunayObjects obj(points);
                obj.compute_useful_order_k_edges_and_triangles(points.size());
                FlipIlpEdgeBased flip(GeometricTriangulation(loaded_triangulation_1,points),GeometricTriangulation(loaded_triangulation_2,points),obj);
                flip.solve();

                std::string to_append=p.first+" "+p.second+" "+std::to_string((int)flip.get_number_of_flips())+" "+ std::to_string((int)(flip.runtime*1000))+" "+std::to_string((int)flip.get_primal()) ;
                file << to_append<<std::endl;
                result.emplace_back(to_append);
                std::cout<<to_append<<std::endl;
            }



        }
        else if (type==DECOMPOSITION){
            AStarDecompositionFlipDistance a_star_decomposition(points,loaded_triangulation_1,loaded_triangulation_2);
            a_star_decomposition.run_with_combined();
            dists.set_new_start_triangulation(loaded_triangulation_1);
            dists.set_new_target_triangulation(loaded_triangulation_2);
            auto hanke=dists.run_and_return_only_hanke();
            if(hanke!=a_star_decomposition.get_flipdistance()){
                std::cout<<"better: "<<hanke-a_star_decomposition.get_flipdistance()<<std::endl;
            }

            std::string to_append=p.first+" "+p.second+" "+a_star_decomposition.get_data_as_string();
            file << to_append<<std::endl;
            result.emplace_back(to_append);
            //           std::cout<<to_append<<std::endl;
        }

    }
    return result;

}


bool customComparator(const std::tuple<std::string, int, int>& a, const std::tuple<std::string, int, int>& b) {
    if (std::get<2>(a) != std::get<2>(b)) {
        return std::get<2>(a) < std::get<2>(b); // Primary sorting by third entry
    } else {
        return std::get<1>(a) < std::get<1>(b); // Secondary sorting by second entry if third entries are equal
    }
}

std::vector<std::tuple<std::string,int,int>> get_all_folder_with_size_and_order(std::string folder_path){

    std::vector<std::tuple<std::string,int,int>> result;

    for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
        int size=-1; int order=-1;
        if (std::filesystem::is_directory(entry.path())) {
            std::string last_integer, second_last_integer;

            std::string folder_name = std::filesystem::path(entry.path()).filename().string();
            std::istringstream iss(entry.path());
            std::string token;

            while (std::getline(iss, token, '_')) {
                second_last_integer = last_integer;
                last_integer = token;
            }
            if (!last_integer.empty() && !second_last_integer.empty()) {
                size = std::stoi(last_integer);
                order = std::stoi(second_last_integer);
            } else {
                std::cout << "Unable to parse integers" << std::endl;
            }
            result.emplace_back(entry.path(),size,order);
        }

    }
    std::sort(result.begin(), result.end(), customComparator);
    return result;
}

void do_all_products_for_heuristics_random(std::string folder_path,int type){
    std::string outputfile=folder_path+"/results";
    std::string header;
    if(type==SIMPLE){
        std::cout<<"doing_simple"<<std::endl;
        outputfile+="_simple.csv";
        header="file1 file2 flip heuristic runtime closed open size order\n";


    }
    else if (type==EPPSTEIN){
        std::cout<<"doing_eppstein"<<std::endl;
        outputfile+="_eppstein.csv";
        header="file1 file2 flip heuristic runtime closed open size order\n";
    }
    else if (type==COMBINED){
        std::cout<<"doing_combined"<<std::endl;
        outputfile+="_combined.csv";
        header="file1 file2 flip heuristic runtime closed open size order\n";
    }
    else if (type==ROOTCOMBINED){
        std::cout<<"doing_rootcombined"<<std::endl;
        outputfile+="_rootcombined.csv";
        header="file1 file2 flip heuristic runtime closed open size order\n";
    }

    else if (type==HEURISTIC){
        std::cout<<"doing_heuristic"<<std::endl;
        outputfile+="_heuristic.csv";
        header="file1 file2 simple eppstein hanke size order\n";
    }
    else if (type==BFS){
        std::cout<<"doing_BFS"<<std::endl;
        outputfile+="_bfs.csv";
        header="file1 file2 flip runtime closed open size order\n";
    }
    else if (type==ILPEDGE){
        std::cout<<"doing_ILP"<<std::endl;
        outputfile+="_ilp.csv";
        header="file1 file2 flip runtime primal_bound size order\n";
    }

    else if (type==DECOMPOSITION){
        std::cout<<"doing_decompositon"<<std::endl;
        outputfile+="_decomposition.csv";
        header="file1 file2 flip runtime used_hanke  size order\n";
    }

    std::ofstream file(outputfile, std::ios::app);
    file<< header;

    auto instance_sets=get_all_folder_with_size_and_order(folder_path);

    for(auto instance :instance_sets){
        auto instance_folder=std::get<0>(instance);
        auto size=std::get<1>(instance);
        auto order=std::get<2>(instance);

        for(int i=0;i<10;i++){
            std::cout<<"Set_"<<i<<std::endl;
            std::string path_to_current_set=instance_folder+"/set_"+std::to_string(i);
            auto rows=do_one_product_of_instances(path_to_current_set,type);
            for(auto x:rows){
                file<<x<<" "<<size<<" "<<order <<std::endl;
            }
        }
    }



}



void do_all_products_for_size_random(std::string folder_path,int type){
    std::string outputfile=folder_path+"/results";
    std::string header;
    if(type==SIMPLE){
        std::cout<<"doing_simple"<<std::endl;
        outputfile+="_simple.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==EPPSTEIN){
        std::cout<<"doing_eppstein"<<std::endl;
        outputfile+="_eppstein.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==COMBINED){
        std::cout<<"doing_combined"<<std::endl;
        outputfile+="_combined.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }
    else if (type==ROOTCOMBINED){
        std::cout<<"doing_rootcombined"<<std::endl;
        outputfile+="_rootcombined.csv";
        header="file1 file2 flip heuristic runtime closed open\n";
    }

    else if (type==HEURISTIC){
        std::cout<<"doing_heuristic"<<std::endl;
        outputfile+="_heuristic.csv";
        header="file1 file2 simple eppstein hanke\n";
    }
    else if (type==BFS){
        std::cout<<"doing_BFS"<<std::endl;
        outputfile+="_bfs.csv";
        header="file1 file2 flip runtime closed open\n";
    }
    else if (type==ILPEDGE){
        std::cout<<"doing_ILP"<<std::endl;
        outputfile+="_ilp.csv";
        header="file1 file2 flip runtime primal_bound\n";
    }

    else if (type==DECOMPOSITION){
        std::cout<<"doing_decompositon"<<std::endl;
        outputfile+="_decomposition.csv";
        header="file1 file2 flip runtime used_hanke  size order\n";
    }

    std::ofstream file(outputfile, std::ios::app);
    file<< header;

    for(int i=0;i<10;i++){
        std::cout<<"Set_"<<i<<std::endl;
        std::string path_to_current_set=folder_path+"/set_"+std::to_string(i);
        auto rows=do_one_product_of_instances(path_to_current_set,type);
        for(auto x:rows){
            file<<x<<std::endl;
        }
    }

}




void do_random_experiments_fixed_one_set_with_arg(std::string type,std::string dataset,int which ,int which_secondary =1, const std::string& path_to_data="../data"){
    int type_int=-22;
    if(type=="eppstein"){
        type_int=EPPSTEIN;
    }
    if(type=="heuristic"){
        type_int=HEURISTIC;
    }
    if(type=="simple"){
        type_int=SIMPLE;
    }

    if(type=="combined"){
        type_int=COMBINED;
    }
    if(type=="rootcombined"){
        type_int=ROOTCOMBINED;
    }
    if(type=="bfs"){
        type_int=BFS;
    }
    if(type=="ilp"){
        type_int=ILPEDGE;
    }
    if(type=="allflip"){
        type_int=ALLFLIP;
    }

    if(type=="decomposition"){
        type_int=DECOMPOSITION;
    }

    std::cout<<"your path: "+path_to_data<<std::endl;

    if(dataset=="sealevel"){
        std::vector<int> valid={25,30,35,41};
        if( std::find(valid.begin(), valid.end(), which) == valid.end()){
            std::cout<<"wrong args 1 \n";
            return;
        }
        std::vector<int> valid_orders{2,5,100,which};
        if( std::find(valid_orders.begin(), valid_orders.end(), which_secondary) == valid_orders.end()){
            std::cout<<"wrong args 2\n";
            return;
        }
        std::string folder=path_to_data+"/sealevel_experiments_paper/s_"+std::to_string(which)+"/s_"+std::to_string(which)+"_"+std::to_string(which_secondary);
        do_one_product_of_instances_sea_level(folder,type_int);

    } else if(dataset=="random"){
        std::vector<int> valid={10,15,20,25,30,35,40};
        if( std::find(valid.begin(), valid.end(), which) == valid.end()){
            return;
        }
        std::string conv_or_point="";
        if(which_secondary==1){
            conv_or_point="n_";
            std::string folder=path_to_data+"/random_experiments_paper/"+conv_or_point+std::to_string(which);
            do_all_products_for_size_random(folder,type_int);
        }else if(which_secondary==0){
            conv_or_point="c_";
            std::string folder=path_to_data+"/random_experiments_paper/"+conv_or_point+std::to_string(which)+"/set_0";
            do_one_product_of_instances(folder,type_int);
        }
        else{
            std::cout<<"wrong args \n";
        }




    } else{
        std::cout<<"wrong args \n";
    }



}





int main(int argc, char *argv[]){
    if(argc!=7 ){
        std::cerr<<"wrong arguments, it should be: path_points path_t1 path_t2 which_algo output_path output_file_name"<<std::endl;
        return 17;
    }

    std::string path_vert= argv[1];
    std::string path_t1= argv[2];
    std::string path_t2= argv[3];
    std::string type= argv[4];
    std::string output_path= argv[5];
    std::string output_file_name= argv[6];
    int type_int=-1;

    std::vector<Triangle> loaded_triangulation_1= load_triangulation(path_t1);
    std::vector<Triangle> loaded_triangulation_2= load_triangulation(path_t2);
    auto points=load_vertices(path_vert);

    if(type=="eppstein"){
        type_int=EPPSTEIN;
    }
    if(type=="heuristic"){
        type_int=HEURISTIC;
    }
    if(type=="simple"){
        type_int=SIMPLE;
    }

    if(type=="combined"){
        type_int=COMBINED;
    }
    if(type=="rootcombined"){
        type_int=ROOTCOMBINED;
    }
    if(type=="bfs"){
        type_int=BFS;
    }
    if(type=="ilp"){
        type_int=ILPEDGE;
    }
    if(type=="allflip"){
        type_int=ALLFLIP;
    }

    if(type=="decomposition"){
        type_int=DECOMPOSITION;
    }

    std::cout <<type_int<<std::endl;

    if(type_int==SIMPLE){
        AStarFlipDistance a_star_epp(points,loaded_triangulation_1,loaded_triangulation_2);
        a_star_epp.run_simple();

        nlohmann::json j2 = {
                {"vertex_path", argv[1]},
                {"triangulation_1_path", argv[2]},
                {"triangulation_2_path", argv[3]},
                {"algorithm", argv[4]},
                {"flip distance", a_star_epp.get_flipdistance()},
                {"lowerbound_first_heuristic_value", a_star_epp.get_initial_heuristic()},
                {"runtime", a_star_epp.get_runtime()},
                {"closed_nodes", a_star_epp.get_nr_of_closed_nodes()},
                {"opened_nodes", a_star_epp.get_nr_of_opened_nodes()},

        };
        std::ofstream o((output_path+"/"+output_file_name));
        o << std::setw(4) << j2 << std::endl;

        std::cout<<"doing_simple"<<std::endl;

    }
    else if (type_int==EPPSTEIN){
        AStarFlipDistance a_star_epp(points,loaded_triangulation_1,loaded_triangulation_2);
        a_star_epp.run_simple();

        nlohmann::json j2 = {
                {"vertex_path", argv[1]},
                {"triangulation_1_path", argv[2]},
                {"triangulation_2_path", argv[3]},
                {"algorithm", argv[4]},
                {"flip distance", a_star_epp.get_flipdistance()},
                {"lowerbound_first_heuristic_value", a_star_epp.get_initial_heuristic()},
                {"runtime", a_star_epp.get_runtime()},
                {"closed_nodes", a_star_epp.get_nr_of_closed_nodes()},
                {"opened_nodes", a_star_epp.get_nr_of_opened_nodes()},

        };
        std::ofstream o((output_path+"/"+output_file_name));
        o << std::setw(4) << j2 << std::endl;

        std::cout<<"doing_eppstein"<<std::endl;

    }

    else if (type_int==COMBINED){
        AStarFlipDistance a_star_epp(points,loaded_triangulation_1,loaded_triangulation_2);
        a_star_epp.run_combined();

        nlohmann::json j2 = {
                {"vertex_path", argv[1]},
                {"triangulation_1_path", argv[2]},
                {"triangulation_2_path", argv[3]},
                {"algorithm", argv[4]},
                {"flip distance", a_star_epp.get_flipdistance()},
                {"lowerbound_first_heuristic_value", a_star_epp.get_initial_heuristic()},
                {"runtime", a_star_epp.get_runtime()},
                {"closed_nodes", a_star_epp.get_nr_of_closed_nodes()},
                {"opened_nodes", a_star_epp.get_nr_of_opened_nodes()},

        };
        std::ofstream o((output_path+"/"+output_file_name));
        o << std::setw(4) << j2 << std::endl;

        std::cout<<"doing_combined"<<std::endl;

    }
    else if (type_int==ROOTCOMBINED){
        AStarFlipDistance a_star_epp(points,loaded_triangulation_1,loaded_triangulation_2);
        a_star_epp.run_root_combined();

        nlohmann::json j2 = {
                {"vertex_path", argv[1]},
                {"triangulation_1_path", argv[2]},
                {"triangulation_2_path", argv[3]},
                {"algorithm", argv[4]},
                {"flip distance", a_star_epp.get_flipdistance()},
                {"lowerbound_first_heuristic_value", a_star_epp.get_initial_heuristic()},
                {"runtime", a_star_epp.get_runtime()},
                {"nodes_closed", a_star_epp.get_nr_of_closed_nodes()},
                {"nodes_opened", a_star_epp.get_nr_of_opened_nodes()},

        };
        std::ofstream o((output_path+"/"+output_file_name));
        o << std::setw(4) << j2 << std::endl;

    }
    else if (type_int==HEURISTIC){
        HeuristicDistanceCalculator dists(points,loaded_triangulation_1,loaded_triangulation_2);
        dists.set_new_start_triangulation(loaded_triangulation_1);
        dists.set_new_target_triangulation(loaded_triangulation_2);
        dists.run();
//        std::to_string(dists.simple_heuristic)+" "+std::to_string(dists.eppstein_heuristic)+" "+std::to_string(dists.hanke_heuristic);

        nlohmann::json j2 = {
                {"vertex_path", argv[1]},
                {"triangulation_1_path", argv[2]},
                {"triangulation_2_path", argv[3]},
                {"algorithm", argv[4]},
                {"simple heuristic", dists.simple_heuristic},
                {"eppstein heuristic", dists.eppstein_heuristic},
                {"hanke upperbound", dists.hanke_heuristic},
        };
        std::ofstream o((output_path+"/"+output_file_name));
        o << std::setw(4) << j2 << std::endl;

        std::cout<<"doing_heuristic"<<std::endl;

    }

    else if (type_int==BFS){
        std::cout<<"doing_heuristic"<<std::endl;
        BiDirectionalBfs bfs(points,loaded_triangulation_1,loaded_triangulation_2);
        bfs.run();

        nlohmann::json j2 = {
                {"vertex_path", argv[1]},
                {"triangulation_1_path", argv[2]},
                {"triangulation_2_path", argv[3]},
                {"algorithm", argv[4]},
                {"flip distance", bfs._flipdistance},
                {"runtime", bfs._runtime},
        };
        std::ofstream o((output_path+"/"+output_file_name));
        o << std::setw(4) << j2 << std::endl;

    }
    else if (type_int==ILPEDGE){
        std::cout<<"doing_ILP"<<std::endl;
        HeuristicDistanceCalculator dists(points,loaded_triangulation_1,loaded_triangulation_2);
        dists.run();

        HigherOrderDelaunayObjects obj(points);
        obj.compute_useful_order_k_edges_and_triangles(points.size());

        FlipIlpEdgeBased flip(GeometricTriangulation(loaded_triangulation_1, points),
                              GeometricTriangulation(loaded_triangulation_2, points), obj,dists.hanke_heuristic);
        flip.solve();

        nlohmann::json j2 = {
                {"vertex_path", argv[1]},
                {"triangulation_1_path", argv[2]},
                {"triangulation_2_path", argv[3]},
                {"algorithm", argv[4]},
                {"flip distance(dual bound)", flip.get_number_of_flips()},
                {"runtime", flip.get_runtime()},
                {"primal_bound", flip.get_primal()},
        };
        std::ofstream o((output_path+"/"+output_file_name));
        o << std::setw(4) << j2 << std::endl;

//        header="file1 file2 flip runtime primal_bound\n";

    }




    return 0;


}