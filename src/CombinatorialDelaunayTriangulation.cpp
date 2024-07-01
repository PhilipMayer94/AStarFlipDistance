

#include <AStarFlipDistance/Geometry.hpp>
#include "AStarFlipDistance/CombinatorialDelaunayTriangulation.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;

typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef Delaunay::Point Point;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_on_sphere_traits_2<K> Traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits> DToS2;
typedef Traits::Point_3 Point_3;


int getIndexFromPoint(Point_3 p, std::vector<Point_s> _vertices) {
    double x = p.x();
    double y = p.y();
    double z = p.z();
    for (int i = 0; i < _vertices.size(); i++) {
        Point_s v = _vertices[i];
        if (std::abs(x - v.x) < 0.0000001 &&
            std::abs(y - v.y) < 0.0000001 &&
            std::abs(z - v.z) < 0.0000001
                ) {
            return i;
        }
    }
    return -1;
}


bool triangle_is_not_in_vector(int i, std::vector<int> vec) {
    if (vec.empty()) {
        return true;
    }
    if (std::find(vec.begin(), vec.end(), i) != vec.end()) {
        return false;
    } else {
        return true;
    }
}

bool do_triangles_share_exactly_two_points(Triangle t, Triangle d) {
    int arr1[3] = {t.vertices[0], t.vertices[1], t.vertices[2]};
    int arr2[3] = {d.vertices[0], d.vertices[1], d.vertices[2]};
    int len = sizeof(arr1) / sizeof(arr1[0]); //get array length
    std::sort(arr1, arr1 + len);
    std::sort(arr2, arr2 + len);
    int i = 0;

    for (int n: t.vertices) {
        for (int m: d.vertices) {
            if (n == m) {
                i = i + 1;
            }
        }
    }
    if (i == 2) {
        return true;
    } else {
        return false;
    }

}


template<>
void CombinatorialDelaunayTriangulation<Point_2>::triangulate(std::vector<Point_2> const &_vertices) {
    has_boundary = true;
    std::vector<std::pair<Point, unsigned>> points;
    for (int i = 0; i < _vertices.size(); i++) {
        points.emplace_back(std::make_pair(Point(_vertices[i].x, _vertices[i].y), i));
    }
    Delaunay dt;
    dt.insert(points.begin(), points.end());

    int i = 0;

    for (int j = 0; j < _vertices.size(); j++) {
        triangles_adjacent_to_vertices.emplace_back(std::vector<int>());
    }
    for (auto fIter = dt.finite_faces_begin(); fIter != dt.finite_faces_end(); fIter++) {
        auto current_face = *fIter;
        int x = current_face.vertex(0)->info();
        int y = current_face.vertex(1)->info();
        int z = current_face.vertex(2)->info();
        triangles.emplace_back(Triangle(x, y, z));
        neighbors_of_triangles.emplace_back(std::vector<int>());
        shared_edge_of_triangle_with_neighbor.emplace_back(std::vector<Edge>());
        triangles_adjacent_to_vertices[x].push_back(i);
        triangles_adjacent_to_vertices[y].push_back(i);
        triangles_adjacent_to_vertices[z].push_back(i);
        i++;
    }
    for (i = 0; i < vertices.size(); i++) {
        for (int j = 0; j < triangles_adjacent_to_vertices[i].size(); j++) {
            for (int l = j + 1; l < triangles_adjacent_to_vertices[i].size(); l++) {
                int t1 = triangles_adjacent_to_vertices[i][j];
                int t2 = triangles_adjacent_to_vertices[i][l];


                if (do_triangles_share_exactly_two_points(triangles[t1], triangles[t2])) {

                    if (triangle_is_not_in_vector(t2, neighbors_of_triangles[t1])) {
                        neighbors_of_triangles[t1].push_back(t2);
                        Edge e = get_shared_edge(t1, t2);
                        shared_edge_of_triangle_with_neighbor[t1].push_back(e);
                    }
                    if (triangle_is_not_in_vector(t1, neighbors_of_triangles[t2])) {
                        neighbors_of_triangles[t2].push_back(t1);
                        Edge e = get_shared_edge(t1, t2);
                        shared_edge_of_triangle_with_neighbor[t2].push_back(e);
                    }
                }
            }
        }
    }
    boundary_size = vertices.size() * 2 - 2 - dt.number_of_faces();
}

void is_delaunay_assertion(std::vector<double> t, std::vector<Point_s> const &vertices) {
    t[3] = t[3] * 0.95;
    bool isOk = true;
    for (auto v: vertices) {
        if (is_point_in_circle(v, t)) {
            std::cout << " r: " << t[3];
            v.print();
            std::cout << "problempoint" << std::endl;

            isOk = false;
        }
    }
    assert(isOk);

}


template<>
void CombinatorialDelaunayTriangulation<Point_s>::triangulate(std::vector<Point_s> const &_vertices) {
    for (int j = 0; j < _vertices.size(); j++) {
        triangles_adjacent_to_vertices.emplace_back(std::vector<int>());
    }

    std::vector<Point_3> points;

    for (auto &v: _vertices) {
        points.emplace_back(v.x, v.y, v.z);
    }


    Traits traits(Point_3(0, 0, 0), 1);
    DToS2 dtos(traits);
    dtos.insert(points.begin(), points.end());
    int i = 0;
    for (auto fIter = dtos.solid_faces_begin(); fIter != dtos.solid_faces_end(); fIter++) {
        auto current_face = *fIter;
        auto x = getIndexFromPoint(current_face.vertex(0)->point(), _vertices);
        auto y = getIndexFromPoint(current_face.vertex(1)->point(), _vertices);
        auto z = getIndexFromPoint(current_face.vertex(2)->point(), _vertices);
        auto t = circle_through_three_points(vertices[x], vertices[y], vertices[z]);


        is_delaunay_assertion(t, vertices);


        triangles.emplace_back(Triangle(x, y, z));

        neighbors_of_triangles.emplace_back(std::vector<int>());
        shared_edge_of_triangle_with_neighbor.emplace_back(std::vector<Edge>());
        triangles_adjacent_to_vertices[x].push_back(i);
        triangles_adjacent_to_vertices[y].push_back(i);
        triangles_adjacent_to_vertices[z].push_back(i);
        i++;
    }

    for (int i = 0; i < vertices.size(); i++) {
        for (int j = 0; j < triangles_adjacent_to_vertices[i].size(); j++) {
            for (int l = j + 1; l < triangles_adjacent_to_vertices[i].size(); l++) {
                int t1 = triangles_adjacent_to_vertices[i][j];
                int t2 = triangles_adjacent_to_vertices[i][l];


                if (do_triangles_share_exactly_two_points(triangles[t1], triangles[t2])) {

                    if (triangle_is_not_in_vector(t2, neighbors_of_triangles[t1])) {
                        neighbors_of_triangles[t1].push_back(t2);
                        Edge e = get_shared_edge(t1, t2);
                        shared_edge_of_triangle_with_neighbor[t1].push_back(e);
                    }
                    if (triangle_is_not_in_vector(t1, neighbors_of_triangles[t2])) {
                        neighbors_of_triangles[t2].push_back(t1);
                        Edge e = get_shared_edge(t1, t2);
                        shared_edge_of_triangle_with_neighbor[t2].push_back(e);
                    }
                }
            }
        }
    }
    if (dtos.number_of_ghost_faces() == 0) {
        std::cout << "triangulation has no boundary" << std::endl;
        has_boundary = false;
        boundary_size = 0;
    } else {
        has_boundary = true;
        boundary_size = vertices.size() * 2 - 2 - dtos.number_of_solid_faces();
    }

}
