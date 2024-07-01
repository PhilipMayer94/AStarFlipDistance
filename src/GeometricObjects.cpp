#include <queue>
#include "AStarFlipDistance/GeometricObjects.hpp"
#include "AStarFlipDistance/HigherOrderDelaunayObjects.hpp"

class UpperTriangularMatrix {
private:
    int size;
    std::vector<int> matrix; // 1D array to store upper triangular elements

    // Convert 2D index to 1D index for upper triangular matrix
    int getIndex(int row, int col) {
        // Ensure row <= col as it's an upper triangular matrix
        if (row > col) {
            std::swap(row, col);
        }
        return row * size - ((row - 1) * row) / 2 + col - row;
    }

public:
    UpperTriangularMatrix(int n) : size(n) {
        // Initialize the 1D array with zeros or default values
        int totalElements = (n * (n + 1)) / 2;
        matrix.assign(totalElements, 0);
    }

    // Set value at row, col (0-based indexing)
    void setValue(int row, int col, int value) {
        if (row >= size || col >= size) {
            // Handle out of bounds error
            return;
        }
        int index = getIndex(row, col);
        matrix[index] = value;
    }

    // Get value at row, col (0-based indexing)
    int getValue(int row, int col) {
        if (row >= size || col >= size) {
            // Handle out of bounds error
            return -1; // Or any default value
        }
        int index = getIndex(row, col);
        return matrix[index];
    }
};

GeometricObjects::GeometricObjects(std::vector<Point_2> vertices, bool compute_distance) {
    _vertices = std::move(vertices);
    n = (int) _vertices.size();
    std::vector<std::vector<int>> edge_matrix(n, std::vector<int>(n, -1));
    _edges_matrix = edge_matrix;
    std::vector<std::vector<std::vector<int>>> triangle_matrix(n, std::vector<std::vector<int>>(n, std::vector<int>(n,
                                                                                                                    -1)));
    _triangle_matrix = triangle_matrix;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int ind = (int) _edge_list.size();
            _edge_list.emplace_back(i, j);
            _edges_matrix[i][j] = ind;
        }
    }
    _triangle_list.clear();
    _triangle_is_allowed_list.clear();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                int ind = (int) _triangle_list.size();
                Triangle tmp(i, j, k);
                bool is_empty = is_triangle_empty_as_indices_in_set(i, j, k, _vertices);
                _triangle_list.emplace_back(tmp);
                _triangle_is_allowed_list.emplace_back(is_empty);
                _triangle_matrix[i][j][k] = ind;
                _triangle_matrix[i][k][j] = ind;
                _triangle_matrix[j][k][i] = ind;
                _triangle_matrix[j][i][k] = ind;
                _triangle_matrix[k][j][i] = ind;
                _triangle_matrix[k][i][j] = ind;

            }
        }
    }

    m = (int) _edge_list.size();
    std::vector<std::vector<bool>> quad_matrix(m, std::vector<bool>(m, 0));
    _allowed_flips = quad_matrix;

    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            bool allowed = false;


            Edge e_1 = _edge_list[i];
            Edge e_2 = _edge_list[j];


            if (e_1.src != e_2.src && e_1.src != e_2.dst && e_1.dst != e_2.src && e_1.dst != e_2.dst) {

                bool convex = is_quadrailateral_convex(_vertices[e_1.src], _vertices[e_1.dst], _vertices[e_2.src],
                                                       _vertices[e_2.dst]);
                Triangle t1(e_1.src, e_1.dst, e_2.src);
                Triangle t2(e_1.src, e_1.dst, e_2.dst);




                bool triangles_empty =
                           _triangle_is_allowed_list[_triangle_matrix[t1.vertices[0]][t1.vertices[1]][t1.vertices[2]]]
                        && _triangle_is_allowed_list[_triangle_matrix[t2.vertices[0]][t2.vertices[1]][t2.vertices[2]]];
//
//                bool triangles_empty_2 =
//                        is_triangle_empty_as_indices_in_set(t1.vertices[0], t1.vertices[1], t1.vertices[2], _vertices)
//                        &&
//                        is_triangle_empty_as_indices_in_set(t2.vertices[0], t2.vertices[1], t2.vertices[2], _vertices);

                allowed = convex && triangles_empty;
            } else {
                allowed = false;
            }

            if (!do_arcs_intersect(_vertices[e_1.src], _vertices[e_1.dst], _vertices[e_2.src], _vertices[e_2.dst])) {
                allowed = false;
            }

            _allowed_flips[i][j] = allowed;

        }
    }
    if (compute_distance) {
        compute_distance_matrix();
    }
    //restrict_flips_to_order(10);



}

void GeometricObjects::restrict_flips_to_order(int order) {
    int sizee = _vertices.size();
    HigherOrderDelaunayObjects hod(_vertices);
    hod.compute_useful_order_k_edges_and_triangles(order);
    std::vector<std::vector<std::vector<int>>> threeDVector(sizee, std::vector<std::vector<int>>(sizee,
                                                                                                 std::vector<int>(sizee,
                                                                                                                  0)));
    for (auto t: hod.get_useful_triangles()) {
        threeDVector[t.vertices[0]][t.vertices[1]][t.vertices[2]] = 1;
    }
    for (int i = 0; i < _allowed_flips.size(); i++) {
        for (int j = 0; j < _allowed_flips.size(); j++) {
            Edge e1 = get_edge(i);
            Edge e2 = get_edge(j);
            Triangle t1(e1.src, e1.dst, e2.src);
            Triangle t2(e1.src, e1.dst, e2.dst);
            Triangle t3(e1.src, e2.src, e2.dst);
            Triangle t4(e1.dst, e2.src, e2.dst);
            std::vector<Triangle> triangles = {t1, t2, t3, t4};
            int counter = 0;
            for (auto t: triangles) {
                counter = counter + threeDVector[t.vertices[0]][t.vertices[1]][t.vertices[2]];
            }
            if (counter != 4) {
                _allowed_flips[i][j] = false;
            }


        }
    }

}

void GeometricObjects::print_stuff() {

    for (const auto &row: _edges_matrix) {
        for (const auto &element: row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
    for (auto e: _edge_list) {
        e.print();
        std::cout << " ";
    }
    std::cout << std::endl;

    for (auto e: _triangle_list) {
        e.print();
        std::cout << " ";
    }
    std::cout << std::endl;
    for (auto e: _triangle_is_allowed_list) {
        std::cout << e << " ";
    }


    std::cout << std::endl;
    std::cout << _triangle_list.size() << " size";
    std::cout << std::endl;
    std::cout << _vertices.size() << " size";
    std::cout << std::endl;

    for (const auto &row: _allowed_flips) {
        for (const auto &element: row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }


}

int GeometricObjects::get_number_of_triangles() {
    return _triangle_list.size();
}

Triangle GeometricObjects::get_triangle(int i) {
    return _triangle_list[i];
}

int GeometricObjects::get_triangle(Triangle t) {
    return _triangle_matrix[t.vertices[0]][t.vertices[1]][t.vertices[2]];
}

bool GeometricObjects::is_flip_valid(int e_0, int e_1, int f_0, int f_1) {
    int edge_1 = _edges_matrix[std::min(e_0, e_1)][std::max(e_0, e_1)];
    int edge_2 = _edges_matrix[std::min(f_0, f_1)][std::max(f_0, f_1)];
    if (std::min(edge_1, edge_2) < 0) {
        return false;
    }
    return (_allowed_flips[std::min(edge_1, edge_2)][std::max(edge_1, edge_2)]);
}

bool GeometricObjects::is_flip_valid(int edge_1, int edge_2) {
    int smaller = edge_1;
    int larger = edge_2;
    if (edge_1 > edge_2) {
        smaller = edge_2;
        larger = edge_1;
    }

    return (_allowed_flips[smaller][larger]);
}

int GeometricObjects::get_number_of_vertices() {
    return _vertices.size();
}

int GeometricObjects::get_edge(int e_0, int e_1) {
    return _edges_matrix[std::min(e_0, e_1)][std::max(e_0, e_1)];;
}

Edge GeometricObjects::get_edge(int e_index) {
    return _edge_list[e_index];
}

int GeometricObjects::get_number_of_edges() {
    return _edge_list.size();
}

int GeometricObjects::get_triangle(int i, int j, int k) {
    return _triangle_matrix[i][j][k];
}

void GeometricObjects::compute_distance_matrix() {
    std::vector<int> visited(_edge_list.size(), -1);
    _edge_distance_matrix = std::vector<std::vector<triangle_int>>(_edge_list.size(),
                                                              std::vector<triangle_int>(_edge_list.size(), {10000}));
    std::queue<std::pair<int, int>> queue;
    for (int start = 0; start < _edge_list.size(); start++) {
        queue = std::queue<std::pair<int, int>>();
        queue.emplace(start, 0);
        visited[start] = start;
        while (!queue.empty()) {
            auto current = queue.front();
            queue.pop();
            _edge_distance_matrix[start][current.first] = static_cast<triangle_int>(current.second);

            for (int next = 0; next < _edge_list.size(); next++) {
                if (is_flip_valid(current.first, next) && visited[next] != start) {
                    queue.emplace(next, current.second + 1);
                    visited[next] = start;
                }
            }
        }
    }
    std::vector<triangle_int> edge_distance_flattened;
    for (int i = 0; i < _edge_distance_matrix.size(); i++) {
        for (int j = 0; j < _edge_distance_matrix.size(); j++) {
            edge_distance_flattened.emplace_back(_edge_distance_matrix[i][j]);
        }
    }
    _edge_distance_flattened = edge_distance_flattened;
    _edge_matrix_size = _edge_distance_matrix.size();

}

int GeometricObjects::get_distance(int edge_ind_1, int edge_ind_2) {

    //return static_cast<int>( _edge_distance_matrix[edge_ind_1][edge_ind_2]);

    return static_cast<int>(_edge_distance_flattened[edge_ind_1 * _edge_matrix_size + edge_ind_2]);
}

std::vector<Point_2> GeometricObjects::get_vertices() {
    return _vertices;
}
