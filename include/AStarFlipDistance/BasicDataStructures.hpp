#ifndef A_STAR_FOR_FLIPDISTANCE__BASICDATASTRUCTURES_HPP
#define A_STAR_FOR_FLIPDISTANCE__BASICDATASTRUCTURES_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>




struct Point_2 {
    double x, y;

    Point_2() = default;

    Point_2(double p1, double p2) {
        x = p1;
        y = p2;
    }

    bool operator==(const Point_2 &p) const {
        if (this->x == p.x && this->y == p.y) {
            return true;
        }
        return false;
    }

    bool roughly_equals(const Point_2 &p) const {
        double eps=0.000001;
        if (std::abs(this->x - p.x)<eps && std::abs(this->y-p.y)<eps) {
            return true;
        }
        return false;
    }

    Point_2 operator+(const Point_2 &obj) const {
        Point_2 tmp{};
        tmp.x = x + obj.x;
        tmp.y = y + obj.y;
        return tmp;

    }

    Point_2 operator-(const Point_2 &obj) const {
        Point_2 tmp{};
        tmp.x = x - obj.x;
        tmp.y = y - obj.y;
        return tmp;

    }

    void print() const {
        std::cout << "(" << x << ", " << y << ")";
    }


};

struct Point_s {
public:
    double x{}, y{}, z{}, lat{MAXFLOAT}, lon{MAXFLOAT}, r{};

    Point_s() = default;

    Point_s(double x1, double x2, double x3, bool normalize = true) {
        x = x1;
        y = x2;
        z = x3;


        if (normalize) {
            r = sqrt(x * x + y * y + z * z);
            x = x / r;
            y = y / r;
            z = z / r;
            r = 1;
        }
        //lat = asin(z / r);
        //lon = atan2(y, x);
    }

    Point_s(double longitude, double latitude) {

        r = 1;
        lat = latitude*M_PI/180;
        lon = longitude*M_PI/180;

        //lat in -90 --- 90
        x = cos(lat) * cos(lon);
        y = cos(lat) * sin(lon);
        z = sin(lat);

        //lat in 0 --- 180
        //x = sin(lat) * cos(lon);
        //y = sin(lat) * sin(lon);
        //z = cos(lat);

        r= sqrt(x*x+y*y+z*z);
        x = x / r;
        y = y / r;
        z = z / r;
        r = 1;

    }

    void setXYZ(double _x, double _y, double _z){
        x=_x;
        y=_y;
        z=_z;
    }

    double get_norm(){
        return sqrt(x*x+y*y+z*z);
    }

    Point_s get_a_normalized_point(){
        double x_res=x/get_norm();
        double y_res=y/get_norm();
        double z_res=z/get_norm();
        return{x_res,y_res,z_res};
    }

    void set_lat_lon(){
        if(lat==MAXFLOAT){
            r = sqrt(x * x + y * y + z * z);
            lat = asin(z / r);
            lon = atan2(y, x);
        }

    }

    Point_s operator+(const Point_s &obj) const {
        Point_s tmp;
        tmp.x = x + obj.x;
        tmp.y = y + obj.y;
        tmp.z = z + obj.z;
        return tmp;

    }

    Point_s operator-(const Point_s &obj) const {
        Point_s tmp;
        tmp.x = x - obj.x;
        tmp.y = y - obj.y;
        tmp.z = z - obj.z;
        return tmp;

    }

    void print() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")";
    }

    bool operator==(const Point_s &p) const {
        if (this->x == p.x && this->y == p.y && this->z == p.z) {
            return true;
        }
        return false;
    }
};

template <typename T>
struct PointCloud
{
    struct Point
    {
        T x, y, z;
    };

    using coord_t = T;  //!< The type of each coordinate

    std::vector<Point> pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline T kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0)
            return pts[idx].x;
        else if (dim == 1)
            return pts[idx].y;
        else
            return pts[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }

    void addPoint(Point_s q){
        Point p;
        p.x=q.x;
        p.y=q.y;
        p.z=q.z;
        pts.emplace_back(p);

    }
    void addPoint(Point_2 q){
        Point p;
        p.x=q.x;
        p.y=q.y;
        pts.emplace_back(p);

    }
};

struct Triangle {
    std::array<int, 3> vertices{};


    Triangle(int t1, int t2, int t3) {
        int tmp1[3];
        tmp1[0] = t1;
        tmp1[1] = t2;
        tmp1[2] = t3;
        int len = sizeof(tmp1) / sizeof(tmp1[0]);
        std::sort(tmp1, tmp1 + len);
        vertices[0] = tmp1[0];
        vertices[1] = tmp1[1];
        vertices[2] = tmp1[2];

        //vertices[0] = t1;
        //vertices[1] = t2;
        //vertices[2] = t3;


    }

    bool operator==(Triangle t) const {
        int tmp1[3];
        int tmp2[3];
        tmp1[0] = t.vertices[0];
        tmp1[1] = t.vertices[1];
        tmp1[2] = t.vertices[2];

        tmp2[0] = vertices[0];
        tmp2[1] = vertices[1];
        tmp2[2] = vertices[2];
        int len = sizeof(tmp1) / sizeof(tmp1[0]); //get array length
        std::sort(tmp1, tmp1 + len);
        std::sort(tmp2, tmp2 + len);

        if (tmp1[0] == tmp2[0] && tmp1[1] == tmp2[1] && tmp1[2] == tmp2[2]) {
            return true;
        }
        return false;

    }

    int third_point(int v_1, int v_2) {
        for (int i = 0; i < 3; i++) {
            if (vertices[i] != v_1 && vertices[i] != v_2) {
                return vertices[i];
            }
        }
        return -1;

    }


    bool operator<(const Triangle &t) const {
        int tmp1[3];
        int tmp2[3];
        tmp2[0] = t.vertices[0];
        tmp2[1] = t.vertices[1];
        tmp2[2] = t.vertices[2];

        tmp1[0] = vertices[0];
        tmp1[1] = vertices[1];
        tmp1[2] = vertices[2];
        int len = sizeof(tmp1) / sizeof(tmp1[0]); //get array length
        std::sort(tmp1, tmp1 + len);
        std::sort(tmp2, tmp2 + len);

        if (tmp1[0] < tmp2[0]) { return true; }
        if (tmp1[0] > tmp2[0]) { return false; }
        if (tmp1[1] < tmp2[1]) { return true; }
        if (tmp1[1] > tmp2[1]) { return false; }
        if (tmp1[2] < tmp2[2]) { return true; }
        if (tmp1[2] > tmp2[2]) { return false; }

        return false;
    }


    void print() const {
        int tmp1[3];
        tmp1[0] = vertices[0];
        tmp1[1] = vertices[1];
        tmp1[2] = vertices[2];
        int len = sizeof(tmp1) / sizeof(tmp1[0]);
        std::sort(tmp1, tmp1 + len);


        std::cout << "(" << tmp1[0] << ", " << tmp1[1] << ", " << tmp1[2] << ")";
    }

};

struct triangle_hash {
    std::size_t operator()(const Triangle &t) const {
        int tmp1[3];
        tmp1[0] = t.vertices[0];
        tmp1[1] = t.vertices[1];
        tmp1[2] = t.vertices[2];

        int hash = 23;
        hash = hash * 31 + tmp1[0] + 1;
        hash = hash * 31 + tmp1[1] + 1;
        hash = hash * 31 + tmp1[2] + 1;

        return hash;

    }
};

struct triangle_pair_hash {
    std::size_t operator()(const std::pair<Triangle, int> &p) const {
        Triangle t = p.first;
        int tmp1[3];
        tmp1[0] = t.vertices[0];
        tmp1[1] = t.vertices[1];
        tmp1[2] = t.vertices[2];
        int len = sizeof(tmp1) / sizeof(tmp1[0]); //get array length

        std::sort(tmp1, tmp1 + len);

        int hash = 23;
        hash = hash * 31 + tmp1[0] + 1;
        hash = hash * 31 + tmp1[1] + 1;
        hash = hash * 31 + tmp1[2] + 1;

        return hash;

    }
};


struct Edge {
    int src;
    int dst;

    Edge(int u, int v) {
        src = std::min(u, v);
        dst = std::max(u, v);
    }

    bool operator<(const Edge &e) const {
        if (src < e.src) {
            return true;
        }
        if (src == e.src && dst < e.dst) {
            return true;
        }
        return false;
    }

    bool operator==(const Edge &p) const {
        if (this->src == p.src && this->dst == p.dst) {
            return true;
        }
        return false;
    }

    void print() const {
        std::cout << "(" << src << ", " << dst << ")";
    }
};

struct edge_hash {
    std::size_t operator()(const Edge &p) const {
        int tmp1[2];
        tmp1[0] = p.src;
        tmp1[1] = p.dst;
        int len = sizeof(tmp1) / sizeof(tmp1[0]); //get array length

        std::sort(tmp1, tmp1 + len);

        int hash = 23;
        hash = hash * 31 + tmp1[0] + 1;
        hash = hash * 31 + tmp1[1] + 1;
        return hash;

    }
};


#endif //HIGHERORDERDELAUNAY_BASICDATASTRUCTURES_HPP
