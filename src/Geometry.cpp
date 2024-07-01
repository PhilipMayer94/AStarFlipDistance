
#include <AStarFlipDistance/Geometry.hpp>





bool is_quadrailateral_convex(Point_2 a, Point_2 b, Point_2 c, Point_2 d) {
    double signABC = (a.y - b.y) * (c.x - b.x) - (a.x - b.x) * (c.y - b.y);
    double signABD = (a.y - b.y) * (d.x - b.x) - (a.x - b.x) * (d.y - b.y);
    double signBCD = (b.y - c.y) * (d.x - c.x) - (b.x - c.x) * (d.y - c.y);
    double signCAD = (c.y - a.y) * (d.x - a.x) - (c.x - a.x) * (d.y - a.y);

    int counterplus=0;
    int counterminus=0;

    if(signABC>=0){counterplus++;}else{counterminus++;}
    if(signABD>=0){counterplus++;}else{counterminus++;}
    if(signBCD>=0){counterplus++;}else{counterminus++;}
    if(signCAD>=0){counterplus++;}else{counterminus++;}

   return counterplus==3 || counterminus==3;
}



double determinante(Point_s v1,Point_s v2,Point_s v3) {
    return(scalarProduct(get_scaled(100,v1),crossProduct(get_scaled(100,v2),get_scaled(100,v3))));
}

std::vector<double> get_interior_angles_of_triangle(Point_s A, Point_s B,Point_s C){
    double a = distance(C,B);
    double b = distance(C,A);
    double c = distance(A,B);


    // From Cosine law
    auto alpha = (double) acos(   ( cos(a)- (cos(c)*cos(b) ) )
                                       /  (sin(c)*sin(b))    );

    auto beta = (double) acos(    ( cos(b)- (cos(a)*cos(c) ) )
                                       /  (sin(a)*sin(c))    );

    auto gamma = (double) acos(   ( cos(c)- (cos(a)*cos(b) ) )
                                       /  (sin(a)*sin(b))    );

    // Converting to degree
    //alpha = (float) (alpha * 180 / Math.PI);
    //beta = (float) (beta * 180 / Math.PI);
    //gamma = (float) (gamma * 180 / Math.PI);

    return {alpha,beta,gamma};
}

double get_area_of_triangle(Point_s A, Point_s B, Point_s C){
    std::vector<double> angles= get_interior_angles_of_triangle(A,B,C);
    return angles[0]+angles[1]+angles[2]-M_PI;
}


bool straddles(Point_s a1,Point_s a2, Point_s b1,Point_s b2){
    return(determinante(a1,b1,b2)*determinante(a2,b1,b2)<0);
}



double distance(Point_2 a, Point_2 b) {
    return (sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)));
}

double distance(Point_s _p, Point_s _q) {
    double x = 0;

    if(_p==_q){
        return 0;
    }


    Point_s p=_p.get_a_normalized_point();
    p.r=1;
    Point_s q=_q.get_a_normalized_point();
    q.r=1;
    if (p.r != 0 || q.r != 0) {
        x = p.x / p.r * q.x / q.r + p.y / p.r * q.y / q.r + p.z / p.r * q.z / q.r;
    } else {
        x = p.x * q.x + p.y * q.y + p.z * q.z;
    }

    return acos(x);
}

double distance(Point_2 a, std::vector<double> &circle) {
    Point_2 tmp(circle[0], circle[1]);
    return distance(a, tmp);


}

double distance(Point_s a, std::vector<double> &circle) {
    Point_s tmp(circle[0], circle[1], circle[2]);
    return distance(a, tmp);
}


double area_of_triangle(Point_2 p1, Point_2 p2, Point_2 p3){

    double a = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    double b = sqrt(pow(p3.x - p2.x, 2) + pow(p3.y - p2.y, 2));
    double c = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2));
    double s = (a + b + c) / 2;
    return sqrt(s * (s - a) * (s - b) * (s - c));
}

bool is_point_in_circle(Point_2 a, std::vector<double> &circle) {
    if (distance(a, circle) < circle[2]) {
        return true;
    }
    return false;
}

bool is_triangle_empty_as_indices_in_set(int t1,int t2, int t3, std::vector<Point_2> &points){
    double EPS=0.001;
    for(int i=0;i<points.size();i++){
        auto p=points[i];
        if(i==t1|| i ==t2|| i==t3 ){
            continue;
        }
        if(is_point_in_triangle(p,points[t1],points[t2],points[t3])){
            return false;
        }
    }
    return true;
}




bool is_triangle_empty(Point_2 t1,Point_2 t2, Point_2 t3, std::vector<Point_2> &points){
    double EPS=0.001;
    for(auto p: points){
        if(distance(t1,p)<EPS||distance(t2,p)<EPS||distance(t3,p)<EPS){
            continue;
        }
        if(is_point_in_triangle(p,t1,t2,t3)){
            return false;
        }
    }
    return true;
}

bool is_triangle_empty(Point_s t1,Point_s t2, Point_s t3, std::vector<Point_s> &points){
    double EPS=0.001;
    int i=0;
    for(auto p: points){
        if(t1==p||t2==p||t3==p){
            i++;
            continue;
        }
        else{
            if(is_point_in_triangle(p,t1,t2,t3)){
                return false;
            }
        }
        i++;
    }
    return true;
}



bool is_point_in_circle(Point_s a, std::vector<double> &circle) {
    if (distance(a, circle) < circle[3]) {
        return true;
    }
    return false;
}


bool is_point_in_triangle2(Point_s p, Point_s t_1, Point_s t_2, Point_s t_3){
    Point_s n0= crossProduct(t_1,t_2);
    Point_s n1=crossProduct(t_2,t_3);
    Point_s n2=crossProduct(t_3,t_1);
    if(scalarProduct(n0,p)>0){
        if(scalarProduct(n1,p)>0&&scalarProduct(n2,p)>0){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        if(scalarProduct(n1,p)<=0&&scalarProduct(n2,p)<=0){
            return true;
        }
        else{
            return false;
        }

    }
    return false;
}

bool is_point_in_triangle(Point_s p, Point_s t_1, Point_s t_2, Point_s t_3) {
    double EPSILON = 0.000001;

    Point_s ray = p;

    Point_s v0 = t_1;
    Point_s v1 = t_2;
    Point_s v2 = t_3;

    Point_s e1{}, e2{};
    Point_s h{}, s{}, q{};

    double a, f, u, v;

    e1 = v1-v0;
    e2 = v2- v0;
    h = crossProduct(ray, e2);
    a = scalarProduct(e1, h);

    if (a > -EPSILON && a < EPSILON) {
        return false;
    }

    f = 1.0 / a;
    s = Point_s(0, 0, 0,false)- v0;
    u = f * (scalarProduct(s, h));
    if (u < 0.0 || u > 1.0) {
        return false;

    }
    q = crossProduct(s, e1);

    v = f * scalarProduct(ray, q);
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }
    double t = f * scalarProduct(e2, q);
    if (t > EPSILON) {
        return true;
    } else {
        return false;
    }
}

bool is_left_oriented(Point_2 v, Point_2 e_1, Point_2 e_2) {
    double pos = (e_2.x - e_1.x) * (v.y - e_1.y) - (e_2.y - e_1.y) * (v.x - e_1.x);
    if (pos < 0) {
        return true;
    }
    if (pos > 0) {
        return false;
    }
    return false;
}

bool is_left_oriented(Point_s p, Point_s e_1, Point_s e_2) {
    Point_s normal = crossProduct(e_1, e_2);
    double scalar = scalarProduct(normal, p);
    if (scalar > 0) {
        return true;
    } else {
        return false;
    }
}

bool is_point_in_triangle(Point_2 p, Point_2 t_1, Point_2 t_2, Point_2 t_3) {
    double denominator = ((t_2.y - t_3.y) * (t_1.x - t_3.x) + (t_3.x - t_2.x) * (t_1.y - t_3.y));

    double a = ((t_2.y - t_3.y) * (p.x - t_3.x) + (t_3.x - t_2.x) * (p.y - t_3.y)) / denominator;

    double b = ((t_3.y - t_1.y) * (p.x - t_3.x) + (t_1.x - t_3.x) * (p.y - t_3.y)) / denominator;

    double c = 1 - a - b;

    return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;
}


bool do_arcs_intersect_including_endpoints(Point_2 p_1, Point_2 p_2, Point_2 q_1, Point_2 q_2) {
    Point_2 p1 = p_1;
    Point_2 q1 = p_2;
    Point_2 p2 = q_1;
    Point_2 q2 = q_2;

    if ((p_1.roughly_equals( q_1) && p_2.roughly_equals( q_2)) || (p_2.roughly_equals(  q_1 )&& p_1.roughly_equals(  q_2)) ) {
        return true;
    }
    if (p_1 == q_1 || p_1 == q_2 || p_2 == q_1 || p_2 == q_2 || q_1 == q_2 || p_1 == p_2) {
        return false;
    }


    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false;
}

bool do_arcs_intersect(Point_2 p_1, Point_2 p_2, Point_2 q_1, Point_2 q_2) {
    Point_2 p1 = p_1;
    Point_2 q1 = p_2;
    Point_2 p2 = q_1;
    Point_2 q2 = q_2;

    if (p_1 == q_1 || p_1 == q_2 || p_2 == q_1 || p_2 == q_2 || q_1 == q_2 || p_1 == p_2) {
        return false;
    }


    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false;
}

bool do_arcs_intersect(Point_s p1, Point_s p2, Point_s q1, Point_s q2) {

    if (p1 == q1 || p1 == q2 || p2 == q1 || p2 == q2 || q1 == q2 || p1 == p2) {
        return false;
    }


    Point_s l = crossProduct(crossProduct(p1, p2), crossProduct(q1, q2));
    Point_s l1(l.x, l.y, l.z);
    Point_s tmp(0, 0, 0, false);
    Point_s cand = getLineIntersection(p1, p2, tmp, l1);

    double thetaStart = sqrt(normSquared((cand- p1)));
    double thetaEnd = sqrt(normSquared((cand- p2)));
    double thetaAll = sqrt(normSquared((p1- p2)));

    Point_s tmp2(0, 0, 0, false);
    Point_s cand2 = getLineIntersection(q1, q2, tmp2, l1);
    double thetaStart2 = sqrt(normSquared((cand2- q1)));
    double thetaEnd2 = sqrt(normSquared((cand2- q2)));
    double thetaAll2 = sqrt(normSquared((q1- q2)));
    bool intersects = true;


    if (cand2.x / cand.x < 0) { return false; }
    if (std::abs((thetaAll - thetaStart - thetaEnd)) >= 0.0000001 ||
        std::abs((thetaAll2 - thetaStart2 - thetaEnd2)) >= 0.0000001) {
        intersects = false;
    }
    return intersects;
}

bool do_arcs_intersect_straddle(Point_s p1, Point_s q1, Point_s p2, Point_s q2){
    return straddles(p1,q1,p2,q2);
}

Point_s getLineIntersection(Point_s p_1, Point_s p_2, Point_s q1, Point_s q2) {
    Point_s p1 = p_1;
    Point_s d1 = (p_2- p_1);
    Point_s p2 = q1;
    Point_s d2 = (q2- q1);
    Point_s n = crossProduct(d1, d2);
    Point_s n2 = crossProduct(d2, n);


    double fac = scalarProduct((p2- p1), n2);
    double bottom = scalarProduct(d1, n2);
    fac = fac / bottom;
    Point_s result = p1+get_scaled(fac, d1);
    return result;
}

std::vector<double> circle_through_three_points(Point_2 b, Point_2 c, Point_2 d) {
    double temp = c.x * c.x + c.y * c.y;
    double bc = (b.x * b.x + b.y * b.y - temp) / 2;
    double cd = (temp - d.x * d.x - d.y * d.y) / 2;
    double det = (b.x - c.x) * (c.y - d.y) - (c.x - d.x) * (b.y - c.y);

    if (std::abs(det) < 1.0e-10) {
        std::cout << "error cocircular";
        std::vector<double> result = {0, 0, -1};
        return result;
    }
    double cx = (bc * (c.y - d.y) - cd * (b.y - c.y)) / det;
    double cy = ((b.x - c.x) * cd - (c.x - d.x) * bc) / det;
    double radius = sqrt((cx - b.x) * (cx - b.x) + (cy - b.y) * (cy - b.y));

    std::vector<double> result = {cx, cy, radius};
    return result;
}

std::vector<double> circle_through_three_points(Point_s a, Point_s b, Point_s c) {
    Point_s center = circumcenter_spherical_triangle(a, b, c);
    double r = distance(center, b);

    std::vector<double> result = {center.x, center.y, center.z, r};
    return result;
}

Point_s circumcenter_spherical_triangle(Point_s a, Point_s b, Point_s c) {
    Point_s d3 = circumcenter_3D_triangle(a, b, c);
    Point_s result(d3.x, d3.y, d3.z);
    return result;
}

Point_s circumcenter_3D_triangle(Point_s a, Point_s b, Point_s c) {

    Point_s cMa = (c- a);
    Point_s bMa = (b- a);
    Point_s test=crossProduct(bMa, cMa);
    Point_s tmp1 = get_scaled(normSquared(cMa), crossProduct(crossProduct(bMa, cMa), bMa));
    Point_s tmp2 = get_scaled(normSquared(bMa), crossProduct(crossProduct(cMa, bMa), cMa));
    Point_s topCombined = tmp1+ tmp2;


    double bottom = 1 / (2 * normSquared(crossProduct(bMa, cMa)));

    Point_s combined = get_scaled(bottom, topCombined);


    Point_s result=combined+a;
    return result;

}

void test() {
    Triangle t(1, 1, 1);
}


int orientation(Point_2 p, Point_2 q, Point_2 r) {
    double val = (q.y - p.y) * (r.x - q.x) -
                 (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; // colinear
    return (val > 0) ? 1 : 2;
}

bool onSegment(Point_2 p, Point_2 q, Point_2 r) {
    if (q.x <= fmax(p.x, r.x) && q.x >= fmax(p.x, r.x) &&
        q.y <= fmax(p.y, r.y) && q.y >= fmax(p.y, r.y))
        return true;

    return false;
}


Point_s crossProduct(Point_s p1, Point_s p2) {
    double x = p1.y * p2.z - p1.z * p2.y;
    double y = p1.z * p2.x - p1.x * p2.z;
    double z = p1.x * p2.y - p1.y * p2.x;
    Point_s result(x, y, z,false);

    return result;
}

double scalarProduct(Point_s p1, Point_s p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}


double normSquared(Point_s a) {
    return (a.x * a.x + a.y * a.y + a.z * a.z);
}

Point_s get_scaled(double scalar, Point_s point) {
    double x = scalar * point.x;
    double y = scalar * point.y;
    double z = scalar * point.z;
    Point_s result(x, y, z, false);
    return result;
}


bool edge_in_polygon(int e_1, int e_2,  const std::vector<Point_2>& polygon) {

    for (int i = 0; i < polygon.size(); i++) {

        if ((e_1 == i && e_2 == (i + 1) % polygon.size()) ||
            (e_2 == i && e_1 == (i + 1) % polygon.size())) {
            return true;
        }

        if (do_arcs_intersect(polygon[e_1], polygon[e_2], polygon[i], polygon[(i + 1) % polygon.size()])) {
            return false;
        }

    }
    if (!has_edge_angle_inside_polygon(e_1, e_2, polygon)) {
        return true;
    }
    return false;
}

    bool edge_in_polygon(int e_1, int e_2,  const std::vector<Point_s>& polygon){


    for(int i=0;i<polygon.size();i++){

        if((e_1==i&&e_2==(i+1)%polygon.size()) ||
        (e_2==i&&e_1==(i+1)%polygon.size())){
            return true;
        }

        if(do_arcs_intersect(polygon[e_1],polygon[e_2],polygon[i],polygon[(i+1)%polygon.size()])){
            return false;
        }

    }

    if(has_edge_angle_inside_polygon(e_1,e_2,polygon)){
        return true;
    }
    return false;

}




double angle_between_anchored_vectors(Point_2 anchor,Point_2 main, Point_2 target){
    if(main==target){
        return 2*M_PI;
    }
    double targetX=target.x-anchor.x;
    double targetY=target.y-anchor.y;
    double preEdgeX=main.x-anchor.x;
    double preEdgeY=main.y-anchor.y;

    double targetA=atan2(targetY, targetX);

    double preEdgeA=atan2(preEdgeY, preEdgeX);

    targetA=targetA-preEdgeA;

    targetA=normalizeAngle(targetA);
    return targetA;
}

bool has_edge_angle_inside_polygon(int e_1, int e_2,  const std::vector<Point_2>& polygon){
    int root=e_1;
    int target=e_2;
    int succ=(root+1)%polygon.size();
    int pred=root-1;
    if(pred<0){pred=polygon.size()-1;}

    double angle_target= angle_between_anchored_vectors(polygon[root],polygon[pred],polygon[target]);
    double angle_succ= angle_between_anchored_vectors(polygon[root],polygon[pred],polygon[succ]);

    if(angle_target<angle_succ){
        return true;
    }
    else{
        return false;
    }
}


bool has_edge_angle_inside_polygon(int e_1, int e_2,  const std::vector<Point_s>& polygon){
    int root=e_1;
    int target=e_2;
    int succ=(root+1)%polygon.size();
    int pred=root-1;
    if(pred<0){pred=polygon.size()-1;}

    double angle_target= angle_between_anchored_vectors(polygon[pred],polygon[root],polygon[target]);
    double angle_succ= angle_between_anchored_vectors(polygon[pred],polygon[root],polygon[succ]);

    if(angle_target<angle_succ){
        return true;
    }
    else{
        return false;
    }
}

double angle_between_anchored_vectors( Point_s previous, Point_s anchor, Point_s target){
    if(target==anchor){
        return 2*M_PI;
    }
    if(target==previous){
        return 0;
    }
    Point_s b0=anchor;
    Point_s b1=previous-anchor;
    Point_s b2=target-anchor;
    Point_s c1= crossProduct(b0,b1);
    Point_s c2= crossProduct(b0,b2);

    double up= scalarProduct(c1,c2);
    double down= sqrt(normSquared(c1))*sqrt(normSquared(c2));
    double angle=acos(up/down);
    //! this is counter clockwise!!
    if(!is_left_oriented(target,previous,anchor)){
        angle=2*M_PI-angle;
    }
    return angle;
}

double normalizeAngle(double angle) {

    if (angle>2*M_PI) {
        angle=angle-2*M_PI;
    }
    if(angle<0) {
        angle=angle+2*M_PI;
    }

    return angle;
}

double calculateAngle(Point_2 anchorPred, Point_2 anchor, Point_2 v) {
    // Check if v is equal to anchorPred
    if (v == anchorPred) {
        return 2 * M_PI;
    } else {
        // Calculate differences in coordinates
        double targetX = v.x - anchor.x;
        double targetY = v.y - anchor.y;
        double preEdgeX = anchorPred.x - anchor.x;
        double preEdgeY = anchorPred.y - anchor.y;

        // Calculate angles using atan2
        double targetA = atan2(targetY, targetX);
        double preEdgeA = atan2(preEdgeY, preEdgeX);

        // Calculate relative angle and normalize
        targetA -= preEdgeA;
        targetA = normalizeAngle(targetA);

        return targetA;
    }
}























