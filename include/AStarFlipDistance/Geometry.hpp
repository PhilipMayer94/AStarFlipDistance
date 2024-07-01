
#ifndef A_STAR_FOR_FLIPDISTANCE__GEOMETRY_HPP
#define A_STAR_FOR_FLIPDISTANCE__HPP

#include "AStarFlipDistance/BasicDataStructures.hpp"
#include <cmath>
#include <vector>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>


void test();




bool is_quadrailateral_convex(Point_s p1, Point_s p2, Point_s p3, Point_s p4);

bool is_quadrailateral_convex(Point_2 a, Point_2 b, Point_2 c, Point_2 d);


bool is_triangle_empty_as_indices_in_set(int t1,int t2, int t3, std::vector<Point_2> &points);

bool is_triangle_empty(Point_2 t1,Point_2 t2, Point_2 t3, std::vector<Point_2> &points);

bool is_triangle_empty(Point_s t1,Point_s t2, Point_s t3, std::vector<Point_s> &points);

double distance(Point_2 a, Point_2 b);

double distance(Point_s a, Point_s b);

double distance(Point_2 a, std::vector<double> &circle);

double distance(Point_s a, std::vector<double> &circle);

double area_of_triangle(Point_2 p1, Point_2 p2, Point_2 p3);


bool edge_in_polygon(int e_1, int e_2,  const std::vector<Point_2>& polygon);

bool edge_in_polygon(int e_1, int e_2,  const std::vector<Point_s>& polygon);

double angle_between_anchored_vectors(Point_2 anchor,Point_2 main, Point_2 target);

bool has_edge_angle_inside_polygon(int e_1, int e_2,  const std::vector<Point_2>& polygon);

double angle_between_anchored_vectors(Point_s anchor,Point_s main, Point_s target);

bool has_edge_angle_inside_polygon(int e_1, int e_2,  const std::vector<Point_s>& polygon);


std::vector<double> get_interior_angles_of_triangle(Point_s A, Point_s B,Point_s C);

double get_area_of_triangle(Point_s A, Point_s B, Point_s C);

bool is_point_in_circle(Point_2 a, std::vector<double> &circle);

bool is_point_in_circle(Point_s a, std::vector<double> &circle);

bool is_point_in_triangle(Point_2 p, Point_2 t_1, Point_2 t_2, Point_2 t_3);

bool is_point_in_triangle(Point_s p, Point_s t_1, Point_s t_2, Point_s t_3);


std::vector<double> circle_through_three_points(Point_2 a, Point_2 b, Point_2 c);

std::vector<double> circle_through_three_points(Point_s a, Point_s b, Point_s c);

Point_s circumcenter_spherical_triangle(Point_s a, Point_s b, Point_s c);

Point_s circumcenter_3D_triangle(Point_s a, Point_s b, Point_s c);

bool is_left_oriented(Point_2 v, Point_2 e_1, Point_2 e_2);

bool is_left_oriented(Point_s v, Point_s e_1, Point_s e_2);

bool do_arcs_intersect_including_endpoints(Point_2 p_1, Point_2 p_2, Point_2 q_1, Point_2 q_2);

bool do_arcs_intersect(Point_2 p1, Point_2 p2, Point_2 q1, Point_2 q2);

bool do_arcs_intersect(Point_s p1, Point_s q1, Point_s p2, Point_s q2);

bool do_arcs_intersect_straddle(Point_s p1, Point_s q1, Point_s p2, Point_s q2);


Point_s getLineIntersection(Point_s p_1, Point_s p_2, Point_s q1, Point_s q2);


int orientation(Point_2 p, Point_2 q, Point_2 r);

bool onSegment(Point_2 p, Point_2 q, Point_2 r);

double scalarProduct(Point_s p1, Point_s p2);

Point_s crossProduct(Point_s p1, Point_s p2);

/*Point_s addPoints(Point_s a, Point_s b);

Point_s subPoints(Point_s a, Point_s b);*/

double normSquared(Point_s a);

Point_s get_scaled(double scalar, Point_s point);

double determinante(Point_s v1,Point_s v2,Point_s v3);


double normalizeAngle(double angle);

// Calculates the angle between (a, v) and (a, anchorPred)
double calculateAngle(Point_2 anchorPred, Point_2 anchor, Point_2 v);





#endif //HIGHERORDERDELAUNAY_GEOMETRY_HPP
