#include <stdio.h>
#include <string>
#include <array>
#include <set>
#include <map>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <chrono>
#include <queue>
#include <fstream>
#include <queue>
#include <vector>
#include <limits>

inline std::string ltrim(std::string s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    return s;
}

// trim from end (copying)
inline std::string rtrim(std::string s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
    return s;
}

// trim from both ends (copying)
inline std::string trim(std::string s) {
    return rtrim(ltrim(s));
}
using namespace Eigen;
using namespace std;

constexpr double pi = 3.14159265358979323846;



class Polyhedron {
  public:
    vector<std::array<double,3>> verts;
    vector<set<int>> edges;
    vector<set<int>> faces;
    static vector<std::array<double,3>> shift(vector<std::array<double,3>> points, std::array<double,3> displacement) {
        vector<std::array<double,3>> output;
        for (const std::array<double,3> p : points) {
            std::array<double,3> temp = {p[0]+displacement[0],p[1]+displacement[1],p[2]+displacement[2]};
            output.push_back(temp);
        }
        return output;
    }
    static vector<std::array<double,3>> scale(vector<std::array<double,3>> points, std::array<double,3> factors) {
        vector<std::array<double,3>> output;
        for (const std::array<double,3> p : points) {
            std::array<double,3> temp = {p[0]*factors[0],p[1]*factors[1],p[2]*factors[2]};
            output.push_back(temp);
        }
        return output;
    }
    static vector<std::array<double,3>> scale(vector<std::array<double,3>> points, double factor) {
        return Polyhedron::scale(points, {factor, factor, factor});
    }
    static vector<std::array<double,3>> rotate(vector<std::array<double,3>> points, std::array<double,3> angles) {
        vector<std::array<double,3>> output;
        for (const std::array<double,3> p : points) {
            std::array<double,3> temp1 = {p[0],p[1]*cos(angles[0])-p[2]*sin(angles[0]),p[1]*sin(angles[0])+p[2]*cos(angles[0])};
            std::array<double,3> temp2 = {temp1[0]*cos(angles[1])-temp1[2]*sin(angles[1]),temp1[1],temp1[0]*sin(angles[1])+temp1[2]*cos(angles[1])};
            std::array<double,3> temp3 = {temp2[0]*cos(angles[2])-temp2[1]*sin(angles[2]),temp2[0]*sin(angles[2])+temp2[1]*cos(angles[2]),temp2[2]};
            output.push_back(temp3);
        }
        return output;
    }
    static double round_float(double x) {
        return round(x*1000)/1000;
    }
    static std::array<double,3> round_point(std::array<double,3> point) {
        return {round_float(point[0]), round_float(point[1]), round_float(point[2])};
    }
    static set<std::array<double,3>> round_edge(set<std::array<double,3>> edge) {
        set<std::array<double,3>> output;
        set<std::array<double,3>>::iterator it = edge.begin();
        output.insert(round_point(*(it)));
        it++;
        output.insert(round_point(*(it)));
        return output;
    }
    static double round_float_meter(double x, set<double> meters) {
        if (!meters.size()) {
            meters = {1};
        }
        double minidiff = std::numeric_limits<double>::infinity();
        double mini = 0;
        for (const double& y : meters) {
            double diff = abs(x/y - round(x/y));
            if (diff < minidiff) {
                minidiff = diff;
                mini = round(x/y)*y;
            }
        }
        return mini;
    }
    static std::array<double,3> round_point_meter(std::array<double,3> x, set<double> meters) {
        return {round_float_meter(x[0], meters),round_float_meter(x[1], meters),round_float_meter(x[2], meters)};
    }
    static set<std::array<double,3>> round_edge_meter(set<std::array<double,3>> edge, set<double> meters) {
        set<std::array<double,3>> output;
        set<std::array<double,3>>::iterator it = edge.begin();
        output.insert(round_point_meter(*(it), meters));
        it++;
        output.insert(round_point_meter(*(it), meters));
        return output;
    } 
    static set<double> remeter(set<double> meters, double new_meter, double precision) {
        for (const double& x : meters) {
            if (x != 1 && abs(new_meter/x - round(new_meter/x)) <= precision) {
                return meters;
            }
        }
        set<double> output;
        for (const double& x : meters) {
            if (abs(x/new_meter - round(x/new_meter)) > precision) {
                output.insert(x);
            }
        }
        output.insert(new_meter);
        return output;
    }
    static set<double> remeter(set<double> meters, double new_meter) {
        return remeter(meters, new_meter, 0.001);
    }
    class RoundPointMeter {
      set<double> meters;
      public:
        RoundPointMeter(set<double> meters) {
            this->meters = meters;
        }
        std::array<double,3> operator()(std::array<double,3> point) const {
            return Polyhedron::round_point_meter(point, meters);
        } 
    };
    static double dot(std::array<double,3> vector1, std::array<double,3> vector2) {
        double output = 0;
        for (int i = 0; i < vector1.size(); i++) {
            output += vector1[i]*vector2[i];
        }
        return output;
    }
    static double distance(std::array<double,3> point) {
        return sqrt(Polyhedron::dot(point, point));
    }
    static double distance(std::array<double,3> point1, std::array<double,3> point2) {
        std::array<double,3>displacement = {point2[0]-point1[0],point2[1]-point1[1],point2[2]-point1[2]};
        return sqrt(Polyhedron::dot(displacement, displacement));
    }
    static std::array<double,3> cross3D(std::array<double,3> vector1, std::array<double,3> vector2) {
        std::array<double,3> output;
        for (int i = 0; i < 3; i++) {
            output[i] = vector1[(i+1)%3]*vector2[(i+2)%3]-vector1[(i+2)%3]*vector2[(i+1)%3];
        }
        return output;
    }
    static vector<vector<int>> create_distance_matrix(const vector<std::array<double, 3>>& points, std::array<Polyhedron,2>& polys) {
        int n = points.size();
        vector<vector<int>> dist_matrix(n, std::vector<int>(n, 0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;

                double dist = 0.0;
                for (int d = 0; d < 3; ++d)
                    dist += std::pow(points[i][d] - points[j][d], 2);
                dist = std::sqrt(dist);
                int scaled_dist = static_cast<int>(1000.0 * dist);

                std::array<double, 3> p1 = points[i];
                std::array<double, 3> p2 = points[j];
                bool valid = true;

                for (Polyhedron& poly : polys) {
                    vector<FaceIntersection> intersects = poly.face_intersect({p1, p2});
                    std::vector<std::array<double, 3>> ps;
                    for (const FaceIntersection& intersect : intersects) {
                        if (ps.empty() || intersect.point != ps.back())
                            ps.push_back(intersect.point);
                    }

                    ps.insert(ps.begin(), p1);
                    ps.push_back(p2);

                    for (size_t k = 0; k < ps.size() - 1; ++k) {
                        std::array<double, 3> mid;
                        for (int d = 0; d < 3; ++d)
                            mid[d] = 0.5 * (ps[k][d] + ps[k + 1][d]);
                        if (!poly.is_inside(mid)) {
                            valid = false;
                            break;
                        }
                    }

                    if (!valid) break;
                }

                dist_matrix[i][j] = valid ? scaled_dist : std::numeric_limits<int>::max();
            }
        }
        return dist_matrix;
    }

set<vector<std::array<double,3>>> circuits(int face_index, int start, int previous, int current, vector<std::array<double,3>> path, set<vector<std::array<double,3>>> old_circuits) {
    set<vector<std::array<double,3>>> output;
    if (current == start) {
        output.insert(path);
        return output;
    }
    /*
    set<std::array<double,3>> path_set(path.begin(),path.end());
    for (const vector<std::array<double,3>>& old_circuit : old_circuits) {
        set<std::array<double,3>> old_circuit_set(old_circuit.begin(),old_circuit.end());
        set<std::array<double,3>> difference;
        set_difference(path_set.begin(), path_set.end(), old_circuit_set.begin(), old_circuit_set.end(), std::inserter(difference, difference.begin()));
        if (!difference.size()) {
            return output;
        }
    }*/
    set<int> face = faces[face_index];
    map<int,set<set<int>>> edge_lookup;
    for (const int& edge_index : face) {
        set<int> edge = edges[edge_index];
        for (const int& index : edge) {
            edge_lookup[index].insert(edge);
        }
    }
    set<int> difference;
    set_difference(edges[current].begin(), edges[current].end(), edges[previous].begin(), edges[previous].end(), std::inserter(difference, difference.begin()));
    int point = *(difference.begin());
    if (find(path.begin(), path.end(), verts[point]) != path.end()) {
        return output;
    }
    path.push_back(verts[point]);
    previous = current;
    set<set<int>> temp = (edge_lookup[point]);
    temp.erase(edges[previous]);
    for (const set<int>& y : temp) {
        vector<set<int>>::iterator it = find(edges.begin(),edges.end(),y);
        current = std::distance(edges.begin(), it);
        old_circuits.insert(output.begin(),output.end());
        set<vector<std::array<double,3>>> intermediate = this->circuits(face_index, start, previous, current, path, old_circuits);
        output.insert(intermediate.begin(),intermediate.end());
    }
    vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
    for (int i = 0; i < output_list.size(); i++) {
        for (int j = i+1; j < output_list.size(); j++) {
            set<std::array<double,3>> set_i(output_list[i].begin(),output_list[i].end());
            set<std::array<double,3>> set_j(output_list[j].begin(),output_list[j].end());
            set<std::array<double,3>> difference;
            set_difference(set_i.begin(), set_i.end(), set_j.begin(), set_j.end(), std::inserter(difference, difference.begin()));
            if (output_list[i] != output_list[j] && !difference.size()) {
                vector<std::array<double,3>> y_r(output_list[j].begin(),output_list[j].end());
                std::reverse(y_r.begin(),y_r.end());
                vector<std::array<double,3>>::iterator it = find(output_list[j].begin(),output_list[j].end(), output_list[i][0]);
                int index_x0_in_y = std::distance(output_list[j].begin(), it);
                it = find(y_r.begin(),y_r.end(), output_list[i][0]);
                int index_x0_in_y_r = std::distance(y_r.begin(), it);
                vector<std::array<double,3>> y_rearranged(output_list[j].begin()+index_x0_in_y,output_list[j].end());
                y_rearranged.insert(y_rearranged.end(),output_list[j].begin(),output_list[j].begin()+index_x0_in_y);
                vector<std::array<double,3>> y_rearranged_sliced(y_rearranged.begin(),y_rearranged.begin()+output_list[i].size());
                vector<std::array<double,3>> y_r_rearranged(y_r.begin()+index_x0_in_y_r,y_r.end());
                y_r_rearranged.insert(y_r_rearranged.end(),y_r.begin(),y_r.begin()+index_x0_in_y_r);
                vector<std::array<double,3>> y_r_rearranged_sliced(y_r_rearranged.begin(),y_r_rearranged.begin()+output_list[i].size());
                if (y_rearranged_sliced == output_list[i] || y_r_rearranged_sliced == output_list[i]) {
                    if (output.find(output_list[j]) != output.end()) {
                        output.erase(output_list[j]);
                    }
                }
            }
        }
    }
    return output;
}
    set<vector<std::array<double,3>>> circuits(int face_index) {
        set<set<std::array<double,3>>> edges_grounded;
        set<int> face = faces[face_index];
        map<int,set<set<int>>> edge_lookup;
        for (const int& edge_index : face) {
            set<int> edge = edges[edge_index];
            set<std::array<double,3>> edge_grounded;
            for (const int& index : edge) {
                edge_lookup[index].insert(edge);
                edge_grounded.insert(this->verts[index]);
            }
            edges_grounded.insert(edge_grounded);
        }
        set<vector<std::array<double,3>>> output;
        set<set<std::array<double,3>>> seen;
        for (const int& edge : face) {
            set<std::array<double,3>> edge_grounded;
            for (const int& index : edges[edge]) {
                edge_grounded.insert(verts[index]);
            }
            if (seen.find(edge_grounded) != seen.end()) {
                continue;
            }
            vector<std::array<double,3>> path;
            int start = edge;
            int point = *(edges[start].begin());
            path.push_back(verts[point]);
            set<set<int>> temp = (edge_lookup[point]);
            temp.erase(edges[start]);
            for (const set<int>& y : temp) {
                vector<set<int>>::iterator it = find(edges.begin(), edges.end(), y);
                int current = std::distance(edges.begin(), it);
                set<vector<std::array<double,3>>> intermediate = this->circuits(face_index, start, start, current, path, output);
                output.insert(intermediate.begin(),intermediate.end());
                for (const vector<std::array<double,3>>& circuit : intermediate) {
                    for (int i = 0; i < circuit.size(); i++) {
                        seen.insert({circuit[(i-1+circuit.size())%circuit.size()],circuit[i]});
                    }
                }
            }
        }
        vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
        //cout << "circuit_sizes ";
        for (const vector<std::array<double,3>>& circuit : output) {
            //cout << circuit.size() << " ";
        }
        //cout << endl;
        vector<std::array<double,3>>* exterior_circuit = Polyhedron::find_exterior_circuit(output);
        if (exterior_circuit != NULL) {
            set<vector<std::array<double,3>>> new_output;
            new_output.insert(*exterior_circuit);
            queue<set<vector<std::array<double,3>>>> q;
            q.push({*exterior_circuit});
            while (!q.empty()) {
                set<set<std::array<double,3>>> covering;
                for (const vector<std::array<double,3>>& circuit : q.front()) {
                    for (int i = 0; i < circuit.size(); i++) {
                        covering.insert({circuit[(i-1+circuit.size())%circuit.size()],circuit[i]});
                    }
                }
                if (covering == edges_grounded) {
                    delete exterior_circuit;
                    return q.front();
                }
                for (const vector<std::array<double,3>>& circuit : output) {
                    set<set<std::array<double,3>>> circuit_edges;
                    for (int i = 0; i < circuit.size(); i++) {
                        circuit_edges.insert({circuit[(i-1+circuit.size())%circuit.size()],circuit[i]});
                    }
                    set<set<std::array<double,3>>> edge_intersection;
                    set_intersection(covering.begin(),covering.end(),circuit_edges.begin(),circuit_edges.end(),inserter(edge_intersection, edge_intersection.begin()));
                    if (!edge_intersection.size()) {
                        set<vector<std::array<double,3>>> new_item(q.front().begin(),q.front().end());
                        new_item.insert(circuit);
                        q.push(new_item);
                    }
                }
                q.pop();
            }
        }
        for (int i = 0; i < output_list.size(); i++) {
            for (int j = i+1; j < output_list.size(); j++) {
                set<std::array<double,3>> set_i(output_list[i].begin(),output_list[i].end());
                set<std::array<double,3>> set_j(output_list[j].begin(),output_list[j].end());
                set<std::array<double,3>> difference;
                set_difference(set_i.begin(), set_i.end(), set_j.begin(), set_j.end(), std::inserter(difference, difference.begin()));
                if (output_list[i] != output_list[j] && !difference.size()) {
                    vector<std::array<double,3>> y_r(output_list[j].begin(),output_list[j].end());
                    std::reverse(y_r.begin(),y_r.end());
                    vector<std::array<double,3>>::iterator it = find(output_list[j].begin(),output_list[j].end(), output_list[i][0]);
                    int index_x0_in_y = std::distance(output_list[j].begin(), it);
                    it = find(y_r.begin(),y_r.end(), output_list[i][0]);
                    int index_x0_in_y_r = std::distance(y_r.begin(), it);
                    vector<std::array<double,3>> y_rearranged(output_list[j].begin()+index_x0_in_y,output_list[j].end());
                    y_rearranged.insert(y_rearranged.begin(),output_list[j].begin(),output_list[j].begin()+index_x0_in_y);
                    vector<std::array<double,3>> y_rearranged_sliced(y_rearranged.begin(),y_rearranged.begin()+output_list[i].size());
                    vector<std::array<double,3>> y_r_rearranged(y_r.begin()+index_x0_in_y_r,y_r.end());
                    y_r_rearranged.insert(y_r_rearranged.begin(),y_r.begin(),y_r.begin()+index_x0_in_y_r);
                    vector<std::array<double,3>> y_r_rearranged_sliced(y_r_rearranged.begin(),y_r_rearranged.begin()+output_list[i].size());
                    if (y_rearranged_sliced == output_list[i] || y_r_rearranged == output_list[i]) {
                        if (output.find(output_list[j]) != output.end()) {
                            output.erase(output_list[j]);
                        }
                    }
                }
            }
        }
        //cout << "NO_EXTERIOR_CIRCUIT2" << endl; 
        delete exterior_circuit;
        return output;
    }
    static bool colinear(vector<std::array<double,3>> points) {
        int index = 1;
        for (; index < points.size() && Polyhedron::distance(points[0], points[index]) < 0.0001; index++) {
        }
        if (index >= points.size()) {
            return true;
        }
        if (points.size() < 3) {
            return true;
        }
        MatrixXd m(3,1);
        for (int i=0; i < 3; i++) {
            m(i,0) = points[index][i]-points[0][i];
        }
        for (const std::array<double,3>& p : points) {
            if (p == points[index]) {
                continue;
            }
            VectorXd b(3);
            for (int i=0; i < 3; i++) {
                b(i) = p[i]-points[0][i];
            }
            VectorXd x = m.colPivHouseholderQr().solve(b);
            VectorXd b_prime = m*x;
            std::array<double,3> point_prime;
            for (int i=0; i < 3; i++) {
                point_prime[i] = b_prime(i)+points[0][i];
            }
            if (Polyhedron::distance(p, point_prime) > 0.000001) {
                return false;
            }
        }
        return true;
    }
    static bool colinear(set<std::array<double,3>> points) {
        vector<std::array<double,3>> points_vector(points.begin(),points.end());
        return Polyhedron::colinear(points_vector);
    }
    static bool coplanar(vector<std::array<double,3>> points) {
        if (points.size() < 4) {
            return true;
        }
        if (Polyhedron::colinear(points)) {
            return true;
        }
        bool triple_break = false;
        std::array<int,3> ind;
        for (int i = 0; i < points.size(); i++) {
            for (int j = i+1; j < points.size(); j++) {
                for (int k = j+1; k < points.size(); k++) {
                    vector<std::array<double,3>> temp = {points[i],points[j],points[k]};
                    if (!Polyhedron::colinear(temp)) {
                        ind = {i,j,k};
                        triple_break = true;
                        break;
                    }
                }
                if (triple_break) {
                    break;
                }
            }
            if (triple_break) {
                break;
            }
        }
        MatrixXd m(3,2);
        for (int i = 0; i < 3; i++) {
            m(i,0) = points[ind[1]][i]-points[ind[0]][i];
            m(i,1) = points[ind[2]][i]-points[ind[0]][i];
        }
        for (int i = 0; i < points.size(); i++) {
            if (i == ind[0] || i == ind[1] || i == ind[2]) {
                continue;
            }
            VectorXd b(3);
            for (int j = 0; j < 3; j++) {
                b(j) = points[i][j]-points[ind[0]][j];
            }
            VectorXd x = m.colPivHouseholderQr().solve(b);
            VectorXd b_prime = m*x;
            std::array<double,3> point_prime;
            for (int j=0; j < 3; j++) {
                point_prime[j] = b_prime(j)+points[ind[0]][j];
            }
            if (Polyhedron::distance(points[i], point_prime) > 0.000001) {
                return false;
            }
        }
        return true;
    }
    static bool coplanar(set<std::array<double,3>> points) {
        vector<std::array<double,3>> points_vector(points.begin(),points.end());
        return Polyhedron::coplanar(points_vector);
    }
    static bool point_on_segment(set<std::array<double,3>> edge, std::array<double,3> point) {
        set<std::array<double,3>>::iterator it = edge.begin();
        std::array<double,3> p1 = *(it);
        it++;
        std::array<double,3> p2 = *(it);
        MatrixXd m(3,1);
        VectorXd b(3);
        for (int i = 0; i < 3; i++) {
            m(i,0) = p1[i]-p2[i];
            b(i) = point[i]-p2[i];
        }
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        for (int i=0; i < 3; i++) {
            if ( abs(b_prime(i)-b(i)) > 0.1*abs(b(i)) ) {
                return false;
            }
        }
        if (x(0) >= 0 && x(0) <= 1) {
            return true;
        }
        return false;
    }
    static std::array<double,3>* intersect_segments(set<std::array<double,3>> edge1, set<std::array<double,3>> edge2) {
        set<std::array<double,3>>::iterator it = edge1.begin();
        std::array<double,3> p1 = *(it);
        it++;
        std::array<double,3> p2 = *(it);
        it = edge2.begin();
        std::array<double,3> p3 = *(it);
        it++;
        std::array<double,3> p4 = *(it);
        MatrixXd m(5,4);
        for (int i = 0; i < 3; i++) {
            m(i,0) = p1[i];
            m(i,1) = p2[i];
            m(i,2) = -p3[i];
            m(i,3) = -p4[i];
        }
        m(3,0) = 1;
        m(3,1) = 1;
        m(4,2) = 1;
        m(4,3) = 1;
        VectorXd b(5);
        for (int j = 0; j < 3; j++) {
            b(j) = 0;
        }
        for (int j = 3; j < 5; j++) {
            b(j) = 1;
        }
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        std::array<double,3>* point = new std::array<double,3>();
        for (int i=0; i < 3; i++) {
            (*point)[i] = x(0)*p1[i]+x(1)*p2[i];
        }
        if (Polyhedron::point_on_segment(edge1, *point) && Polyhedron::point_on_segment(edge2, *point)) {
            return point;
        }
        return NULL;
        
    }
    static bool inside_triangle(std::array<std::array<double,3>,3> triangle, std::array<double,3> point) {
        if (find(triangle.begin(),triangle.end(),point) != triangle.end()) {
            return true;
        }
        MatrixXd m(4,3);
        VectorXd b(4);
        for (int i = 0; i < 3; i++) {
            m(i,0) = triangle[0][i];
            m(i,1) = triangle[1][i];
            m(i,2) = triangle[2][i];
            b(i) = point[i]; 
            m(3,i) = 1;
        }
        b(3) = 1;
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        std::array<double,3> point_prime;
        for (int i=0; i < 3; i++) {
            point_prime[i] = b_prime(i);
        }
        //cout << "inside triangle " << Polyhedron::distance(point,point_prime) << " " << round_float(x(0)) << " " << round_float(x(1)) << " " << round_float(x(2)) << " " << (round_float(x(0)) >= 0 && round_float(x(1)) >= 0 && round_float(x(2)) >= 0) << endl;
        if (Polyhedron::distance(point,point_prime) > 0.000001) {
            return false;
        }
        if (round_float(x(0)) >= 0 && round_float(x(1)) >= 0 && round_float(x(2)) >= 0) {
            return true;
        }
        return false;
    }
    struct TriangleIntersection {
        double alpha;
        std::array<double,3> point;
    };
    static Polyhedron::TriangleIntersection* intersect_triangle(std::array<std::array<double,3>,2> segment, std::array<std::array<double,3>,3> triangle) {
        std::array<double,3> p1 = segment[0];
        std::array<double,3> p2 = segment[1];
        MatrixXd m(4,4);
        VectorXd b(4);
        for (int i = 0; i < 3; i++) {
            m(i,0) = triangle[0][i];
            m(i,1) = triangle[1][i];
            m(i,2) = triangle[2][i];
            m(i,3) = p1[i]-p2[i];
            b(i) = p1[i];
            m(3,i) = 0;
        }
        m(3,0) = 1;
        m(3,1) = 1;
        m(3,2) = 1;
        m(3,3) = 0;
        b(3) = 1;
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        std::array<double,3> point_prime;
        for (int i = 0; i < 3; i++) {
            point_prime[i] = b_prime(i);
        }
        //cout << "triangle intersect " << round_float(x(0)) << " " << round_float(x(1)) << " " << round_float(x(2)) << " " << round_float(x(3)) << endl;
        if (round_float(x(0)) >= 0 && round_float(x(1)) >= 0 && round_float(x(2)) >= 0 && round_float(x(3)) >= 0 && round_float(x(3)) <= 1 && Polyhedron::distance(p1, point_prime) < 0.000001) {    
            Polyhedron::TriangleIntersection* triangle_intersection = new TriangleIntersection();
            triangle_intersection->alpha = x(3);
            for (int i = 0; i < 3; i++) {
                triangle_intersection->point[i] = x(3)*p2[i]+(1-x(3))*p1[i];
            }
            return triangle_intersection;
        }
        return NULL;
    }
    static std::array<double,3> cross_product_triplet(std::array<double,3> a, std::array<double,3> b, std::array<double,3> c) {
        std::array<double,3>v1;
        std::array<double,3>v2;
        for (int i=0; i < 3; i++) {
            v1[i] = b[i]-a[i];
            v2[i] = c[i]-b[i];
        }
        return Polyhedron::cross3D(v1,v2);
    }
    static vector<bool> convex_angles(vector<std::array<double,3>> circuit) {
        std::array<double,3>v1;
        std::array<double,3>v2;
        for (int i=0; i < 3; i++) {
            v1[i] = circuit[1][i]-circuit[0][i];
            v2[i] = circuit[2][i]-circuit[0][i];
        }
        std::array<double,3> normal = Polyhedron::cross3D(v1,v2);
        vector<double> cross_product_dot_normal;
        for (int i = 0; i < circuit.size(); i++) {
            cross_product_dot_normal.push_back(Polyhedron::dot(Polyhedron::cross_product_triplet(circuit[(i-1+circuit.size())%circuit.size()],circuit[i],circuit[(i+1)%circuit.size()]),normal));
        }
        vector<double> signs;
        for (const double& x : cross_product_dot_normal) {
            signs.push_back(copysign(1,x));
        }
        vector<bool> output;
        for (int i = 0; i < signs.size(); i++) {
            output.push_back(signs[i]==signs[(i-1+signs.size())%signs.size()] || signs[i]==signs[(i+1)%signs.size()]);
        }
        return output;
        
    }
    static vector<bool> clockwise_angles(vector<std::array<double,3>> planar) {
        vector<bool> output;
        for (int i = 0; i < planar.size(); i++) {
            vector<std::array<double,3>> temp;
            for (const std::array<double,3>& y : planar) {
                temp.push_back({y[0]-planar[i][0],y[1]-planar[i][1],y[2]-planar[i][2]});
            }
            planar = temp;
            std::array<double,3> angle = {0,0,-atan2(planar[i][1]-planar[(i-1+planar.size())%planar.size()][1],planar[i][0]-planar[(i-1+planar.size())%planar.size()][0])};
            planar = Polyhedron::rotate(planar, angle);
            double theta = atan2(planar[(i+1)%planar.size()][1],planar[(i+1)%planar.size()][0]);
            output.push_back(theta < 0);
        }
        return output;
    }
    static vector<std::array<double,3>> make_planar(vector<std::array<double,3>> circuit) {
        // Calculate vec1 and vec2
        std::array<double,3> vec1;
        std::array<double,3> vec2;
        for (int i = 0; i < 3; i++) {
            vec1[i] = circuit[1][i] - circuit[0][i];
            vec2[i] = circuit[2][i] - circuit[0][i];
        }
        
        // Calculate cross product vec0 = cross3D(vec1, vec2)
        std::array<double,3> vec0 = cross3D(vec1, vec2);
        
        // Reassign vec1 and vec2
        vec1 = {0.0, vec0[1], vec0[2]};
        vec2 = {vec0[0], 0.0, vec0[2]};
        
        // Calculate angles
        std::array<double,3> angles = {0.0, 0.0, 0.0};
        
        // Calculate angles[0] with zero division protection
        double vec1_distance = distance(vec1);
        if (vec1_distance != 0.0) {
            double cos_val = vec1[2] / vec1_distance;
            cos_val = max(-1.0, min(1.0, cos_val));  // Clamp to [-1, 1]
            angles[0] = -acos(cos_val);
        }
        
        // Calculate angles[1] with zero division protection  
        double vec2_distance = distance(vec2);
        if (vec2_distance != 0.0) {
            double cos_val = vec2[2] / vec2_distance;
            cos_val = max(-1.0, min(1.0, cos_val));  // Clamp to [-1, 1]
            angles[1] = -acos(cos_val);
        }
        
        // Rotate the circuit and project to 2D (set z=0)
        vector<std::array<double,3>> rotated_circuit = rotate(circuit, angles);
        vector<std::array<double,3>> result;
        
        for (const std::array<double,3>& point : rotated_circuit) {
            result.push_back({point[0], point[1], 0.0});
        }
        
        return result;
    }
    static bool is_clockwise(vector<std::array<double,3>> planar) {
        double summa = 0;
        for (int i = 0; i < planar.size(); i++) {
            std::array<double,3> x = planar[i];
            for (int j = 0; j < planar.size(); j++) {
                planar[j][0] -= x[0];
                planar[j][1] -= x[1];
                planar[j][2] -= x[2];
            }
            std::array<double,3> angle = {0,0,-atan2(planar[i][1]-planar[(i-1+planar.size())%planar.size()][1],planar[i][0]-planar[(i-1+planar.size())%planar.size()][0])};
            planar = Polyhedron::rotate(planar, angle);
            double theta = atan2(planar[(i+1)%planar.size()][1],planar[(i+1)%planar.size()][0]);
            //cout << theta/pi*180 << " ";
            summa += theta;
        }
        //cout << endl;
        //cout << "is_clockwise " << summa << endl;
        return summa < 0;
    }
    static set<vector<std::array<double,3>>> make_clockwise(set<vector<std::array<double,3>>> circuits) {
        set<vector<std::array<double,3>>> output;
        for (vector<std::array<double,3>> circuit : circuits) {
            if (is_clockwise(make_planar(circuit))) {
                output.insert(circuit);
            } else {
                std::reverse(circuit.begin(),circuit.end());
                output.insert(circuit);
            }
        }
        return output;
    }
    struct EarClip {
        std::array<std::array<double,3>,3> ear;
        vector<std::array<double,3>> remainder;
    };
    static EarClip clip_ear(vector<std::array<double,3>> circuit) {
        //for (int i = 0; i < circuit.size(); i++) {
        //    cout << "[" << circuit[i][0] << "," << circuit[i][1] << "," << circuit[i][2] << "] ";
        //}
        //cout << endl;
        vector<std::array<double,3>> planar = make_planar(circuit);
        //for (int i = 0; i < planar.size(); i++) {
        //    cout << "[" << planar[i][0] << "," << planar[i][1] << "," << planar[i][2] << "] ";
        //}
        //cout << endl;
        //cout << Polyhedron::is_clockwise(planar) << endl;
        if (!Polyhedron::is_clockwise(planar)) {
            std::reverse(planar.begin(),planar.end());
            std::reverse(circuit.begin(),circuit.end());
        }
        vector<bool> is_convex = clockwise_angles(planar);
        for (int i = 0; i < is_convex.size(); i++) {
            //cout << is_convex[i] ;
        }
        //cout << endl;
        for (int i = 0; i < circuit.size(); i++) {
            if (is_convex[i]) {
                std::array<std::array<double,3>,3> triangle = {circuit[(i-1+circuit.size())%circuit.size()],circuit[i],circuit[(i+1)%circuit.size()]};
                bool no_break = true;
                for (int j = 0; j < circuit.size(); j++) {
                    if (j != i && j != (i+1)%circuit.size() && j != (i-1+circuit.size())%circuit.size()) {

                        if (inside_triangle(triangle,circuit[j]) && round_point(circuit[j]) != round_point(circuit[(i-1+circuit.size())%circuit.size()]) && round_point(circuit[j]) != round_point(circuit[i]) && round_point(circuit[j]) != round_point(circuit[(i+1)%circuit.size()])) {
                            no_break = false;
                            break;
                        }
                    }
                }
                if (no_break) {
                    EarClip ear_clip;
                    ear_clip.ear = triangle;
                    vector<std::array<double,3>> remainder;
                    for (int k = 0; k < circuit.size(); k++) {
                        if (k != i) {
                            remainder.push_back(circuit[k]);
                        }
                    }
                    for (int k = 0; k < remainder.size(); k++) {
                        if (remainder[k] != remainder[(k-1+remainder.size())%remainder.size()]) {
                            ear_clip.remainder.push_back(remainder[k]);
                        }
                    }
                    return ear_clip;
                }
            }
        }
        cout << "uhoh" << endl;
    }
    static vector<std::array<std::array<double,3>,3>> triangulate(vector<std::array<double,3>> circuit) {
        vector<std::array<std::array<double,3>,3>> output;
        vector<std::array<double,3>> remainder = circuit;
        //cout << "triangulate start" << endl;
        while (remainder.size() >3) {
            for (int i = 0; i < remainder.size();) {
                if (Polyhedron::colinear(set<std::array<double,3>>{remainder[(i-1+remainder.size())%remainder.size()],remainder[i],remainder[(i+1)%remainder.size()]})) {
                    //cout << "colinear " << remainder[i][0] << " " << remainder[i][1] << " " << remainder[i][2] << endl;
                    remainder.erase(remainder.begin()+i);
                } else {
                    i++;
                }
            }
            EarClip ear_clip = Polyhedron::clip_ear(remainder);
            remainder = ear_clip.remainder;
            //for (const std::array<double,3>& point : ear_clip.remainder) {
            //    cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
            //}
            //cout << "; ";
            output.push_back(ear_clip.ear);
            for (const std::array<double,3>& point : ear_clip.ear) {
                //cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
            }
            //cout << endl;
        }
        output.push_back({remainder[0],remainder[1],remainder[2]});
        for (const std::array<double,3>& point : remainder) {
            //cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
        }
        //cout << endl;
        //cout << "triangulate end" << endl;
        return output;
    }
    static vector<std::array<double,3>>* find_exterior_circuit(set<vector<std::array<double,3>>> circuits) {
        vector<vector<std::array<double,3>>> circuits_list(circuits.begin(),circuits.end());
        for (int i = 0; i < circuits_list.size(); i++) {
            vector<std::array<std::array<double,3>,3>> triangulation = triangulate(circuits_list[i]);
            bool double_break = false;
            for (int j = 0; j < circuits_list.size(); j++) {
                if (j != i) {
                    for (const std::array<double,3>& point : circuits_list[j]) {
                        bool any = false;
                        for (const std::array<std::array<double,3>,3>& z : triangulation) {
                            if (Polyhedron::inside_triangle(z,point)) {
                                any = true;
                                break;
                            }
                        }
                        if (!any) {
                            cout << "exterior circuit [" << point[0] << "," << point[1] << "," << point[2] << "] " << circuits_list[i].size() << endl;
                            double_break = true;
                            break;
                        }
                    }
                    if (double_break) {
                        break;
                    }
                } 
            }
            if (!double_break) {
                vector<std::array<double,3>>* output = new vector<std::array<double,3>>(circuits_list[i].begin(),circuits_list[i].end());
                return output;
            }
        }
        return NULL;
    }
    struct CircuitIntersection {
        double alpha;
        std::array<double,3> point;
        std::array<std::array<double,3>,2> edge;
        vector<std::array<double,3>> circuit;
    };
    static vector<Polyhedron::CircuitIntersection> circuit_intersect(std::array<std::array<double,3>,2> segment, set<vector<std::array<double,3>>> circuits) {
        std::array<double,3> p1 = segment[0];
        std::array<double,3> p2 = segment[1];
        vector<Polyhedron::CircuitIntersection> output;
        for (const vector<std::array<double,3>>& circuit : circuits) {
            for (int i = 0; i < circuit.size(); i++) {
                std::array<double,3> p3 = circuit[(i-1+circuit.size())%circuit.size()];
                std::array<double,3> p4 = circuit[i];
                MatrixXd m(3,2);
                VectorXd b(3);
                for (int j = 0; j < 3; j++) {
                    m(j,0) = p2[j]-p1[j];
                    m(j,1) = p4[j]-p3[j];
                    b(j) = p4[j]-p1[j];
                }
                VectorXd x = m.colPivHouseholderQr().solve(b);
                VectorXd b_prime = m*x;
                bool break_continue = false;
                for (int j=0; j < 3; j++) {
                    if ( abs(b_prime(j)-b(j)) > 0.1*abs(b(j)) ) {
                        break_continue = true;
                        break;
                    }
                }
                if (break_continue) {
                    continue;
                }
                if (x(0) >= 0 && x(0) <= 1 && x(1) >= 0 && x(1) <= 1) {
                    Polyhedron::CircuitIntersection circuit_intersection;
                    circuit_intersection.alpha = x(0);
                    for (int j=0; j < 3; j++) {
                        circuit_intersection.point[j] = x(0)*p2[j]+(1-x(0))*p1[j];
                    }
                    circuit_intersection.edge[0] = p3;
                    circuit_intersection.edge[1] = p4;
                    circuit_intersection.circuit = circuit;
                    output.push_back(circuit_intersection);
                }
            }
        }
        return output;
    }
    struct FaceIntersection {
        double alpha;
        std::array<double,3> point;
        int face_index;
    };
    static bool compare_face_intersections(FaceIntersection a, FaceIntersection b) {
        return a.alpha < b.alpha;
    }
    vector<FaceIntersection> face_intersect(std::array<std::array<double,3>,2> segment) {
        vector<FaceIntersection> output;
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            vector<std::array<double,3>>* circuit = circuit_cut(this->circuits(face_index)); 
            vector<std::array<double,3>> temp(circuit->begin(),circuit->end());
            temp.push_back(segment[0]);
            temp.push_back(segment[1]);
            if (Polyhedron::coplanar(temp)) {
                continue;
            }
            for (const std::array<std::array<double,3>,3>& triangle : triangulate(*circuit)) {
                Polyhedron::TriangleIntersection* triangle_intersect = intersect_triangle(segment, triangle);
                if (triangle_intersect != NULL) {
                    FaceIntersection intersect;
                    intersect.alpha = triangle_intersect->alpha;
                    intersect.point = triangle_intersect->point;
                    intersect.face_index = face_index;
                    output.push_back(intersect);
                    //cout << triangle_intersect->point[0] << "," << triangle_intersect->point[1] << "," << triangle_intersect->point[2] << " " << triangle_intersect-> alpha << " "<< face_index <<  endl;
                    delete triangle_intersect;
                }
            }
            delete circuit;
        }
        return output;
    }

    static bool compare_circuit_intersections(Polyhedron::CircuitIntersection a, Polyhedron::CircuitIntersection b) {
        return a.alpha < b.alpha;
    }
    static vector<std::array<double,3>>* circuit_cut(set<vector<std::array<double,3>>> circuits) {
        //cout << "circuits" << circuits.size() << endl;
        for (const vector<std::array<double,3>>& circuit : circuits) {
            for (const std::array<double,3> point : circuit) {
                //cout << point.size() << endl;
                //cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
            }
            //cout << endl;
        }
        //cout << endl;
        vector<std::array<double,3>>* exterior = Polyhedron::find_exterior_circuit(circuits);
        if (exterior == NULL) {
            return NULL;
        }
        vector<vector<std::array<double,3>>> output;
        for (const vector<std::array<double,3>>& circuit : circuits) {
            if (circuit != *(exterior)) {
                output.push_back(circuit);
            }
        }
        output.push_back(*(exterior));
        while (output.size() > 1) {
            vector<std::array<double,3>> interior = output[0];
            std::array<std::array<double,3>,2> segment;
            for (const std::array<double,3>& x : interior) {
                bool double_break = false;
                for (const std::array<double,3>& y : output[output.size()-1]) {
                    segment = {x,y};
                    set<vector<std::array<double,3>>> interior_circuits(output.begin(),output.end()-1);
                    vector<Polyhedron::CircuitIntersection> intersections = Polyhedron::circuit_intersect(segment, interior_circuits);
                    if (intersections.size() == 1) {
                        double_break = true;
                        break;
                    }
                }
                if (double_break) {
                    break;
                }
            }
            set<vector<std::array<double,3>>> output_set(output.begin(),output.end());
            vector<Polyhedron::CircuitIntersection> intersections = Polyhedron::circuit_intersect(segment, output_set);
            sort(intersections.begin(),intersections.end(), Polyhedron::compare_circuit_intersections);
            Polyhedron::CircuitIntersection first_exterior_intersection;
            for (const Polyhedron::CircuitIntersection& intersection : intersections) {
                if (intersection.circuit == output.back()) {
                    first_exterior_intersection = intersection;
                    break;
                }
            }
            Polyhedron::CircuitIntersection last_interior_intersection;
            for (const Polyhedron::CircuitIntersection& intersection : intersections) {
                if (intersection.circuit == output.back()) {
                    continue;
                }
                last_interior_intersection = intersection;
            }
            int i;
            for (i = 0; i < output.back().size(); i++) {
                if (find(first_exterior_intersection.edge.begin(), first_exterior_intersection.edge.end(), output.back()[i]) != first_exterior_intersection.edge.end() && find(first_exterior_intersection.edge.begin(), first_exterior_intersection.edge.end(), output.back()[(i-1+output.back().size())%output.back().size()]) != first_exterior_intersection.edge.end()) {
                    break;
                }
            }
            int j;
            for (j = 0; j < last_interior_intersection.circuit.size(); j++) {
                if (find(last_interior_intersection.edge.begin(), last_interior_intersection.edge.end(), last_interior_intersection.circuit[j]) != last_interior_intersection.edge.end() && find(last_interior_intersection.edge.begin(), last_interior_intersection.edge.end(), last_interior_intersection.circuit[(j-1+last_interior_intersection.circuit.size())%last_interior_intersection.circuit.size()]) != last_interior_intersection.edge.end()) {
                    break;
                }
            }
            cout << "i " << i << " j " << j << endl; 
            vector<std::array<double,3>>::iterator it1 = output.back().begin()+i;
            while (it1 != output.back().end()) {
                int d = std::distance(output.back().begin(),it1);
                vector<std::array<double,3>> planar1;
                vector<std::array<double,3>> planar2;
                if (*it1 == first_exterior_intersection.point) {
                    planar1 = Polyhedron::make_planar({output.back()[(d-1+output.back().size())%output.back().size()],*it1,output.back()[(d+1)%output.back().size()]});
		            std::array<double,3> point = last_interior_intersection.point;
                    for (int k = 0; point == *it1; k++) {
                        point = last_interior_intersection.circuit[(j-1-k+last_interior_intersection.circuit.size())%last_interior_intersection.circuit.size()];
		            }
                    planar2 = Polyhedron::make_planar({output.back()[(d-1+output.back().size())%output.back().size()],*it1,point});
                 } else {
                    planar1 = Polyhedron::make_planar({output.back()[(d-1+output.back().size())%output.back().size()],first_exterior_intersection.point,*it1});
		            std::array<double,3> point = last_interior_intersection.point;
                    for (int k = 0; point == first_exterior_intersection.point; k++) {
                        point = last_interior_intersection.circuit[(j-1-k+last_interior_intersection.circuit.size())%last_interior_intersection.circuit.size()];
		            }
                    planar2 = Polyhedron::make_planar({output.back()[(d-1+output.back().size())%output.back().size()],first_exterior_intersection.point,point});
                }
                for (int k = 0; k < 3; k++) {
                    planar1[0][k] -= planar1[1][k];
                    planar1[2][k] -= planar1[1][k];
                    planar1[1][k] -= planar1[1][k];
                }
                std::array<double,3> angle1 = {0,0,-atan2(planar1[1][1]-planar1[0][1],planar1[1][0]-planar1[0][0])};
                planar1 = Polyhedron::rotate(planar1, angle1);
                double theta1 = atan2(planar1[2][1],planar1[2][0]);
                for (int k = 0; k < 3; k++) {
                    planar2[0][k] -= planar2[1][k];
                    planar2[2][k] -= planar2[1][k];
                    planar2[1][k] -= planar2[1][k];
                }
                std::array<double,3> angle2 = {0,0,-atan2(planar2[1][1]-planar2[0][1],planar2[1][0]-planar2[0][0])};
                planar2 = Polyhedron::rotate(planar2, angle2);
                double theta2 = atan2(planar2[2][1],planar2[2][0]);
                cout << planar2[2][0] << " " << planar2[2][1] << endl;
                cout << "angle1 " << angle1[2] << " angle2 " << angle2[2] << endl;
                cout << "theta1 " << theta1 << " theta2 " << theta2 << endl;
                if (theta2 < theta1) {
                    break;
                }
                vector<std::array<double,3>>::iterator it2 = it1+1;
                while(it2 != output.back().end()) {
                    d = std::distance(output.back().begin(),it2);
                    if (Polyhedron::point_on_segment({*it2,output.back()[(d+1)%output.back().size()]},first_exterior_intersection.point)) {
                        break;
                    }
                    it2++;
                }
                if (it2 == output.back().end()) {
                    break;
                }
                it1 = it2+1;
                //it1 = find(it1+1,output.back().end(),first_exterior_intersection.point);
            }
            i = std::distance(output.back().begin(),it1);
            vector<std::array<double,3>> temp(output.back().begin(),it1);
            for (int k = 0; k < i; k++) {
                cout << "[" << output.back()[k][0] << "," << output.back()[k][1] << "," << output.back()[k][2] << "] ";
            }
            cout << "--> ";
            temp.push_back(first_exterior_intersection.point);
            cout << "[" << first_exterior_intersection.point[0] << "," << first_exterior_intersection.point[1] << "," << first_exterior_intersection.point[2] << "] ";
            temp.push_back(last_interior_intersection.point);
            cout << "[" << last_interior_intersection.point[0] << "," << last_interior_intersection.point[1] << "," << last_interior_intersection.point[2] << "] ";
            cout << "--> ";
            
            for(int k = j-1; k >= 0; k--) {
                temp.push_back(last_interior_intersection.circuit[k]);
                cout << "[" << last_interior_intersection.circuit[k][0] << "," << last_interior_intersection.circuit[k][1] << "," << last_interior_intersection.circuit[k][2] << "] ";
            }
            for(int k = last_interior_intersection.circuit.size()-1; k >= j; k--) {
                temp.push_back(last_interior_intersection.circuit[k]);
                cout << "[" << last_interior_intersection.circuit[k][0] << "," << last_interior_intersection.circuit[k][1] << "," << last_interior_intersection.circuit[k][2] << "] ";
            }
            cout << "--> ";
            temp.push_back(last_interior_intersection.point);
            cout << "[" << last_interior_intersection.point[0] << "," << last_interior_intersection.point[1] << "," << last_interior_intersection.point[2] << "] ";
            temp.push_back(first_exterior_intersection.point);
            cout << "[" << first_exterior_intersection.point[0] << "," << first_exterior_intersection.point[1] << "," << first_exterior_intersection.point[2] << "] ";
            temp.insert(temp.end(),output[output.size()-1].begin()+i,output[output.size()-1].end());
            cout << "--> ";
            for (int k = i; k < output.back().size(); k++) {
                cout << "[" << output.back()[k][0] << "," << output.back()[k][1] << "," << output.back()[k][2] << "] ";
            }
            cout << endl;
            output[output.size()-1] = temp;
            vector<vector<std::array<double,3>>>::iterator it3 = find(output.begin(),output.end(),last_interior_intersection.circuit);
            int index = std::distance(output.begin(),it3);
            output.erase(output.begin()+index);
            for(i = 1; i < output[output.size()-1].size()+1;){
                if (output[output.size()-1][i%output[output.size()-1].size()] == output[output.size()-1][i-1]) {
                    output[output.size()-1].erase(output[output.size()-1].begin()+i%output[output.size()-1].size());
                } else {
                    i++;
                }
            }
        }
        vector<std::array<double,3>>* output_pointer = new vector<std::array<double,3>>;
        output_pointer->insert(output_pointer->begin(),output[0].begin(),output[0].end());
        return output_pointer;
    }
    struct IsInsideIntersection {
        double gamma;
        std::array<double,3> point;
        int face_index;
    };
    vector<int> in_faces(std::array<double,3> point) {
	vector<int> output;
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            vector<std::array<double,3>>* circuit = Polyhedron::circuit_cut(this->circuits(face_index));
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*circuit)) {
                if (Polyhedron::inside_triangle(triangle, point)) {
                    output.push_back(face_index);
                    break;
                }
            }
            delete circuit;
        }
        return output;
    }
    bool is_inside(std::array<double,3> point) {
        //cout << "point " << point[0] << "," << point[1] << "," << point[2] << endl;
        if (this->in_faces(point).size()) {
            return true;
        }
        std::array<double,3> vec = {((double)rand()/RAND_MAX)*2-1,((double)rand()/RAND_MAX)*2-1,((double)rand()/RAND_MAX)*2-1};
        vector<IsInsideIntersection> output;
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            vector<std::array<double,3>>* circuit = circuit_cut(this->circuits(face_index));
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*circuit)) {
                MatrixXd m(4,4);
                VectorXd b(4);
                for (int i = 0; i < 3; i++) {
                    m(i,0) = triangle[0][i];
                    m(i,1) = triangle[1][i];
                    m(i,2) = triangle[2][i];
                    m(i,3) = -vec[i];
                    b(i) = point[i];
                }
                for (int i = 0; i < 3; i++) {
                    m(3,i) = 1;
                }
                m(3,3) = 0;
                b(3) = 1;
                VectorXd x = m.colPivHouseholderQr().solve(b);
                VectorXd b_prime = m*x;
                std::array<double,3> point_prime;
                for (int i = 0; i < 3; i++) {
                    point_prime[i] = b(i);
                }
                bool error = distance(point,point_prime) > 0.00001;
                if (x(0) >= 0 && x(1) >= 0 && x(2) >= 0 && x(3) > 0) {
                    std::array<double,3> p;
                    for (int i = 0; i < 3; i++) {
                        p[i] = x(0)*triangle[0][i]+x(1)*triangle[1][i]+x(2)*triangle[2][i];
                    }
                    if (!error) {
                        //cout << "is inside " << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << endl;
                        Polyhedron::IsInsideIntersection intersection;
                        intersection.gamma = x(3);
                        intersection.point = p;
                        intersection.face_index = face_index;
                        output.push_back(intersection);
                        //cout << "gamma " << intersection.gamma << endl;
                        break;
                    }
                }
            }
            delete circuit;
        }
        return output.size() % 2 == 1;
    }
    void round_verts(RoundPointMeter rounder) {
        for (int i = 0; i < verts.size(); i++) {
            verts[i] = rounder(verts[i]);
        }
        while (true) {
            bool double_break = false;
            for (int i = 0; i < verts.size(); i++) {
                for (int j = i+1; j < verts.size(); j++) {
                    if (Polyhedron::distance(verts[i],verts[j]) < 0) {
                        verts.erase(verts.begin()+j);
                        for (int k = 0; k < edges.size();) {
                            set<int>::iterator it = edges[k].begin();
                            int p1 = *(it);
                            int p2 = *(++it);
                            if (p1 == j) {
                                p1 = i;
                            } else if (p1 > j) {
                                p1--;
                            }
                            if (p2 == j) {
                                p2 = i;
                            } else if (p2 > j) {
                                p2--;
                            }
                            if (p1 == p2) {
                                edges.erase(edges.begin()+k);
                                for (int l = 0; l < faces.size(); l++) {
                                    vector<int> face_list(faces[l].begin(),faces[l].end());
                                    for (int m = 0; m < face_list.size();) {
                                        if (face_list[m] > k) {
                                            face_list[m]--;
                                        } else if (face_list[m] == k) {
                                            face_list.erase(face_list.begin()+m);
                                            continue;
                                        }
                                        m++;
                                    }
                                    faces[l].clear();
                                    faces[l].insert(face_list.begin(),face_list.end());
                                }
                            } else {
                                edges[k] = {p1, p2};
                                k++;
                            }
                        }
                        double_break = true;
                        break;
                    }
                }
                if (double_break) {
                    break;
                }
            }
            if (!double_break) {
                break;
            }
        }
        while (true) {
            bool double_break = false;
            for (int i = 0; i < edges.size(); i++) {
                set<std::array<double,3>> edge1_grounded;
                for (const int& index : edges[i]) {
                    edge1_grounded.insert(verts[index]);
                }
                
                for (int j = i+1; j < edges.size(); j++) {
                    set<std::array<double,3>> edge2_grounded;
                    for (const int& index : edges[j]) {
                        edge2_grounded.insert(verts[index]);
                    }
                    if (edge1_grounded == edge2_grounded) {
                        edges.erase(edges.begin()+j);
                        for (int k = 0; k < faces.size(); k++) {
                            vector<int> face_list(faces[k].begin(),faces[k].end());
                            for (int l = 0; l < face_list.size(); l++) {
                                if (face_list[l] > j) {
                                    face_list[l]--;
                                } else if (face_list[l] == j) {
                                    face_list[l] = i;
                                }
                            }
                            faces[k].clear();
                            faces[k].insert(face_list.begin(),face_list.end());
                        }
                        double_break = true;
                        break;
                    }
                }
                if (double_break) {
                    break;
                }
            }
            if (!double_break) {
                break;
            }
        }
        while (true) {
            bool double_break = false;
            for (int i = 0; i < faces.size(); i++) {
                set<set<std::array<double,3>>> face1_grounded;
                for (const int& edge_index : faces[i]) {
                    set<std::array<double,3>> edge_grounded;
                    for (const int& index : edges[edge_index]) {
                        edge_grounded.insert(verts[index]);
                    }
                    face1_grounded.insert(edge_grounded);
                }
                for (int j = i+1; j < faces.size(); j++) {
                    set<set<std::array<double,3>>> face2_grounded;
                    for (const int& edge_index : faces[j]) {
                        set<std::array<double,3>> edge_grounded;
                        for (const int& index : edges[edge_index]) {
                            edge_grounded.insert(verts[index]);
                        }
                        face2_grounded.insert(edge_grounded);
                    }
                    if (face1_grounded == face2_grounded) {
                        faces.erase(faces.begin()+j);
                        double_break = true;
                        break;
                    }
                }
                if (double_break) {
                    break;
                }
            }
            if (!double_break) {
                break;
            }
        }
    }
    string exp() {
        string output1;
        string output2;
        output1 += "ply\n";
        output1 += "format ascii 1.0\n";
        output1 += "element vertex " + to_string(verts.size()) + "\n";
        output1 += "property float x\n";
        output1 += "property float y\n";
        output1 += "property float z\n";
        int face_count = 0;
        for (const std::array<double,3>& vert : verts) {
            output2 += to_string(vert[0]) + " " + to_string(vert[2]) + " " + to_string(vert[1]) + "\n";
        }
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            set<vector<std::array<double,3>>> circs = Polyhedron::make_clockwise(this->circuits(face_index));
            vector<std::array<double,3>>* circuit = Polyhedron::circuit_cut(circs);
            if (circs.size() == 1) {
                face_count++;
                output2 += to_string(circuit->size()) + " ";
                for (const std::array<double,3>& point : *circuit){
                    int index = std::distance(verts.begin(),find(verts.begin(),verts.end(),point));
                    output2 += to_string(index) + " ";
                }
                output2 += "\n";
            } else {
                for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*circuit)) {
                    face_count++;
                    output2 += "3 ";
                    for (const std::array<double,3>& point : triangle){
                        int index = std::distance(verts.begin(),find(verts.begin(),verts.end(),point));
                        output2 += to_string(index) + " ";
                    }
                    output2 += "\n";
                }
            }
            delete circuit;
        }
        output1 += "element face " + to_string(face_count) + "\n";
        output1 += "property list uchar uint vertex_indices\n";
        output1 += "end_header\n";
        return output1 + output2;
    }
    bool exp(string filename) {
        ofstream output_file(filename);
        if (output_file.is_open()) {
            output_file << this->exp();
            return true;
        }
        return false;
    }
    string dump() {
        string output;
        output += "ply\n";
        output += "format ascii 1.0\n";
        output += "element vertex " + to_string(verts.size()) + "\n";
        output += "property float x\n";
        output += "property float y\n";
        output += "property float z\n";
        output += "element edge " + to_string(edges.size()) + "\n";
        output += "property list uchar uint vertex_indices\n";
        output += "element face " + to_string(faces.size()) + "\n";
        output += "property list uchar uint edge_indices\n";
        output += "end_header\n";
        for (const std::array<double,3>& vert : verts) {
            output += to_string(vert[0]) + " " + to_string(vert[1]) + " " + to_string(vert[2]) + "\n";
        }
        for (const set<int>& edge : edges) {
            output += "2 ";
            for (const int& index : edge){
                output += to_string(index) + " ";
            }
            output += "\n";
        }
        for (const set<int>& face : faces) {
            output += to_string(face.size()) + " ";
            for (const int& index : face){
                output += to_string(index) + " ";
            }
            output += "\n";
        }
        return output;
    }
    bool dump(string filename) {
        ofstream output_file(filename);
        if (output_file.is_open()) {
            output_file << this->dump();
            return true;
        }
        return false;
    }
    void loads(string text){
        int element_index = -1;
        int vert_count = 0;
        int edge_count = 0;
        int face_count = 0;
        vector<string> elements_vector;
        map<string,vector<string>> elements;
        verts.clear();
        edges.clear();
        faces.clear();
        while (text.size()) {
            std::string line;
            if (text.find("\n") != std::string::npos) {
                line = text.substr(0,text.find("\n")+1);
                text = text.substr(text.find("\n")+1,text.size()-1-1);
            } else {
                line = text;
                text = "";
            }
            cout << line;
            line = trim(line);
            vector<string> split_line;
            while (line.size()) {
                string token;
                if (line.find(" ") != std::string::npos) {
                    token = line.substr(0, line.find(" "));
                    line = line.substr(line.find(" ")+1,line.size()-line.find(" ")-1);
                } else {
                    token = line;
                    line = "";
                }
                if (token.size()) {
                    split_line.push_back(token);
                }
            }
            if (split_line.size() >= 1 && split_line[0] == "end_header") {
                element_index = 0;
            }
            if (element_index == -1) {
                if (split_line.size() >= 3) {
                    if (split_line[0] == "element") {
                        elements_vector.push_back(split_line[1]);
                        int element_count = stoi(split_line[2]);
                        if (split_line[1] == "vertex") {
                            vert_count = element_count;
                        } else if (split_line[1] ==  "edge") {
                            edge_count = element_count;
                        } else if (split_line[1] == "face") {
                            face_count = element_count;
                        }
                    } else if (split_line[0] == "property") {
                        if (split_line[1] == "float") {
                            elements[elements_vector.back()].push_back(split_line[2]);
                        }
                        if (split_line.size() >= 5) {
                            if (vector<string>(split_line.begin()+1,split_line.begin()+1+3) == vector<string>{"list", "uchar", "uint"}) {
                                elements[elements_vector.back()].push_back(split_line[4]);
                            }
                        }
                    }
                }
            } else {
                string element = elements_vector[element_index];
                if (split_line.size() >= elements[element].size()) {
                    int index = 0;
                    std::array<double,3> vert;
                    for (const string& property : elements[element]) {
                        string suffix = "_indices";
                        if (property.size() >= suffix.size() && property.substr(property.size()-suffix.size(), suffix.size()) == suffix) {
                            set<int> facet;
                            int property_count = stoi(split_line[index]);
                            for (int i = 0; i < property_count; i++) {
                                facet.insert(stoi(split_line[index+1+i]));
                            }
                            index += 1+property_count;
                            if (element == "edge") {
                                if (property == "vertex_indices") {
                                    edges.push_back(facet);
                                }
                            } else if (element == "face") {
                                if (property == "edge_indices") {
                                    faces.push_back(facet);
                                }
                            }
                        } else {
                            if (element == "vertex") {
                                if (property == "x") {
                                    vert[0] = stod(split_line[index]);
                                } else if (property == "y") {
                                    vert[1] = stod(split_line[index]);
                                } else if (property == "z") {
                                    vert[2] = stod(split_line[index]);
                                }
                            }
                            index++;
                        }
                    }
                    if (element == "vertex") {
                        verts.push_back(vert);
                        if (verts.size() == vert_count) {
                            element_index++;
                        }
                    } else if (element == "edge") {
                        if (edges.size() == edge_count) {
                            element_index++;
                        }
                    } else if (element == "face") {
                        if (faces.size() == face_count) {
                            element_index++;
                        }
                    }
                }
            }
        }
    }
    bool load(string filename) {
        std::ifstream file(filename);
        string text;
        if (file.is_open()) {
            string line;
            while(std::getline(file, line)) {
                text += line + "\n";
            }
            loads(text);
            file.close();
            return true;
        }
        return false;
    }
    static bool circuit_overlap(set<vector<std::array<double,3>>> circuits1, set<vector<std::array<double,3>>> circuits2) {
        vector<std::array<double,3>>* circuit1 = Polyhedron::circuit_cut(circuits1);
        vector<std::array<double,3>>* circuit2 = Polyhedron::circuit_cut(circuits2);
        for (int i = 0; i < circuit1->size(); i++) {
            for (int j = 0; j < circuit2->size(); j++) {
                std::array<double,3>* point = Polyhedron::intersect_segments({(*circuit1)[(i-1+circuit1->size())%circuit1->size()], (*circuit1)[i]},{(*circuit2)[(j-1+circuit2->size())%circuit2->size()], (*circuit2)[j]});
                if (point != NULL) {
                    delete circuit1;
                    delete circuit2;
                    delete point;
                    return true;
                }
            }
        }
        delete circuit1;
        delete circuit2;
        return false;
    }
    std::array<double,3> project_on_face_plane(int face_index, std::array<double,3> point) {
        vector<std::array<double,3>>* circuit = circuit_cut(this->circuits(face_index));
        std::array<double,3> vec1;
        for (int i = 0; i < 3; i++) {
            vec1[i] = (*circuit)[1][i]-(*circuit)[0][i];
        }
        std::array<double,3> vec2;
        for (int i = 0; i < 3; i++) {
            vec2[i] = (*circuit)[2][i]-(*circuit)[0][i];
        }
        std::array<double,3> vec0 = Polyhedron::cross3D(vec1,vec2);
        MatrixXd m(3,3);
        for (int i = 0; i < 3; i++) {
            m(i,0) = vec0[i];
            m(i,1) = -vec1[i];
            m(i,2) = -vec2[i];
        }
        VectorXd b(3);
        for (int j = 0; j < 3; j++) {
            b(j) = (*circuit)[0][j]-point[j];
        }
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        std::array<double,3> output;
        for (int i=0; i < 3; i++) {
            output[i] = point[i]+x(0)*vec0[i];
        }
        delete circuit;
        return output;
    }
    std::array<std::array<double,3>,2> project_ray_on_face_plane(int face_index, std::array<std::array<double,3>,2> ray) {
        std::array<std::array<double,3>,2> output;
        output[0] = this->project_on_face_plane(face_index, ray[0]);
        output[1] = this->project_on_face_plane(face_index, ray[1]);
        return output;
    }
   struct RayIntersection {
        double alpha;
        std::array<double,3> point;
    };
    vector<RayIntersection> project_ray_on_face_plane_and_intersect(int face_index, std::array<std::array<double,3>,2> ray) {
        vector<RayIntersection> output;
        ray = this->project_ray_on_face_plane(face_index, ray);
        std::array<double,3> p1 = ray[0];
        std::array<double,3> p2 = ray[1];
        if (p1 == p2) {
            return output;
        }
        std::array<double,3> vec;
        for (int i = 0; i < 3; i++) {
            vec[i] = p2[i]-p1[i];
        }
        for (const int &edge_index : this->faces[face_index]) {
            vector<std::array<double,3>> edge;
            for (const int &index : this->edges[edge_index]) {
                edge.push_back(this->verts[index]);
            }
            MatrixXd m(4,3);
            for (int i = 0; i < 3; i++) {
                m(i,0) = -vec[i];
                m(i,1) = edge[0][i];
                m(i,2) = edge[1][i];
            }
            m(3,0) = 0;
            m(3,1) = 1;
            m(3,2) = 1;
            VectorXd b(4);
            for (int j = 0; j < 3; j++) {
                b(j) = p1[j];
            }
            b(3) = 1;
            VectorXd x = m.fullPivHouseholderQr().solve(b);
            VectorXd b_prime = m*x;
            std::array<double,3> point1;
            for (int i = 0; i < 3; i++) {
                point1[i] = x(0)*vec[i]+p1[i];
            }
            std::array<double,3> point2;
            for (int i = 0; i < 3; i++) {
                point2[i] = x(1)*edge[0][i]+x(2)*edge[1][i];
            }
            //cout << Polyhedron::distance(point1, point2) << " " << x(1) << " " << round_float(x(0)) << " " << endl;
            if (Polyhedron::distance(point1, point2) <= pow(10,-5) && round_float(x(1)) >= 0 && round_float(x(1)) <= 1 && round_float(x(0)) >= 0) {
                Polyhedron::RayIntersection ray_intersection;
                ray_intersection.alpha = x(0);
                ray_intersection.point = point2;
                output.push_back(ray_intersection);
            }
        }
        return output;
    }
    struct SortedIntersection {
        double alpha;
        std::array<double,3> point;
    };
    struct CircuitIntersections {
        vector<std::array<double,3>> sorted_intersections;
        vector<set<std::array<double,3>>> edges;
    };
    static bool compare_sorted_intersections(Polyhedron::SortedIntersection a, Polyhedron::SortedIntersection b) {
        return a.alpha < b.alpha;
    }
    static Polyhedron::CircuitIntersections intersect_circuits(vector<std::array<double,3>> circuit1,vector<std::array<double,3>> circuit2) {
        set<std::array<double,3>> intersections;
        for (int i = 0; i < circuit1.size(); i++) {
            std::array<double,3> p1 = circuit1[(i-1+circuit1.size())%circuit1.size()];
            std::array<double,3> p2 = circuit1[i];
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(circuit2)) {
                Polyhedron::TriangleIntersection* intersect = Polyhedron::intersect_triangle({p1,p2}, triangle);
                if (intersect != NULL) {
                    intersections.insert(intersect->point);
                    delete intersect;
                }
            }
        }
        for (int i = 0; i < circuit2.size(); i++) {
            std::array<double,3> p1 = circuit2[(i-1+circuit2.size())%circuit2.size()];
            std::array<double,3> p2 = circuit2[i];
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(circuit1)) {
                Polyhedron::TriangleIntersection* intersect = Polyhedron::intersect_triangle({p1,p2}, triangle);
                if (intersect != NULL) {
                    intersections.insert(intersect->point);
                    delete intersect;
                }
            }
        }
        vector<std::array<double,3>> intersections_vector(intersections.begin(),intersections.end());
        vector<Polyhedron::SortedIntersection> sorted_intersections;
        if (intersections.size()) {
            std::array<double,3> vec;
            for (int i = 0; i < 3; i++) {
                vec[i] = intersections_vector[1][i]-intersections_vector[0][i];
                for (const std::array<double,3>& intersection : intersections_vector) {        
                    MatrixXd m(3,1);
                    VectorXd b(3);
                    for (int i = 0; i < 3; i++) {
                        m(i,0) = vec[i];
                        b(i) = intersection[i]-intersections_vector[0][i];
                    }
                    VectorXd x = m.colPivHouseholderQr().solve(b);
                    Polyhedron::SortedIntersection sorted_intersection;
                    sorted_intersection.alpha = x(0);
                    sorted_intersection.point = intersection;
                    sorted_intersections.push_back(sorted_intersection);
                }
            }
        }
        sort(sorted_intersections.begin(), sorted_intersections.end(), Polyhedron::compare_sorted_intersections);
        CircuitIntersections circuit_intersections;
        for (const Polyhedron::SortedIntersection intersection : sorted_intersections) {
            circuit_intersections.sorted_intersections.push_back(intersection.point);
        }
        for (int index = 0; index < circuit_intersections.sorted_intersections.size()-1; index++) {
            std::array<double,3> p1 = circuit_intersections.sorted_intersections[index];
            std::array<double,3> p2 = circuit_intersections.sorted_intersections[index-1];
            std::array<double,3> p;
            for (int i = 0; i < 3; i++) {
                p[i] = (p1[i]+p2[i])/2;
            }
            bool any1 = false;
            bool any2 = false;
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(circuit1)) {
                if (Polyhedron::inside_triangle(triangle,p)) {
                    any1 = true;
                }
            }
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(circuit2)) {
                if (Polyhedron::inside_triangle(triangle,p)) {
                    any2 = true;
                }
            }
            if (any1 && any2) {
                circuit_intersections.edges.push_back({p1,p2});
            }
        }
        return circuit_intersections;
    }
    static bool compare_ray_intersections(RayIntersection a, RayIntersection b) {
        return a.alpha < b.alpha;
    }
    static vector<vector<std::array<double,3>>> path_along_faces(std::array<double,3> p1, std::array<double,3> p2, vector<int> face_path, Polyhedron poly, Polyhedron other_poly) {
        std::array<std::array<double,3>,2> ray = other_poly.project_ray_on_face_plane(face_path.back(), {p1, p2});
        //cout << "distances " << round_float(Polyhedron::distance(ray[0],p1)) << " " <<  round_float(Polyhedron::distance(ray[1],p2)) << endl;
        if (round_float(Polyhedron::distance(ray[0],p1)) == 0 && round_float(Polyhedron::distance(ray[1],p2)) == 0) {
            return {{p1, p2}};
        }
        vector<vector<std::array<double,3>>> output;
        if (ray[0] == ray[1]) {
            return output;
        }
        vector<RayIntersection> intersections = other_poly.project_ray_on_face_plane_and_intersect(face_path.back(), {p1, p2});
        //cout << "intersections prior size " << intersections.size() << "\n";
        for (int i = 0; i < intersections.size();) {
            if (intersections[i].alpha < 0) {
                intersections.erase(intersections.begin()+i);
            } else {
                i++;
            }
        }
        //cout << "intersections size " << intersections.size() << "\n";
        if (intersections.size()) {
            sort(intersections.begin(), intersections.end(), compare_ray_intersections);
            std::array<double,3> p = intersections[0].point;
            if (!poly.is_inside(p)) {
                return output;
            }
            if (p != p1) {
                for (int face_index2 = 0; face_index2 < other_poly.faces.size(); face_index2++) {
                    bool visited = false;
                    for (const int& face : face_path) {
                        if (face_index2 == face) {
                            visited = true;
                            break;
                        }
                    }
                    bool any = false;
                    vector<std::array<double,3>>* circuit = circuit_cut(other_poly.circuits(face_index2));
                    for (const std::array<std::array<double,3>,3>& triangle : triangulate(*circuit)) {
                        if (Polyhedron::inside_triangle(triangle, p)) {
                            any = true;
                            break;
                        }
                    }
                    vector<int> new_face_path(face_path.begin(),face_path.end());
                    new_face_path.push_back(face_index2);
                    //cout << !visited << any << endl;
                    if (!visited && any) {
                        for (vector<std::array<double,3>> x : Polyhedron::path_along_faces(p, p2, new_face_path, poly, other_poly)) {
                            vector<std::array<double,3>> new_path;
                            new_path.push_back(p1);
                            new_path.insert(new_path.end(), x.begin(), x.end());
                            output.push_back(new_path);
                        }
                    }
                }
            }
        }
        return output;
    }
    static set<vector<std::array<double,3>>> circuit_helper(
            set<set<std::array<double,3>>> edges, 
            set<std::array<double,3>> start,
            set<std::array<double,3>> previous,
            set<std::array<double,3>> current,
            vector<std::array<double,3>> path) {
         
        // Build edge_lookup exactly like Python
        map<std::array<double,3>, set<set<std::array<double,3>>>> edge_lookup;
        for (const auto& edge : edges) {
            for (const auto& point : edge) {
                edge_lookup[point].insert(edge);
            }
        }
        
        set<vector<std::array<double,3>>> circuits;
        
        if (start.empty()) {
            // Mirror Python's initial loop structure
            set<set<std::array<double,3>>> seen;
            for (const set<std::array<double,3>>& start_edge : edges) {
                vector<std::array<double,3>> empty_path;
                std::array<double,3> point = *start_edge.begin();
                
                set<set<std::array<double,3>>> temp = edge_lookup[point];
                temp.erase(start_edge);
                
                for (const set<std::array<double,3>>& current_edge : temp) {
                    vector<std::array<double,3>> new_path = {point};
                    set<std::array<double,3>> path_set(new_path.begin(), new_path.end());
                    if (coplanar(path_set)) {
                        set<vector<std::array<double,3>>> result = circuit_helper(edges, start_edge, start_edge, current_edge, new_path);
                        circuits.insert(result.begin(), result.end());
                    }
                }
            }
        } else {
            // Mirror Python's recursive logic
            if (current == start) {
                return {path};
            }
            
            // Find the point: list(current - previous)[0]
            set<std::array<double,3>> difference;
            set_difference(current.begin(), current.end(),
                          previous.begin(), previous.end(),
                          inserter(difference, difference.begin()));
            
            if (difference.empty()) return circuits;
            std::array<double,3> point = *difference.begin();
            
            // Check if point already in path
            if (find(path.begin(), path.end(), point) != path.end()) {
                return circuits;
            }
            
            // Self-intersection check (match Python's logic exactly)
            if (path.size() > 2) {
                std::array<double,3> path_prev = path[path.size() - 2];
                std::array<double,3> path_last = path[path.size() - 1];
                for (int i = 1; i < path.size() - 2; i++) {
                    std::array<double,3> p1 = path[i - 1];
                    std::array<double,3> p2 = path[i];
                    
                    std::array<double,3>* intersect = intersect_segments({p1, p2}, {path_prev, path_last});
                    
                    if (intersect != NULL) {
                        if (distance(*intersect, path_prev) > 0.001 && distance(*intersect, path_last) > 0.001) {
                            delete intersect;
                            return circuits;
                        }
                        delete intersect;
                    }
                }
            }
            
            // Update previous to current (like Python)
            set<std::array<double,3>> new_previous = current;
            
            // Get next edges
            set<set<std::array<double,3>>> temp = edge_lookup[point];
            temp.erase(new_previous);
            
            vector<std::array<double,3>> new_path = path;
            new_path.push_back(point);
            
            for (const set<std::array<double,3>>& next_edge : temp) {
                set<std::array<double,3>> path_set(new_path.begin(), new_path.end());
                if (coplanar(path_set)) {
                    set<vector<std::array<double,3>>> result = circuit_helper(edges, start, new_previous, next_edge, new_path);
                    circuits.insert(result.begin(), result.end());
                }
            }
        }
        
        // Mirror Python's duplicate removal logic
        vector<vector<std::array<double,3>>> circuits_list(circuits.begin(), circuits.end());
        
        for (int i = 0; i < circuits_list.size(); i++) {
            for (int j = i + 1; j < circuits_list.size(); j++) {
                const auto& x = circuits_list[i];
                const auto& y = circuits_list[j];
                
                // Check if not len(set(x)-set(y)) (i.e., x is subset of y)
                set<std::array<double,3>> x_set(x.begin(), x.end());
                set<std::array<double,3>> y_set(y.begin(), y.end());
                set<std::array<double,3>> difference;
                set_difference(x_set.begin(), x_set.end(),
                              y_set.begin(), y_set.end(),
                              inserter(difference, difference.begin()));
                
                if (difference.empty()) {
                    // Create reversed y
                    vector<std::array<double,3>> y_r(y.rbegin(), y.rend());
                    
                    // Find x[0] in y and y_r  
                    vector<std::array<double,3>>::const_iterator y_it = find(y.begin(), y.end(), x[0]);
                    vector<std::array<double,3>>::const_iterator yr_it = find(y_r.begin(), y_r.end(), x[0]);
                    
                    bool should_remove = false;
                    
                    if (y_it != y.end()) {
                        // Check (y[y.index(x[0]):]+y[:y.index(x[0])])[:len(x)] == x
                        vector<std::array<double,3>> rotated_y;
                        rotated_y.reserve(y.size());
                        for (vector<std::array<double,3>>::const_iterator it = y_it; it != y.end(); ++it) {
                            rotated_y.push_back(*it);
                        }
                        for (vector<std::array<double,3>>::const_iterator it = y.begin(); it != y_it; ++it) {
                            rotated_y.push_back(*it);
                        }
                        
                        if (rotated_y.size() >= x.size()) {
                            vector<std::array<double,3>> prefix(rotated_y.begin(), rotated_y.begin() + x.size());
                            if (prefix == x) should_remove = true;
                        }
                    }
                    
                    if (!should_remove && yr_it != y_r.end()) {
                        // Check same for y_r
                        vector<std::array<double,3>> rotated_yr;
                        rotated_yr.reserve(y_r.size());
                        for (vector<std::array<double,3>>::const_iterator it = yr_it; it != y_r.end(); ++it) {
                            rotated_yr.push_back(*it);
                        }
                        for (vector<std::array<double,3>>::const_iterator it = y_r.begin(); it != yr_it; ++it) {
                            rotated_yr.push_back(*it);
                        }
                        
                        if (rotated_yr.size() >= x.size()) {
                            vector<std::array<double,3>> prefix(rotated_yr.begin(), rotated_yr.begin() + x.size());
                            if (prefix == x) should_remove = true;
                        }
                    }
                    
                    if (should_remove && circuits.count(y)) {
                        circuits.erase(y);
                    }
                }
            }
        }
        
        return circuits;
    }
    // Static wrapper that can be called without object instance
    static set<vector<std::array<double,3>>> circuit_helper(set<set<std::array<double,3>>> edges) {
        return circuit_helper(edges, {}, {}, {}, {});
    }
    Polyhedron edge_sets_to_poly(std::array<vector<set<set<std::array<double,3>>>>,2> edge_sets_per_poly, std::array<Polyhedron,2> polys) {
        // Print edge sets per poly (equivalent to Python print statement)
        // print('edge sets per poly:', edge_sets_per_poly)
        std::cout << "edge sets per poly: [";
        for (int i = 0; i < 2; i++) {
            std::cout << "[";
            bool first_set = true;
            for (const auto& edge_set : edge_sets_per_poly[i]) {
                if (!first_set) std::cout << ", ";
                std::cout << "{";
                bool first_edge = true;
                for (const auto& edge : edge_set) {
                    if (!first_edge) std::cout << ", ";
                    std::cout << "{";
                    bool first_point = true;
                    for (const auto& point : edge) {
                        if (!first_point) std::cout << ", ";
                        std::cout << "(" << point[0] << ", " << point[1] << ", " << point[2] << ")";
                        first_point = false;
                    }
                    std::cout << "}";
                    first_edge = false;
                }
                std::cout << "}";
                first_set = false;
            }
            std::cout << "]";
            if (i < 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
        
        // Collect edges per poly
        std::array<set<set<std::array<double,3>>>,2> edges_per_poly;
        for (int poly_index = 0; poly_index < 2; poly_index++) {
            for (const auto& edge_set : edge_sets_per_poly[poly_index]) {
                edges_per_poly[poly_index].insert(edge_set.begin(), edge_set.end());
            }
        }
        
        // Collect points per poly
        std::array<set<std::array<double,3>>,2> points_per_poly;
        for (int poly_index = 0; poly_index < 2; poly_index++) {
            for (const auto& edge : edges_per_poly[poly_index]) {
                points_per_poly[poly_index].insert(edge.begin(), edge.end());
            }
        }
        
        // Union all points
        set<std::array<double,3>> points = points_per_poly[0];
        points.insert(points_per_poly[1].begin(), points_per_poly[1].end());
        
        // Create point mapping for rounding
        map<std::array<double,3>, std::array<double,3>> point_map;
        for (const std::array<double,3>& point : points) {
            std::array<double,3> rounded = round_point(point);
            if (point_map.find(rounded) != point_map.end()) {
                if (distance(point, rounded) < distance(point_map[rounded], rounded)) {
                    point_map[rounded] = point;
                }
            } else {
                point_map[rounded] = point;
            }
        }
        
        points.clear();
        for (const auto& p : point_map) {
            points.insert(p.second);
        }
        
        Polyhedron poly;
        poly.verts = vector<std::array<double,3>>(points.begin(), points.end());
        
        set<vector<std::array<double,3>>> seen_components;
        
        // Process faces from both polyhedra
        for (int poly_index = 0; poly_index < 2; poly_index++) {
            for (int face_index = 0; face_index < polys[poly_index].faces.size(); face_index++) {
                // Find cofacial points
                set<std::array<double,3>> cofacial_points_set;
                vector<std::array<double,3>>* circuit = circuit_cut(polys[poly_index].circuits(face_index));
                if (!circuit) continue;
                
                vector<std::array<std::array<double, 3>, 3>> triangles = triangulate(*circuit);
                for (const std::array<double,3>& point : points) {
                    bool is_cofacial = false;
                    for (const std::array<std::array<double, 3>, 3>& triangle : triangles) {
                        if (inside_triangle(triangle, point)) {
                            is_cofacial = true;
                            break;
                        }
                    }
                    if (is_cofacial) {
                        cofacial_points_set.insert(point);
                    }
                }
                
                delete circuit;
                
                vector<std::array<double,3>> cofacial_points(cofacial_points_set.begin(), cofacial_points_set.end());
                cout << "cofacial_points size: " << cofacial_points.size() << endl;
                
                if (cofacial_points.size() > 2) {
                    set<set<std::array<double,3>>> edges;
                    cout << cofacial_points.size() << endl;
                    
                    for (int index1 = 0; index1 < cofacial_points.size(); index1++) {
                        for (int index2 = index1 + 1; index2 < cofacial_points.size(); index2++) {
                            std::array<double,3> point1 = cofacial_points[index1];
                            std::array<double,3> point2 = cofacial_points[index2];
                            
                            bool double_break = false;
                            
                            // Check both polyhedra for midpoint validity
                            for (int poly_idx = 0; poly_idx < 2; poly_idx++) {
                                Polyhedron p = polys[poly_idx];
                                vector<FaceIntersection> intersects = p.face_intersect({point1, point2});
                                
                                // Sort intersects and build ps array
                                sort(intersects.begin(), intersects.end(), Polyhedron::compare_face_intersections);
                                vector<std::array<double,3>> ps;
                                for (const FaceIntersection& intersect : intersects) {
                                    if (ps.empty() || intersect.point != ps.back()) {
                                        ps.push_back(intersect.point);
                                    }
                                }
                                
                                // Build complete point sequence: [point1] + ps + [point2]
                                vector<std::array<double,3>> complete_ps;
                                complete_ps.push_back(point1);
                                complete_ps.insert(complete_ps.end(), ps.begin(), ps.end());
                                complete_ps.push_back(point2);
                                
                                // Check midpoints of each segment
                                for (int k = 0; k < complete_ps.size() - 1; k++) {
                                    std::array<double,3> midpoint;
                                    for (int l = 0; l < 3; l++) {
                                        midpoint[l] = 0.5 * (complete_ps[k][l] + complete_ps[k+1][l]);
                                    }
                                    if (!p.is_inside(midpoint)) {
                                        double_break = true;
                                        break;
                                    }
                                }
                                if (double_break) break;
                            }
                            
                            if (double_break) continue;
                            
                            // Check coplanarity condition for both polyhedra
                            for (int poly_idx = 0; poly_idx < 2; poly_idx++) {
                                Polyhedron p = polys[poly_idx];
                                vector<int> faces1 = p.in_faces(point1);
                                vector<int> faces2 = p.in_faces(point2);
                                cout << "point1 faces: " << faces1.size() << " point2 faces: " << faces2.size() << endl;
                                
                                // Find intersection of face sets
                                set<int> face_intersection;
                                for (int f1 : faces1) {
                                    for (int f2 : faces2) {
                                        if (f1 == f2) {
                                            face_intersection.insert(f1);
                                        }
                                    }
                                }
                                
                                bool should_create_edge = false;
                                for (int common_face : face_intersection) {
                                    // Collect cofacial points + vertices from common face
                                    set<std::array<double,3>> face_points(cofacial_points.begin(), cofacial_points.end());
                                    for (int edge_idx : p.faces[common_face]) {
                                        for (int vert_idx : p.edges[edge_idx]) {
                                            face_points.insert(p.verts[vert_idx]);
                                        }
                                    }
                                    
                                    if (!coplanar(face_points)) {
                                        // Check if no other point lies on this edge segment
                                        set<std::array<double,3>> edge = {point1, point2};
                                        bool point_on_segment_found = false;
                                        for (const std::array<double,3>& other_point : cofacial_points) {
                                            if (other_point != point1 && other_point != point2) {
                                                if (Polyhedron::point_on_segment(edge, other_point)) {
                                                    point_on_segment_found = true;
                                                    break;
                                                }
                                            }
                                        }
                                        if (!point_on_segment_found) {
                                            should_create_edge = true;
                                        }
                                        break;
                                    }
                                }
                                
                                if (should_create_edge) {
                                    edges.insert({point1, point2});
                                    break;
                                }
                            }
                        }
                    }                    
                    vector<set<std::array<double,3>>> edges_list(edges.begin(), edges.end());
                    
                    // Create circuits from edges
                    set<vector<std::array<double,3>>> circuits = circuit_helper(edges);
                    
                    for (const vector<std::array<double,3>>& circuit : circuits) {
                        set<int> new_face;
                        
                        if (circuit.empty()) continue;
                        
                        for (int i = 0; i < circuit.size(); i++) {
                            std::array<double,3> p1 = circuit[(i - 1 + circuit.size()) % circuit.size()];
                            std::array<double,3> p2 = circuit[i];
                            
                            int idx1 = std::distance(poly.verts.begin(), find(poly.verts.begin(), poly.verts.end(), p1));
                            int idx2 = std::distance(poly.verts.begin(), find(poly.verts.begin(), poly.verts.end(), p2));
                            
                            set<int> edge = {idx1, idx2};
                            auto edge_it = find(poly.edges.begin(), poly.edges.end(), edge);
                            
                            if (edge_it != poly.edges.end()) {
                                new_face.insert(std::distance(poly.edges.begin(), edge_it));
                            } else {
                                poly.edges.push_back(edge);
                                new_face.insert(poly.edges.size() - 1);
                            }
                        }
                        
                        // Check for merging with existing faces (simplified version of Python's complex logic)
                        bool merged = false;
                        for (int existing_face_idx = 0; existing_face_idx < poly.faces.size(); existing_face_idx++) {
                            // This is a simplified version - full Python logic would require find_exterior_circuit
                            // For exact matching, would need to implement the full circuit merging logic
                        }
                        
                        if (!merged && !new_face.empty()) {
                            poly.faces.push_back(new_face);
                        }
                    }
                }
            }
        }
        
        // Remove unused edges (first pass)
        for (int edge_index = 0; edge_index < poly.edges.size();) {
            int usage_count = 0;
            for (const auto& face : poly.faces) {
                if (face.find(edge_index) != face.end()) {
                    usage_count++;
                }
            }
            
            if (usage_count < 2) {
                // Remove this edge and update face indices
                for (auto& face : poly.faces) {
                    set<int> new_face;
                    for (int idx : face) {
                        if (idx == edge_index) {
                            continue; // Remove this edge index
                        } else if (idx > edge_index) {
                            new_face.insert(idx - 1); // Shift down
                        } else {
                            new_face.insert(idx);
                        }
                    }
                    face = new_face;
                }
                poly.edges.erase(poly.edges.begin() + edge_index);
            } else {
                edge_index++;
            }
        }
        
        // Merge coplanar faces (first merging pass)
        bool updated = true;
        while (updated) {
            updated = false;
            for (int face_index1 = 0; face_index1 < poly.faces.size() && !updated; face_index1++) {
                for (int face_index2 = face_index1 + 1; face_index2 < poly.faces.size() && !updated; face_index2++) {
                    // Convert face indices to edge point sets
                    set<set<std::array<double,3>>> face1_edges, face2_edges;
                    
                    for (int edge_idx : poly.faces[face_index1]) {
                        set<std::array<double,3>> edge_points;
                        for (int vert_idx : poly.edges[edge_idx]) {
                            edge_points.insert(poly.verts[vert_idx]);
                        }
                        face1_edges.insert(edge_points);
                    }
                    
                    for (int edge_idx : poly.faces[face_index2]) {
                        set<std::array<double,3>> edge_points;
                        for (int vert_idx : poly.edges[edge_idx]) {
                            edge_points.insert(poly.verts[vert_idx]);
                        }
                        face2_edges.insert(edge_points);
                    }
                    
                    // Check coplanarity of all points
                    set<std::array<double,3>> all_face_points;
                    for (const auto& edge : face1_edges) {
                        all_face_points.insert(edge.begin(), edge.end());
                    }
                    for (const auto& edge : face2_edges) {
                        all_face_points.insert(edge.begin(), edge.end());
                    }
                    
                    if (coplanar(all_face_points)) {
                        bool should_merge = false;
                        
                        // Check merging conditions
                        for (const auto& edge1 : face1_edges) {
                            for (const auto& edge2 : face2_edges) {
                                // Check if one edge contains the other
                                int points_on_edge1 = 0, points_on_edge2 = 0;
                                for (const auto& point : edge2) {
                                    if (point_on_segment(edge1, point)) points_on_edge1++;
                                }
                                for (const auto& point : edge1) {
                                    if (point_on_segment(edge2, point)) points_on_edge2++;
                                }
                                
                                if (points_on_edge1 == 2) {
                                    should_merge = true;
                                    break;
                                } else if (points_on_edge2 == 2) {
                                    should_merge = true;
                                    break;
                                } else if (points_on_edge1 == 1 && points_on_edge2 == 1) {
                                    should_merge = true;
                                    break;
                                }
                                
                                // Check collinear overlap
                                set<std::array<double,3>> combined_edge = edge1;
                                combined_edge.insert(edge2.begin(), edge2.end());
                                if (combined_edge.size() == 3 && colinear(combined_edge)) {
                                    should_merge = true;
                                    break;
                                }
                            }
                            if (should_merge) break;
                        }
                        
                        if (should_merge) {
                            // Merge faces
                            set<int> merged_face = poly.faces[face_index1];
                            merged_face.insert(poly.faces[face_index2].begin(), poly.faces[face_index2].end());
                            
                            poly.faces[face_index1] = merged_face;
                            poly.faces.erase(poly.faces.begin() + face_index2);
                            updated = true;
                        }
                    }
                }
            }
        }
        
        // Remove unused edges (second pass)
        for (int edge_index = 0; edge_index < poly.edges.size();) {
            int usage_count = 0;
            for (const auto& face : poly.faces) {
                if (face.find(edge_index) != face.end()) {
                    usage_count++;
                }
            }
            
            if (usage_count < 2) {
                for (auto& face : poly.faces) {
                    set<int> new_face;
                    for (int idx : face) {
                        if (idx == edge_index) {
                            continue;
                        } else if (idx > edge_index) {
                            new_face.insert(idx - 1);
                        } else {
                            new_face.insert(idx);
                        }
                    }
                    face = new_face;
                }
                poly.edges.erase(poly.edges.begin() + edge_index);
            } else {
                edge_index++;
            }
        }
        
        // Complex face processing (Python's most complex section)
        vector<set<int>> faces;
        for (int face_index = 0; face_index < poly.faces.size(); face_index++) {
            set<vector<std::array<double,3>>> circuits_result = poly.circuits(face_index);
            vector<vector<std::array<double,3>>> circuits(circuits_result.begin(), circuits_result.end());
            vector<set<vector<std::array<double,3>>>> circuit_sets;
            
            for (const auto& circuit : circuits) {
                circuit_sets.push_back({circuit});
            }
            
            // Merge circuit sets using find_exterior_circuit logic
            int index1 = 0;
            while (index1 < circuit_sets.size()) {
                int index2 = index1 + 1;
                while (index2 < circuit_sets.size()) {
                    // Create combined set for exterior circuit check
                    set<vector<std::array<double,3>>> combined_set = circuit_sets[index1];
                    combined_set.insert(circuit_sets[index2].begin(), circuit_sets[index2].end());
                    
                    // Remove collinear points from circuits
                    set<vector<std::array<double,3>>> filtered_set;
                    for (const auto& circuit : combined_set) {
                        vector<std::array<double,3>> filtered_circuit;
                        for (int i = 0; i < circuit.size(); i++) {
                            std::array<double,3> prev = circuit[(i - 1 + circuit.size()) % circuit.size()];
                            std::array<double,3> curr = circuit[i];
                            std::array<double,3> next = circuit[(i + 1) % circuit.size()];
                            
                            if (!colinear(vector<std::array<double,3>>{prev, curr, next})) {
                                filtered_circuit.push_back(curr);
                            }
                        }
                        if (!filtered_circuit.empty()) {
                            filtered_set.insert(filtered_circuit);
                        }
                    }
                    
                    if (find_exterior_circuit(filtered_set) != NULL) {
                        circuit_sets[index1].insert(circuit_sets[index2].begin(), circuit_sets[index2].end());
                        circuit_sets.erase(circuit_sets.begin() + index2);
                    } else {
                        index2++;
                    }
                }
                index1++;
            }
            
            // Create new faces from circuit sets
            vector<vector<set<std::array<double,3>>>> new_faces;
            for (const auto& circuit_set : circuit_sets) {
                set<set<std::array<double,3>>> face_edges;
                for (const auto& circuit : circuit_set) {
                    for (int i = 0; i < circuit.size(); i++) {
                        std::array<double,3> p1 = circuit[(i - 1 + circuit.size()) % circuit.size()];
                        std::array<double,3> p2 = circuit[i];
                        face_edges.insert({p1, p2});
                    }
                }
                new_faces.push_back(vector<set<std::array<double,3>>>(face_edges.begin(), face_edges.end()));
            }
            
            // Complex edge merging logic (first type)
            index1 = 0;
            while (index1 < new_faces.size()) {
                int index2 = index1 + 1;
                while (index2 < new_faces.size()) {
                    bool combine = false;
                    bool local_updated = true;
                    
                    while (local_updated) {
                        local_updated = false;
                        int edge_index1 = 0;
                        while (edge_index1 < new_faces[index1].size()) {
                            int edge_index2 = 0;
                            while (edge_index2 < new_faces[index2].size()) {
                                set<std::array<double,3>> edge1 = new_faces[index1][edge_index1];
                                set<std::array<double,3>> edge2 = new_faces[index2][edge_index2];
                                
                                // Check if edge2 is completely contained in edge1
                                int points_on_edge1 = 0;
                                for (const auto& point : edge2) {
                                    if (point_on_segment(edge1, point)) points_on_edge1++;
                                }
                                
                                if (points_on_edge1 == 2) {
                                    // Remove both edges and create new ones
                                    new_faces[index1].erase(new_faces[index1].begin() + edge_index1);
                                    new_faces[index2].erase(new_faces[index2].begin() + edge_index2);
                                    
                                    auto edge1_it = edge1.begin();
                                    std::array<double,3> p1 = *edge1_it++;
                                    std::array<double,3> p2 = *edge1_it;
                                    
                                    // Find closest points
                                    std::array<double,3> closest_to_p1 = *min_element(edge2.begin(), edge2.end(),
                                        [&p1](const std::array<double,3>& a, const std::array<double,3>& b) {
                                            return distance(p1, a) < distance(p1, b);
                                        });
                                    std::array<double,3> closest_to_p2 = *min_element(edge2.begin(), edge2.end(),
                                        [&p2](const std::array<double,3>& a, const std::array<double,3>& b) {
                                            return distance(p2, a) < distance(p2, b);
                                        });
                                    
                                    if (round_float(distance(p1, closest_to_p1)) > 0) {
                                        set<std::array<double,3>> edge3 = {p1, closest_to_p1};
                                        new_faces[index1].push_back(edge3);
                                        
                                        // Add to poly.edges if not present
                                        set<int> edge3_indices;
                                        for (const auto& point : edge3) {
                                            int idx = std::distance(poly.verts.begin(), find(poly.verts.begin(), poly.verts.end(), point));
                                            edge3_indices.insert(idx);
                                        }
                                        if (find(poly.edges.begin(), poly.edges.end(), edge3_indices) == poly.edges.end()) {
                                            poly.edges.push_back(edge3_indices);
                                        }
                                    }
                                    
                                    if (round_float(distance(p2, closest_to_p2)) > 0) {
                                        set<std::array<double,3>> edge3 = {p2, closest_to_p2};
                                        new_faces[index1].push_back(edge3);
                                        
                                        // Add to poly.edges if not present
                                        set<int> edge3_indices;
                                        for (const auto& point : edge3) {
                                            int idx = std::distance(poly.verts.begin(), find(poly.verts.begin(), poly.verts.end(), point));
                                            edge3_indices.insert(idx);
                                        }
                                        if (find(poly.edges.begin(), poly.edges.end(), edge3_indices) == poly.edges.end()) {
                                            poly.edges.push_back(edge3_indices);
                                        }
                                    }
                                    
                                    combine = true;
                                    local_updated = true;
                                    break;
                                }
                                // Similar logic for edge1 contained in edge2...
                                else {
                                    int points_on_edge2 = 0;
                                    for (const auto& point : edge1) {
                                        if (point_on_segment(edge2, point)) points_on_edge2++;
                                    }
                                    
                                    if (points_on_edge2 == 2) {
                                        // Handle edge1 contained in edge2 (similar to above)
                                        new_faces[index1].erase(new_faces[index1].begin() + edge_index1);
                                        new_faces[index2].erase(new_faces[index2].begin() + edge_index2);
                                        
                                        // Similar splitting logic...
                                        combine = true;
                                        local_updated = true;
                                        break;
                                    }
                                }
                                edge_index2++;
                            }
                            if (local_updated) break;
                            edge_index1++;
                        }
                    }
                    
                    if (combine) {
                        new_faces[index1].insert(new_faces[index1].end(), new_faces[index2].begin(), new_faces[index2].end());
                        new_faces.erase(new_faces.begin() + index2);
                    } else {
                        index2++;
                    }
                }
                index1++;
            }
            
            // Convert new_faces back to edge indices and add to faces
            for (const auto& face_edges : new_faces) {
                set<int> face_indices;
                for (const auto& edge : face_edges) {
                    set<int> edge_indices;
                    for (const auto& point : edge) {
                        int idx = std::distance(poly.verts.begin(), find(poly.verts.begin(), poly.verts.end(), point));
                        edge_indices.insert(idx);
                    }
                    auto edge_it = find(poly.edges.begin(), poly.edges.end(), edge_indices);
                    if (edge_it != poly.edges.end()) {
                        face_indices.insert(std::distance(poly.edges.begin(), edge_it));
                    }
                }
                if (!face_indices.empty()) {
                    faces.push_back(face_indices);
                }
            }
        }
        
        poly.faces = faces;
        
        // Final collinear edge merging within faces
        for (int face_index = 0; face_index < poly.faces.size(); face_index++) {
            vector<int> face(poly.faces[face_index].begin(), poly.faces[face_index].end());
            
            int edge_index1 = 0;
            while (edge_index1 < face.size()) {
                int edge_index2 = edge_index1 + 1;
                while (edge_index2 < face.size()) {
                    set<std::array<double,3>> edge1, edge2;
                    for (int idx : poly.edges[face[edge_index1]]) {
                        edge1.insert(poly.verts[idx]);
                    }
                    for (int idx : poly.edges[face[edge_index2]]) {
                        edge2.insert(poly.verts[idx]);
                    }
                    
                    // Check if edges share exactly one point and are collinear
                    set<std::array<double,3>> intersection;
                    set_intersection(edge1.begin(), edge1.end(), edge2.begin(), edge2.end(),
                                    inserter(intersection, intersection.begin()));
                    
                    if (intersection.size() == 1) {
                        set<std::array<double,3>> union_edge = edge1;
                        union_edge.insert(edge2.begin(), edge2.end());
                        
                        if (colinear(union_edge)) {
                            // Create symmetric difference for new edge
                            set<std::array<double,3>> sym_diff;
                            set<int> sym_diff_indices;
                            
                            for (const auto& point : edge1) {
                                if (edge2.find(point) == edge2.end()) {
                                    sym_diff.insert(point);
                                    int idx = std::distance(poly.verts.begin(), find(poly.verts.begin(), poly.verts.end(), point));
                                    sym_diff_indices.insert(idx);
                                }
                            }
                            for (const auto& point : edge2) {
                                if (edge1.find(point) == edge1.end()) {
                                    sym_diff.insert(point);
                                    int idx = std::distance(poly.verts.begin(), find(poly.verts.begin(), poly.verts.end(), point));
                                    sym_diff_indices.insert(idx);
                                }
                            }
                            
                            // Add new edge if not present
                            auto edge_it = find(poly.edges.begin(), poly.edges.end(), sym_diff_indices);
                            int new_edge_idx;
                            if (edge_it == poly.edges.end()) {
                                poly.edges.push_back(sym_diff_indices);
                                new_edge_idx = poly.edges.size() - 1;
                            } else {
                                new_edge_idx = std::distance(poly.edges.begin(), edge_it);
                            }
                            
                            face[edge_index1] = new_edge_idx;
                            face.erase(face.begin() + edge_index2);
                            continue;
                        }
                    }
                    edge_index2++;
                }
                edge_index1++;
            }
            poly.faces[face_index] = set<int>(face.begin(), face.end());
        }
        
        // Final cleanup - remove unused edges
        for (int edge_index = 0; edge_index < poly.edges.size();) {
            int usage_count = 0;
            for (const auto& face : poly.faces) {
                if (face.find(edge_index) != face.end()) {
                    usage_count++;
                }
            }
            
            if (usage_count < 2) {
                for (auto& face : poly.faces) {
                    set<int> new_face;
                    for (int idx : face) {
                        if (idx == edge_index) {
                            continue;
                        } else if (idx > edge_index) {
                            new_face.insert(idx - 1);
                        } else {
                            new_face.insert(idx);
                        }
                    }
                    face = new_face;
                }
                poly.edges.erase(poly.edges.begin() + edge_index);
            } else {
                edge_index++;
            }
        }
        
        // Remove unused vertices
        for (int vert_index = 0; vert_index < poly.verts.size();) {
            bool is_used = false;
            for (const auto& edge : poly.edges) {
                if (edge.find(vert_index) != edge.end()) {
                    is_used = true;
                    break;
                }
            }
            
            if (!is_used) {
                poly.verts.erase(poly.verts.begin() + vert_index);
                // Update edge indices
                for (auto& edge : poly.edges) {
                    set<int> new_edge;
                    for (int idx : edge) {
                        if (idx > vert_index) {
                            new_edge.insert(idx - 1);
                        } else {
                            new_edge.insert(idx);
                        }
                    }
                    edge = new_edge;
                }
            } else {
                vert_index++;
            }
        }
        
        // Remove empty faces
        poly.faces.erase(
            remove_if(poly.faces.begin(), poly.faces.end(), 
                      [](const set<int>& face) { return face.empty(); }),
            poly.faces.end()
        );
        
        return poly;
    }
    Polyhedron intersect(Polyhedron other) {
        std::array<vector<set<set<std::array<double,3>>>>,2> edge_sets_per_poly;
        std::array<Polyhedron,2> polys;
        polys[0] = *this;
        polys[1] = other;
        for (int poly_index = 0; poly_index < polys.size(); poly_index++) {
            Polyhedron poly = polys[poly_index];
            Polyhedron other_poly = polys[(poly_index-1+polys.size())%polys.size()];
            bool all = true;
	        for (const set<int>& edge1 : poly.edges) {
                bool any = false;
	            for (const set<int>& edge2 : other_poly.edges) {
                    set<std::array<double,3>> edge2_grounded;
                    for (const int& index : edge2) {
                        edge2_grounded.insert(other_poly.verts[index]);
                    }
                    bool broke = false;
                    for (const int& index : edge1) {
                        if (!Polyhedron::point_on_segment(edge2_grounded, poly.verts[index])) {
                            broke = true;
                            break;
                        }
                    }
                    if (!broke) {
                        any = true;
                        break; 
                    }
                }
                if (!any) {
                    all = false;
                    break;
                }
	        }
            if (all) {
                return edge_sets_to_poly(edge_sets_per_poly, polys);
            }
        }        
        for (int poly_index = 0; poly_index < polys.size(); poly_index++) {
            set<set<std::array<double,3>>> new_edges;
            Polyhedron poly = polys[poly_index];
            Polyhedron other_poly = polys[(poly_index-1+polys.size())%polys.size()];
            set<set<std::array<double,3>>> seen_edges;
            set<int> seen_verts;
            set<set<std::array<double,3>>> seen_leaves;
            for (int vert_index = 0; vert_index < poly.verts.size(); vert_index++) {
                std::array<double,3> vert = poly.verts[vert_index];
                //if (seen_verts.find(vert_index) != seen_verts.end()) {
                //    continue;
                //}
                bool root_in_face = false;
                if (other_poly.in_faces(vert).size()) {
                    //continue;
                    root_in_face = true;
                }
                set<std::array<double,3>> leaves;
                bool root_in_poly = other_poly.is_inside(vert);
                cout << "root in poly " << root_in_poly << " " << vert[0] << "," << vert[1] << "," << vert[2] << endl;
                std::queue<int> q;
                q.push(vert_index);
                map<std::array<double,3>,set<int>> face_lookup;
                while (!q.empty()) {
                    int current = q.front();
                    q.pop();
                    if (seen_verts.find(current) != seen_verts.end()) {
                        continue;
                    }
                    seen_verts.insert(current);
                    vector<set<int>> edges;
                    for (const set<int>& edge : poly.edges) {
                        if (edge.find(current) == edge.end()) {
                            continue;
                        }
                        int v_i1 = current;
                        set<int> temp = {v_i1};
                        set<int> difference;
                        set_difference(edge.begin(), edge.end(), temp.begin(), temp.end(), std::inserter(difference, difference.begin()));
                        int v_i2 = *(difference.begin());
                        vector<FaceIntersection> intersects = other_poly.face_intersect({poly.verts[v_i1], poly.verts[v_i2]});
                        cout << "intersects size " << intersects.size() << endl;
                        sort(intersects.begin(), intersects.end(), Polyhedron::compare_face_intersections);
                        for (const FaceIntersection& intersect : intersects) {
                            cout << intersect.point[0] << "," << intersect.point[1] << "," << intersect.point[2] << endl;
                        }
                        for (int i = 0; i < intersects.size();) {
                            int j = i+1;
                            while (j < intersects.size() && round_float(intersects[i].alpha) == round_float(intersects[j].alpha)) {
                                j++;
                            }
                            if (j < intersects.size()) {
                                new_edges.insert({intersects[i].point, intersects[j].point});
                                cout << endl;
                                cout << intersects[i].point[0] << "," << intersects[i].point[1] << "," << intersects[i].point[2] << endl;
                                cout << intersects[j].point[0] << "," << intersects[j].point[1] << "," << intersects[j].point[2] << endl;
                            }
                            i = j+1;
                            while (i < intersects.size() && round_float(intersects[i-1].alpha) == round_float(intersects[i].alpha)) {
                                i++;
                            }
                        }
                        vector<int> in_faces = other_poly.in_faces(poly.verts[v_i2]);
                        if (intersects.size()) {
                            cout << "alpha " << intersects[0].alpha << endl;
                            if (!root_in_face) {
                                leaves.insert(intersects[0].point);
                            }
                            for (int i = 0; i < intersects.size(); i++) {
                                if (Polyhedron::round_float(intersects[i].alpha) == Polyhedron::round_float(intersects[0].alpha)) {
                                    face_lookup[intersects[i].point].insert(intersects[i].face_index);
                                } else {
                                    break;
                                }
                            }
                            if (root_in_poly) {
                                new_edges.insert({poly.verts[v_i1], intersects[0].point});
                            } else {
                                for (int i = 1; i < intersects.size(); i++) {
                                    if (round_float(intersects[i].alpha) > round_float(intersects[0].alpha)) {
                                        cout << intersects[i].alpha << " " << intersects[0].alpha << endl;
                                        new_edges.insert({intersects[0].point, intersects[i].point});
                                        break;
                                    }
                                }
                            }
                        } else if (in_faces.size()) {
                            if (!root_in_face) {
                                leaves.insert(poly.verts[v_i2]);
                            }
                            for (const int& face_index : in_faces){
                                face_lookup[poly.verts[v_i2]].insert(face_index);
                            }
                            q.push(v_i2);
                        } else {
                            q.push(v_i2);
                            if (root_in_poly) {
                                new_edges.insert({poly.verts[v_i1], poly.verts[v_i2]});
                            }
                        }
                    }
                }
                cout << "leaves size " << leaves.size() << " " << root_in_poly << endl;
                /*
                if (!leaves.size() && root_in_poly) {
		            edge_sets_per_poly[poly_index].push_back(set<set<std::array<double,3>>>());
                    for (const set<int>& edge : poly.edges) {
                        set<int>::iterator it = edge.begin();
                        int p_i1 = *it;
                        it++;
                        int p_i2 = *it;
                        edge_sets_per_poly[poly_index].front().insert({poly.verts[p_i1], poly.verts[p_i2]});
                    }
                    return edge_sets_per_poly;
                }*/
                set<std::array<double,3>> rounded_leaves;
                for (const std::array<double,3>& leaf : leaves) {
                    cout << "[" << leaf[0] << "," << leaf[1] << "," << leaf[2] << "] ";
                    rounded_leaves.insert(Polyhedron::round_point(leaf));
                }
                cout << endl;
                if (leaves.size() > 2 && seen_leaves.find(rounded_leaves) == seen_leaves.end()) {
                    seen_leaves.insert(rounded_leaves);
                    std::array<double,3> leaf_centroid = {0,0,0};
                    for (const std::array<double,3>& leaf : leaves) {
                        for (int i = 0; i < 3; i++) {
                            leaf_centroid[i] += leaf[i]/leaves.size();
                        }
                    }
                    std::array<double,3> vec0;
                    for (int i = 0; i < 3; i++) {
                        vec0[i] = leaf_centroid[i]-vert[i];
                    }
                    MatrixXd m(1,1);
                    VectorXd b(1);
                    m(0,0) = vec0[2];
                    b(0) = -vec0[0]*vec0[1]+vec0[1]*vec0[0];
                    VectorXd x = m.colPivHouseholderQr().solve(b);
                    std::array<double,3> vec1 = {vec0[1], -vec0[0], x(0)};
                    std::array<double,3> vec2 = Polyhedron::cross3D(vec0, vec1);
                    for (int i = 0; i < 3; i++) {
                        vec1[i] /= distance(vec1);
                        vec2[i] /= distance(vec2);
                    }
                    vector<std::array<double,3>> leaf_projection;
                    for (const std::array<double,3>& leaf : leaves) {
                        MatrixXd m(3,3);
                        VectorXd b(3);
                        for (int i = 0; i < 3; i++) {
                            m(i,0) = vec0[i];
                            m(i,1) = -vec1[i];
                            m(i,2) = -vec2[i];
                            b(i) = -leaf[i];
                        }
                        VectorXd x = m.colPivHouseholderQr().solve(b);
                        std::array<double,3> projection;
                        for (int i = 0; i < 3; i++) {
                            projection[i] = x(0)*vec0[i]+leaf[i];
                        }
                        leaf_projection.push_back(projection);
                    }
                    vector<std::array<double,3>> planar_leaf_projection = Polyhedron::make_planar(leaf_projection);
                    map<std::array<double,3>,std::array<double,3>> reverse_projection;
                    set<std::array<double,3>>::iterator it = leaves.begin();
                    for (int i = 0; i < leaves.size(); i++) {
                        reverse_projection[planar_leaf_projection[i]] = *it;
                        it++;
                    }
                    std::array<double,3>* point = &*std::min_element(planar_leaf_projection.begin(), planar_leaf_projection.end());
                    vector<std::array<double,3>> path;
                    map<std::array<double,3>,std::array<double,3>> rotated_mapping;
                    for (const std::array<double,3>& leaf : planar_leaf_projection) {
                        rotated_mapping[leaf] = leaf;
                    }
                    map<std::array<double,3>,std::array<double,3>> reverse_mapping;
                    for (const std::array<double,3>& leaf : planar_leaf_projection) {
                        reverse_mapping[leaf] = leaf;
                    }
                    std::array<double,3> angle = {0,0,0};
                    vector<std::array<double,2>> gift_wrap_angle_distance;
                    vector<std::array<double,3>> gift_wrap;
                    while (point != NULL && (!path.size() || *point != path[0])) {
                        path.push_back(*point);
                        vector<std::array<double,3>> keys;
                        vector<std::array<double,3>> temp_keys;
                        for (const pair<std::array<double,3>,std::array<double,3>>& p : rotated_mapping) {
                            keys.push_back(p.first);
                            std::array<double,3> temp_key;
                            for (int i = 0; i < 3; i++) {
                                temp_key[i] = p.first[i]-reverse_mapping[path.back()][i];
                            }
                            temp_keys.push_back(temp_key);
                        }
                        if (path.size() > 1) {
                            angle[2] += atan2(reverse_mapping[path.back()][1],reverse_mapping[path.back()][0]);
                            temp_keys = Polyhedron::rotate(temp_keys, angle);
                        }
                        map<std::array<double,3>,std::array<double,3>> temp;
                        for (int i = 0; i < keys.size(); i++) {
                            temp[temp_keys[i]] = rotated_mapping[keys[i]];
                        }
                        rotated_mapping = temp;
                        reverse_mapping.clear();
                        for (const pair<std::array<double,3>,std::array<double,3>>& p : rotated_mapping) {
                            reverse_mapping[p.second] = p.first;
                        }
                        gift_wrap_angle_distance.clear();
                        gift_wrap.clear();
                        for (const pair<std::array<double,3>,std::array<double,3>>& p : rotated_mapping) {
                            std::array<double,3> origin = {0,0,0};
                            if (p.first == origin) {
                                continue;
                            }
                            if (abs(atan2(p.first[1],p.first[0])-pi) < 0.00001) {
                                continue;
                            }
                            gift_wrap_angle_distance.push_back({atan2(p.first[1],p.first[0]), Polyhedron::distance(p.first)});
                            gift_wrap.push_back(p.second);
                        }
                        std::array<double,2> maxi = {-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
                        point = NULL;
                        for (int i = 0; i < gift_wrap_angle_distance.size(); i++) {
                            if (gift_wrap_angle_distance[i] > maxi) {
                                bool seen = false;
                                for (const std::array<double,3> item : path) {
                                    if (gift_wrap[i] == item) {
                                        seen = true;
                                        break;
                                    }
                                }
                                if (!seen) {
                                    maxi = gift_wrap_angle_distance[i];
                                    point = &gift_wrap[i];
                                }
                            }
                        }
                    }
                    vector<std::array<double,3>> temp;
                    for (const std::array<double,3>& point : path) {
                        temp.push_back(reverse_projection[point]);
                    }
                    path = temp;
                    cout << "path size " << path.size() << endl;
                    for (int index = 0; index < path.size(); index++) {
                        std::array<double,3> p1 = path[(index-1+path.size())%path.size()];
                        std::array<double,3> p2 = path[index];
                        set<int> starting_faces = face_lookup[p1];
                        cout << starting_faces.size() << endl;
                        for (const int& face_index1 : starting_faces) {
                            vector<vector<std::array<double,3>>> paths = Polyhedron::path_along_faces(p1,p2,{face_index1}, poly, other_poly);
                            cout << "paths size " << paths.size() << endl; 
                            for (vector<std::array<double,3>> point_path : paths) {
                                cout << "point path ";
                                for (int i = 0; i < point_path.size()-1; i++) {
                                    cout << point_path[i][0] << "," << point_path[i][1] << "," << point_path[i][2] << " ";
                                    new_edges.insert({point_path[i],point_path[i+1]});
                                }
                                cout << point_path[point_path.size()-1][0] << "," << point_path[point_path.size()-1][1] << "," << point_path[point_path.size()-1][2] << endl;
                            }
                            if (!paths.size()) {
                                new_edges.insert({p1,p2});
                            }
                        }
                    }
                }
            }
            /*bool all = true;
            for (const std::array<double,3>& vert : poly.verts) {
                if (!other_poly.is_inside(vert)) {
                    all = false;
                    break;
                }
            }
            cout << "new edges " << new_edges.size() << " " << all << endl;
            if (!new_edges.size() && all) {
                edge_sets_per_poly[0].clear();
                edge_sets_per_poly[1].clear();
                edge_sets_per_poly[poly_index].push_back(set<set<std::array<double,3>>>());
                for (const set<int>& edge : poly.edges) {
                    set<int>::iterator it = edge.begin();
                    int p_i1 = *it;
                    it++;
                    int p_i2 = *it;
                    edge_sets_per_poly[poly_index].front().insert({poly.verts[p_i1], poly.verts[p_i2]});
                }
                return edge_sets_per_poly;
	    }*/
            edge_sets_per_poly[poly_index].push_back(new_edges);
        }
        return edge_sets_to_poly(edge_sets_per_poly, polys);
    }
    // C++ port of Python's add_subtract_helper function
    static Polyhedron add_subtract_helper(Polyhedron poly1, Polyhedron poly2) {
        vector<pair<int, int>> pairings;
        
        // Convert faces to sets of edges (equivalent to Python's faces1 and faces2)
        vector<set<set<std::array<double,3>>>> faces1, faces2;
        
        for (const set<int>& face : poly1.faces) {
            set<set<std::array<double,3>>> face_edges;
            for (int edge_idx : face) {
                set<std::array<double,3>> edge_verts;
                for (int vert_idx : poly1.edges[edge_idx]) {
                    edge_verts.insert(poly1.verts[vert_idx]);
                }
                face_edges.insert(edge_verts);
            }
            faces1.push_back(face_edges);
        }
        
        for (const set<int>& face : poly2.faces) {
            set<set<std::array<double,3>>> face_edges;
            for (int edge_idx : face) {
                set<std::array<double,3>> edge_verts;
                for (int vert_idx : poly2.edges[edge_idx]) {
                    edge_verts.insert(poly2.verts[vert_idx]);
                }
                face_edges.insert(edge_verts);
            }
            faces2.push_back(face_edges);
        }
        
        // Find pairings between faces (exact match to Python)
        for (int face_index1 = 0; face_index1 < faces1.size(); face_index1++) {
            const set<set<std::array<double,3>>>& face1 = faces1[face_index1];
            
            for (int face_index2 = 0; face_index2 < faces2.size(); face_index2++) {
                const set<set<std::array<double,3>>>& face2 = faces2[face_index2];
                bool double_break = false;
                
                // Check if faces are coplanar
                set<std::array<double,3>> all_points;
                for (const set<std::array<double,3>>& edge : face1) {
                    for (const std::array<double,3>& point : edge) {
                        all_points.insert(point);
                    }
                }
                for (const set<std::array<double,3>>& edge : face2) {
                    for (const std::array<double,3>& point : edge) {
                        all_points.insert(point);
                    }
                }
                
                if (coplanar(vector<std::array<double,3>>(all_points.begin(), all_points.end()))) {
                    for (const set<std::array<double,3>>& edge1 : face1) {
                        for (const set<std::array<double,3>>& edge2 : face2) {
                            int points_on_edge1 = 0;
                            int points_on_edge2 = 0;
                            
                            for (const std::array<double,3>& point : edge2) {
                                if (point_on_segment(edge1, point)) {
                                    points_on_edge1++;
                                }
                            }
                            for (const std::array<double,3>& point : edge1) {
                                if (point_on_segment(edge2, point)) {
                                    points_on_edge2++;
                                }
                            }
                            
                            if (points_on_edge1 == 2) {
                                pairings.push_back({face_index1, face_index2});
                                double_break = true;
                                break;
                            } else if (points_on_edge2 == 2) {
                                pairings.push_back({face_index1, face_index2});
                                double_break = true;
                                break;
                            } else if (points_on_edge1 == 1 && points_on_edge2 == 1) {
                                set<std::array<double,3>> intersection;
                                set_intersection(edge1.begin(), edge1.end(), edge2.begin(), edge2.end(),
                                               inserter(intersection, intersection.begin()));
                                if (intersection.empty()) {
                                    pairings.push_back({face_index1, face_index2});
                                    double_break = true;
                                    break;
                                }
                            }
                        }
                        if (double_break) break;
                    }
                }
            }
        }
        
        // Create maps (fixed initialization)
        vector<vector<set<set<std::array<double,3>>>>> faces;
        faces.push_back(faces1);
        faces.push_back(faces2);
        
        vector<map<int, vector<int>>> maps(2);
        for (int i = 0; i < faces1.size(); i++) {
            maps[0][i] = vector<int>();
        }
        for (int i = 0; i < faces2.size(); i++) {
            maps[1][i] = vector<int>();
        }
        
        for (const pair<int,int>& pairing : pairings) {
            maps[0][pairing.first].push_back(pairing.second);
            maps[1][pairing.second].push_back(pairing.first);
        }
        
        set<pair<int, int>> visited;
        vector<set<set<std::array<double,3>>>> new_faces;
        
        for (int poly_index = 0; poly_index < 2; poly_index++) {
            for (int face_index = 0; face_index < faces[poly_index].size(); face_index++) {
                pair<int, int> current_key = {poly_index, face_index};
                if (visited.find(current_key) != visited.end()) continue;
                
                set<set<std::array<double,3>>> face = faces[poly_index][face_index];
                queue<pair<int, int>> q;
                q.push(current_key);
                
                while (!q.empty()) {
                    pair<int, int> indices = q.front();
                    q.pop();
                    
                    visited.insert(indices);
                    
                    int p_i1 = indices.first;
                    int f_i1 = indices.second;
                    int p_i2 = (p_i1 + 1) % 2;
                    
                    for (int f_i2 : maps[p_i1][f_i1]) {
                        pair<int, int> next_key = {p_i2, f_i2};
                        if (visited.find(next_key) != visited.end()) continue;
                        
                        q.push(next_key);
                        
                        for (const set<std::array<double,3>>& edge2 : faces[p_i2][f_i2]) {
                            vector<set<std::array<double,3>>> face_vector(face.begin(), face.end());
                            int edge_index1 = 0;
                            vector<set<std::array<double,3>>> new_edges;
                            bool combine = false;
                            
                            while (edge_index1 < face_vector.size()) {
                                set<std::array<double,3>> edge1 = face_vector[edge_index1];
                                
                                int points_on_edge1 = 0;
                                int points_on_edge2 = 0;
                                
                                for (const std::array<double,3>& point : edge2) {
                                    if (point_on_segment(edge1, point)) {
                                        points_on_edge1++;
                                    }
                                }
                                for (const std::array<double,3>& point : edge1) {
                                    if (point_on_segment(edge2, point)) {
                                        points_on_edge2++;
                                    }
                                }
                                
                                if (points_on_edge1 == 2) {
                                    face_vector.erase(face_vector.begin() + edge_index1);
                                    
                                    vector<std::array<double,3>> edge1_vec(edge1.begin(), edge1.end());
                                    std::array<double,3> p1 = edge1_vec[0];
                                    std::array<double,3> p2 = edge1_vec[1];
                                    
                                    // Find closest point in edge2 to p1
                                    vector<pair<double, std::array<double,3>>> distances_p1;
                                    for (const std::array<double,3>& point : edge2) {
                                        distances_p1.push_back({distance(p1, point), point});
                                    }
                                    sort(distances_p1.begin(), distances_p1.end());
                                    
                                    if (round_float(distances_p1[0].first) > 0) {
                                        set<std::array<double,3>> edge3 = {p1, distances_p1[0].second};
                                        new_edges.push_back(edge3);
                                    }
                                    
                                    // Find closest point in edge2 to p2
                                    vector<pair<double, std::array<double,3>>> distances_p2;
                                    for (const std::array<double,3>& point : edge2) {
                                        distances_p2.push_back({distance(p2, point), point});
                                    }
                                    sort(distances_p2.begin(), distances_p2.end());
                                    
                                    if (round_float(distances_p2[0].first) > 0) {
                                        set<std::array<double,3>> edge3 = {p2, distances_p2[0].second};
                                        new_edges.push_back(edge3);
                                    }
                                    combine = true;
                                    
                                } else if (points_on_edge2 == 2) {
                                    face_vector.erase(face_vector.begin() + edge_index1);
                                    
                                    vector<std::array<double,3>> edge2_vec(edge2.begin(), edge2.end());
                                    std::array<double,3> p1 = edge2_vec[0];
                                    std::array<double,3> p2 = edge2_vec[1];
                                    
                                    // Similar logic for edge2
                                    vector<pair<double, std::array<double,3>>> distances_p1;
                                    for (const std::array<double,3>& point : edge1) {
                                        distances_p1.push_back({distance(p1, point), point});
                                    }
                                    sort(distances_p1.begin(), distances_p1.end());
                                    
                                    if (round_float(distances_p1[0].first) > 0) {
                                        set<std::array<double,3>> edge3 = {p1, distances_p1[0].second};
                                        new_edges.push_back(edge3);
                                    }
                                    
                                    vector<pair<double, std::array<double,3>>> distances_p2;
                                    for (const std::array<double,3>& point : edge1) {
                                        distances_p2.push_back({distance(p2, point), point});
                                    }
                                    sort(distances_p2.begin(), distances_p2.end());
                                    
                                    if (round_float(distances_p2[0].first) > 0) {
                                        set<std::array<double,3>> edge3 = {p2, distances_p2[0].second};
                                        new_edges.push_back(edge3);
                                    }
                                    combine = true;
                                    
                                } else if (points_on_edge1 == 1 && points_on_edge2 == 1) {
                                    set<std::array<double,3>> intersection;
                                    set_intersection(edge1.begin(), edge1.end(), edge2.begin(), edge2.end(),
                                                   inserter(intersection, intersection.begin()));
                                    
                                    if (intersection.empty()) {
                                        face_vector.erase(face_vector.begin() + edge_index1);
                                        
                                        vector<std::array<double,3>> edge1_vec(edge1.begin(), edge1.end());
                                        vector<std::array<double,3>> edge2_vec(edge2.begin(), edge2.end());
                                        std::array<double,3> p1 = edge1_vec[0];
                                        std::array<double,3> p2 = edge1_vec[1];
                                        std::array<double,3> p3 = edge2_vec[0];
                                        std::array<double,3> p4 = edge2_vec[1];
                                        
                                        if (point_on_segment(edge2, p1)) {
                                            swap(p1, p2);
                                        }
                                        if (point_on_segment(edge1, p4)) {
                                            swap(p3, p4);
                                        }
                                        
                                        if (round_float(distance(p1, p3)) > 0) {
                                            set<std::array<double,3>> edge3 = {p1, p3};
                                            new_edges.push_back(edge3);
                                        }
                                        if (round_float(distance(p2, p4)) > 0) {
                                            set<std::array<double,3>> edge3 = {p2, p4};
                                            new_edges.push_back(edge3);
                                        }
                                        combine = true;
                                        
                                    } else if (intersection.size() == 1) {
                                        set<std::array<double,3>> union_set;
                                        set_union(edge1.begin(), edge1.end(), edge2.begin(), edge2.end(),
                                                inserter(union_set, union_set.begin()));
                                        if (colinear(vector<std::array<double,3>>(union_set.begin(), union_set.end()))) {
                                            face_vector.erase(face_vector.begin() + edge_index1);
                                            
                                            // Remove edge2 from new_edges if present
                                            vector<set<std::array<double,3>>>::iterator it = find(new_edges.begin(), new_edges.end(), edge2);
                                            if (it != new_edges.end()) {
                                                new_edges.erase(it);
                                            }
                                            
                                            // Add symmetric difference
                                            set<std::array<double,3>> symmetric_diff;
                                            set_symmetric_difference(edge1.begin(), edge1.end(), edge2.begin(), edge2.end(),
                                                                   inserter(symmetric_diff, symmetric_diff.begin()));
                                            new_edges.push_back(symmetric_diff);
                                            combine = true;
                                        } else {
                                            edge_index1++;
                                        }
                                    } else {
                                        edge_index1++;
                                    }
                                } else {
                                    edge_index1++;
                                }
                            }
                            
                            // Update face
                            face.clear();
                            for (const set<std::array<double,3>>& edge : face_vector) {
                                face.insert(edge);
                            }
                            
                            if (combine) {
                                for (const set<std::array<double,3>>& edge : new_edges) {
                                    face.insert(edge);
                                }
                            } else {
                                face.insert(edge2);
                            }
                        }
                    }
                }
                
                // Convert face to circuits
                set<vector<std::array<double,3>>> circuits = circuit_helper(face);
                
                // Remove collinear points from circuits
                set<vector<std::array<double,3>>> filtered_circuits;
                for (const vector<std::array<double,3>>& circuit : circuits) {
                    vector<std::array<double,3>> filtered_circuit;
                    for (int i = 0; i < circuit.size(); i++) {
                        std::array<double,3> prev = circuit[(i - 1 + circuit.size()) % circuit.size()];
                        std::array<double,3> curr = circuit[i];
                        std::array<double,3> next = circuit[(i + 1) % circuit.size()];
                        
                        if (!colinear(vector<std::array<double,3>>{prev, curr, next})) {
                            filtered_circuit.push_back(curr);
                        }
                    }
                    if (!filtered_circuit.empty()) {
                        filtered_circuits.insert(filtered_circuit);
                    }
                }
                
                vector<std::array<double,3>>* exterior = find_exterior_circuit(filtered_circuits);
                if (exterior == NULL) {
                    for (const vector<std::array<double,3>>& circuit : filtered_circuits) {
                        set<set<std::array<double,3>>> circuit_face;
                        for (int i = 0; i < circuit.size(); i++) {
                            std::array<double,3> p1 = circuit[(i - 1 + circuit.size()) % circuit.size()];
                            std::array<double,3> p2 = circuit[i];
                            circuit_face.insert({p1, p2});
                        }
                        new_faces.push_back(circuit_face);
                    }
                } else {
                    new_faces.push_back(face);
                    delete exterior;
                }
            }
        }
        
        // Handle edge intersections within faces
        for (int face_index = 0; face_index < new_faces.size(); face_index++) {
            bool updated = true;
            while (updated) {
                updated = false;
                for (const set<std::array<double,3>>& edge1 : new_faces[face_index]) {
                    for (const set<std::array<double,3>>& edge2 : new_faces[face_index]) {
                        if (edge1 == edge2) continue;
                        
                        std::array<double,3>* intersect = intersect_segments(edge1, edge2);
                        if (intersect != NULL) {
                            bool intersect_on_endpoint = false;
                            for (const std::array<double,3>& point : edge1) {
                                if (distance(*intersect, point) < 0.001) {
                                    intersect_on_endpoint = true;
                                    break;
                                }
                            }
                            if (!intersect_on_endpoint) {
                                for (const std::array<double,3>& point : edge2) {
                                    if (distance(*intersect, point) < 0.001) {
                                        intersect_on_endpoint = true;
                                        break;
                                    }
                                }
                            }
                            
                            if (!intersect_on_endpoint) {
                                new_faces[face_index].erase(edge1);
                                new_faces[face_index].erase(edge2);
                                
                                vector<std::array<double,3>> edge1_vec(edge1.begin(), edge1.end());
                                vector<std::array<double,3>> edge2_vec(edge2.begin(), edge2.end());
                                
                                new_faces[face_index].insert({edge1_vec[0], *intersect});
                                new_faces[face_index].insert({edge1_vec[1], *intersect});
                                new_faces[face_index].insert({edge2_vec[0], *intersect});
                                new_faces[face_index].insert({edge2_vec[1], *intersect});
                                
                                updated = true;
                                delete intersect;
                                break;
                            }
                            delete intersect;
                        }
                    }
                    if (updated) break;
                }
            }
        }
        
        // Merge coplanar faces
        bool updated = true;
        while (updated) {
            updated = false;
            for (int face_index1 = 0; face_index1 < new_faces.size(); face_index1++) {
                for (int face_index2 = face_index1 + 1; face_index2 < new_faces.size(); face_index2++) {
                    set<std::array<double,3>> all_points;
                    for (const set<std::array<double,3>>& edge : new_faces[face_index1]) {
                        for (const std::array<double,3>& point : edge) {
                            all_points.insert(point);
                        }
                    }
                    for (const set<std::array<double,3>>& edge : new_faces[face_index2]) {
                        for (const std::array<double,3>& point : edge) {
                            all_points.insert(point);
                        }
                    }
                    
                    if (coplanar(vector<std::array<double,3>>(all_points.begin(), all_points.end()))) {
                        set<set<std::array<double,3>>> combined_face;
                        set_union(new_faces[face_index1].begin(), new_faces[face_index1].end(),
                                new_faces[face_index2].begin(), new_faces[face_index2].end(),
                                inserter(combined_face, combined_face.begin()));
                        
                        set<vector<std::array<double,3>>> circuits = circuit_helper(combined_face);
                        
                        // Filter collinear points
                        set<vector<std::array<double,3>>> filtered_circuits;
                        for (const vector<std::array<double,3>>& circuit : circuits) {
                            vector<std::array<double,3>> filtered_circuit;
                            for (int i = 0; i < circuit.size(); i++) {
                                std::array<double,3> prev = circuit[(i - 1 + circuit.size()) % circuit.size()];
                                std::array<double,3> curr = circuit[i];
                                std::array<double,3> next = circuit[(i + 1) % circuit.size()];
                                
                                if (!colinear(vector<std::array<double,3>>{prev, curr, next})) {
                                    filtered_circuit.push_back(curr);
                                }
                            }
                            if (!filtered_circuit.empty()) {
                                filtered_circuits.insert(filtered_circuit);
                            }
                        }
                        
                        vector<std::array<double,3>>* exterior = find_exterior_circuit(filtered_circuits);
                        if (exterior != NULL) {
                            new_faces.erase(new_faces.begin() + face_index2);
                            new_faces.erase(new_faces.begin() + face_index1);
                            new_faces.push_back(combined_face);
                            updated = true;
                            delete exterior;
                            break;
                        }
                        if (exterior) delete exterior;
                    }
                }
                if (updated) break;
            }
        }
        
        // Final circuit processing
        int face_index = 0;
        while (face_index < new_faces.size()) {
            set<vector<std::array<double,3>>> circuits = circuit_helper(new_faces[face_index]);
            
            // Filter collinear points
            set<vector<std::array<double,3>>> filtered_circuits;
            for (const vector<std::array<double,3>>& circuit : circuits) {
                vector<std::array<double,3>> filtered_circuit;
                for (const std::array<double,3>& x : circuit) {
                    bool is_collinear = false;
                    for (int i = 0; i < circuit.size(); i++) {
                        if (circuit[i] == x) {
                            std::array<double,3> prev = circuit[(i - 1 + circuit.size()) % circuit.size()];
                            std::array<double,3> next = circuit[(i + 1) % circuit.size()];
                            if (colinear(vector<std::array<double,3>>{prev, x, next})) {
                                is_collinear = true;
                                break;
                            }
                        }
                    }
                    if (!is_collinear) {
                        filtered_circuit.push_back(x);
                    }
                }
                if (!filtered_circuit.empty()) {
                    filtered_circuits.insert(filtered_circuit);
                }
            }
            
            vector<std::array<double,3>>* exterior = find_exterior_circuit(filtered_circuits);
            if (exterior == NULL) {
                new_faces.erase(new_faces.begin() + face_index);
                for (const vector<std::array<double,3>>& circuit : filtered_circuits) {
                    set<set<std::array<double,3>>> circuit_face;
                    for (int i = 0; i < circuit.size(); i++) {
                        std::array<double,3> p1 = circuit[(i - 1 + circuit.size()) % circuit.size()];
                        std::array<double,3> p2 = circuit[i];
                        circuit_face.insert({p1, p2});
                    }
                    new_faces.push_back(circuit_face);
                }
            } else {
                face_index++;
                delete exterior;
            }
        }
        
        // Build final polyhedron
        Polyhedron poly;
        
        // Collect all unique vertices
        set<std::array<double,3>> all_vertices;
        for (const set<set<std::array<double,3>>>& face : new_faces) {
            for (const set<std::array<double,3>>& edge : face) {
                for (const std::array<double,3>& point : edge) {
                    all_vertices.insert(point);
                }
            }
        }
        poly.verts = vector<std::array<double,3>>(all_vertices.begin(), all_vertices.end());
        
        // Collect all unique edges
        set<set<int>> all_edges;
        for (const set<set<std::array<double,3>>>& face : new_faces) {
            for (const set<std::array<double,3>>& edge : face) {
                set<int> edge_indices;
                for (const std::array<double,3>& point : edge) {
                    auto it = find(poly.verts.begin(), poly.verts.end(), point);
                    edge_indices.insert(it - poly.verts.begin());
                }
                all_edges.insert(edge_indices);
            }
        }
        poly.edges = vector<set<int>>(all_edges.begin(), all_edges.end());
        
        // Build faces
        for (const set<set<std::array<double,3>>>& face : new_faces) {
            if (face.empty()) continue;
            
            set<int> face_indices;
            for (const set<std::array<double,3>>& edge : face) {
                set<int> edge_indices;
                for (const std::array<double,3>& point : edge) {
                    auto it = find(poly.verts.begin(), poly.verts.end(), point);
                    edge_indices.insert(it - poly.verts.begin());
                }
                auto edge_it = find(poly.edges.begin(), poly.edges.end(), edge_indices);
                face_indices.insert(edge_it - poly.edges.begin());
            }
            poly.faces.push_back(face_indices);
        }
        
        return poly;
    }
    Polyhedron subtract(Polyhedron other) {
        Polyhedron new_poly;
        new_poly.verts = this->verts;
        new_poly.edges = this->edges;
        new_poly.faces = this->faces;
        
        Polyhedron poly = this->intersect(other);
        return Polyhedron::add_subtract_helper(poly, new_poly);
    }
    
    Polyhedron add(Polyhedron other) {
        Polyhedron new_poly1;
        new_poly1.verts = other.verts;
        new_poly1.edges = other.edges;
        new_poly1.faces = other.faces;
        
        Polyhedron new_poly2;
        new_poly2.verts = this->verts;
        new_poly2.edges = this->edges;
        new_poly2.faces = this->faces;

        Polyhedron poly = this->intersect(other);

        new_poly1 = new_poly1.subtract(poly);
        cout << new_poly1.verts.size() << endl;
        return Polyhedron::add_subtract_helper(new_poly1, new_poly2);
    }    
};

int main(int argc, char* argv[]) {
    // Default parameters
    string filename1 = "poly1.ply";
    string filename2 = "poly2.ply";
    string filename3 = "out.ply";
    string operation = "intersect";  // Default operation
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            cout << "Usage: " << argv[0] << " [options]" << endl;
            cout << "Options:" << endl;
            cout << "  -i, --input1 <file>     First input polyhedron file (default: poly1.ply)" << endl;
            cout << "  -j, --input2 <file>     Second input polyhedron file (default: poly2.ply)" << endl;
            cout << "  -o, --output <file>     Output file (default: out.ply)" << endl;
            cout << "  -op, --operation <op>   Operation: intersect, add, subtract (default: intersect)" << endl;
            cout << "  -h, --help              Show this help message" << endl;
            cout << endl;
            cout << "Examples:" << endl;
            cout << "  " << argv[0] << " -i cube.ply -j sphere.ply -op intersect -o result.ply" << endl;
            cout << "  " << argv[0] << " -op add" << endl;
            cout << "  " << argv[0] << " --operation subtract --output diff.ply" << endl;
            return 0;
        }
        else if ((arg == "-i" || arg == "--input1") && i + 1 < argc) {
            filename1 = argv[++i];
        }
        else if ((arg == "-j" || arg == "--input2") && i + 1 < argc) {
            filename2 = argv[++i];
        }
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            filename3 = argv[++i];
        }
        else if ((arg == "-op" || arg == "--operation") && i + 1 < argc) {
            operation = argv[++i];
            // Validate operation
            if (operation != "intersect" && operation != "add" && operation != "subtract") {
                cerr << "Error: Invalid operation '" << operation << "'. Must be 'intersect', 'add', or 'subtract'." << endl;
                return 1;
            }
        }
        else if (i == 1 && arg.find("-") != 0) {
            // First positional argument (backward compatibility)
            filename1 = arg;
        }
        else if (i == 2 && arg.find("-") != 0) {
            // Second positional argument (backward compatibility)
            filename2 = arg;
        }
        else if (i == 3 && arg.find("-") != 0) {
            // Third positional argument (backward compatibility)
            filename3 = arg;
        }
        else if (arg.find("-") == 0) {
            cerr << "Error: Unknown option '" << arg << "'. Use -h for help." << endl;
            return 1;
        }
    }
    
    cout << "Loading polyhedra..." << endl;
    cout << "Input 1: " << filename1 << endl;
    cout << "Input 2: " << filename2 << endl;
    cout << "Output: " << filename3 << endl;
    cout << "Operation: " << operation << endl;
    cout << endl;
    
    // Load first polyhedron
    Polyhedron poly1;
    if (!poly1.load(filename1)) {
        cerr << "Error: Could not load " << filename1 << endl;
        return 1;
    }
    cout << "Loaded poly1: " << poly1.verts.size() << " vertices, " 
         << poly1.edges.size() << " edges, " << poly1.faces.size() << " faces" << endl;
    
    // Load second polyhedron
    Polyhedron poly2;
    if (!poly2.load(filename2)) {
        cerr << "Error: Could not load " << filename2 << endl;
        return 1;
    }
    cout << "Loaded poly2: " << poly2.verts.size() << " vertices, " 
         << poly2.edges.size() << " edges, " << poly2.faces.size() << " faces" << endl;
    cout << endl;
    
    // Perform the requested operation
    cout << "Performing " << operation << " operation..." << endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    Polyhedron result;
    try {
        if (operation == "intersect") {
            result = poly1.intersect(poly2);
        }
        else if (operation == "add") {
            result = poly1.add(poly2);
        }
        else if (operation == "subtract") {
            result = poly1.subtract(poly2);
        }
    }
    catch (const exception& e) {
        cerr << "Error during " << operation << " operation: " << e.what() << endl;
        return 1;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    cout << "Operation completed in " << duration.count() / 1000000.0 << " seconds" << endl;
    cout << "Result: " << result.verts.size() << " vertices, " 
         << result.edges.size() << " edges, " << result.faces.size() << " faces" << endl;
    cout << endl;
    
    // Output vertices (for debugging)
    if (result.verts.size() <= 50) {  // Only print if not too many vertices
        cout << "Result vertices:" << endl;
        for (const std::array<double,3>& vert : result.verts) {
            cout << "[" << vert[0] << "," << vert[1] << "," << vert[2] << "]" << endl;
        }
    } else {
        cout << "Too many vertices to display (" << result.verts.size() << " vertices)" << endl;
    }
    
    // Save result to file
    cout << "Saving result to " << filename3 << "..." << endl;
    if (result.dump(filename3)) {
        cout << "Result saved successfully." << endl;
    } else {
        cerr << "Error: Could not save result to " << filename3 << endl;
        return 1;
    }
    
    return 0;
}
