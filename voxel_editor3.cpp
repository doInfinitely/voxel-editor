#include <SDL.h>
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
            set<std::array<double,3>> union_set;
            for (const std::array<double,3>& x : path) {
                union_set.insert(Polyhedron::round_point(x));
            }
            for (const int& x : y) {
                union_set.insert(Polyhedron::round_point(verts[x]));
            }
            if (Polyhedron::coplanar(union_set)) {
                vector<set<int>>::iterator it = find(edges.begin(),edges.end(),y);
                current = std::distance(edges.begin(), it);
                old_circuits.insert(output.begin(),output.end());
                set<vector<std::array<double,3>>> intermediate = this->circuits(face_index, start, previous, current, path, old_circuits);
                output.insert(intermediate.begin(),intermediate.end());
            }
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
        cout << "circuit_sizes ";
        for (const vector<std::array<double,3>>& circuit : output) {
            cout << circuit.size() << " ";
        }
        cout << endl;
        vector<std::array<double,3>>* exterior_circuit = Polyhedron::find_exterior_circuit(output);
        cout << "NO_EXTERIOR_CIRCUIT1" << endl; 
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
        cout << "NO_EXTERIOR_CIRCUIT2" << endl; 
        delete exterior_circuit;
        return output;
    }
    static bool colinear(vector<std::array<double,3>> points) {
        if (points.size() < 3) {
            return true;
        }
        MatrixXd m(3,1);
        for (int i=0; i < 3; i++) {
            m(i,0) = points[1][i]-points[0][i];
        }
        for (const std::array<double,3>& p : points) {
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
        }
        for (int i = 0; i < 3; i++) {
            m(3,i) = 1;
        }
        b(3) = 1;
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        std::array<double,3> point_prime;
        for (int i=0; i < 3; i++) {
            point_prime[i] = b_prime(i);
        }
        //cout << Polyhedron::distance(point,point_prime) << " " << round_float(x(0)) << " " << round_float(x(1)) << " " << round_float(x(2)) << " " << (round_float(x(0)) >= 0 && round_float(x(1)) >= 0 && round_float(x(2)) >= 0) << endl;
        if (Polyhedron::distance(point,point_prime) > 0.000001) {
            return false;
        }
        if (round_float(x(0)) >= 0 && round_float(x(1)) >= 0 && round_float(x(2)) >= 0) {
            return true;
        }
        return false;
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
        std::array<double,3> p1 = {Polyhedron::distance(circuit[0],circuit[1]),0,0};
        std::array<double,3> angle = {0,0,asin((circuit[1][1]-circuit[0][1])/p1[0])};
        p1 = Polyhedron::rotate({p1}, angle)[0];
        angle[1] = asin(min(max((circuit[1][2]-circuit[0][2])/p1[0],-1.0),1.0));
        p1 = Polyhedron::rotate({p1}, {0,angle[1],0})[0];
        angle[1] = -angle[1];
        angle[2] = -angle[2];
        p1 = Polyhedron::rotate({{circuit[1][0]-circuit[0][0],circuit[1][1]-circuit[0][1],circuit[1][2]-circuit[0][2]}}, angle)[0];
        std::array<double,3> p2 = Polyhedron::rotate({{circuit[2][0]-circuit[0][0],circuit[2][1]-circuit[0][1],circuit[2][2]-circuit[0][2]}}, angle)[0];
        angle[0] = -atan2(p2[2]-p1[2],p2[1]-p1[1]);
        std::array<double,3> start = {circuit[0][0],circuit[0][1],circuit[0][2]};
        for (std::array<double,3>& x : circuit) {
            for (int i = 0; i < 3; i++) {
                x[i] -= start[i];
            }
            //cout << "[" << x[0] << "," << x[1] << "," << x[2] << "] ";
        }
        //cout << "make_planr" <<endl;
        vector<std::array<double,3>> output = Polyhedron::rotate(circuit, {0,angle[1],angle[2]});
        output = Polyhedron::rotate(output, {angle[0],0,0});
        for (std::array<double,3>& x : output) {
            x[2] = 0;
        }
        //cout << "[" << angle[0] << "," << angle[1] << "," << angle[2] << "] " << endl;
        //for (std::array<double,3>& x : output) {
        //    cout << "[" << x[0] << "," << x[1] << "," << x[2] << "] ";
        //}
        //cout << "make_planar" <<endl;
        return output;
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
            if (Polyhedron::is_clockwise(Polyhedron::make_planar(circuit))) {
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
        vector<std::array<double,3>> planar = Polyhedron::make_planar(circuit);
        //for (int i = 0; i < planar.size(); i++) {
        //    cout << "[" << planar[i][0] << "," << planar[i][1] << "," << planar[i][2] << "] ";
        //}
        //cout << endl;
        //cout << Polyhedron::is_clockwise(planar) << endl;
        if (!Polyhedron::is_clockwise(planar)) {
            std::reverse(planar.begin(),planar.end());
            std::reverse(circuit.begin(),circuit.end());
        }
        vector<bool> is_convex = Polyhedron::clockwise_angles(planar);
        for (int i = 0; i < is_convex.size(); i++) {
            cout << is_convex[i] ;
        }
        cout << endl;
        for (int i = 0; i < circuit.size(); i++) {
            if (is_convex[i]) {
                std::array<std::array<double,3>,3> triangle = {circuit[(i-1+circuit.size())%circuit.size()],circuit[i],circuit[(i+1)%circuit.size()]};
                bool no_break = true;
                for (int j = 0; j < circuit.size(); j++) {
                    if (j != i && j != (i+1)%circuit.size() && j != (i-1+circuit.size())%circuit.size()) {

                        if (Polyhedron::inside_triangle(triangle,circuit[j]) && round_point(circuit[j]) != round_point(circuit[(i-1+circuit.size())%circuit.size()]) && round_point(circuit[j]) != round_point(circuit[i]) && round_point(circuit[j]) != round_point(circuit[(i+1)%circuit.size()])) {
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
        cout << "triangulate start" << endl;
        while (remainder.size() >3) {
            for (int i = 0; i < remainder.size();) {
                if (Polyhedron::colinear(set<std::array<double,3>>{remainder[(i-1+remainder.size())%remainder.size()],remainder[i],remainder[(i+1)%remainder.size()]})) {
                    cout << "colinear " << remainder[i][0] << " " << remainder[i][1] << " " << remainder[i][2] << endl;
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
                cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
            }
            cout << endl;
        }
        output.push_back({remainder[0],remainder[1],remainder[2]});
        for (const std::array<double,3>& point : remainder) {
            cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
        }
        cout << endl;
        cout << "triangulate end" << endl;
        return output;
    }
    static vector<std::array<double,3>>* find_exterior_circuit(set<vector<std::array<double,3>>> circuits) {
        vector<vector<std::array<double,3>>> circuits_list(circuits.begin(),circuits.end());
        for (int i = 0; i < circuits_list.size(); i++) {
            vector<std::array<std::array<double,3>,3>> triangulation = Polyhedron::triangulate(circuits_list[i]);
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
    static bool intersections_comp(Polyhedron::CircuitIntersection a, Polyhedron::CircuitIntersection b) {
        return a.alpha < b.alpha;
    }
    static vector<std::array<double,3>>* circuit_cut(set<vector<std::array<double,3>>> circuits) {
        cout << "circuits" << circuits.size() << endl;
        for (const vector<std::array<double,3>>& circuit : circuits) {
            for (const std::array<double,3> point : circuit) {
                cout << point.size() << endl;
                cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
            }
            cout << endl;
        }
        cout << endl;
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
            sort(intersections.begin(),intersections.end(), intersections_comp);
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
    bool is_inside(std::array<double,3> point) {
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*Polyhedron::circuit_cut(this->circuits(face_index)))) {
                if (Polyhedron::inside_triangle(triangle, point)) {
                    return true;
                }
            }
        }
        std::array<double,3> vec = {(double)(rand()/RAND_MAX),(double)(rand()/RAND_MAX),(double)(rand()/RAND_MAX)};
        vector<Polyhedron::IsInsideIntersection> output;
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            vector<std::array<double,3>> circuit = *Polyhedron::circuit_cut(this->circuits(face_index));
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(circuit)) {
                MatrixXd m(4,4);
                VectorXd b(4);
                for (int i = 0; i < 3; i++) {
                    m(i,0) = triangle[0][i];
                    m(i,1) = triangle[1][i];
                    m(i,2) = triangle[2][i];
                    m(i,3) = -vec[i];
                    b(i) = point[i] ;
                }
                for (int i = 0; i < 3; i++) {
                    m(3,i) = 1;
                }
                m(3,3) = 0;
                b(3) = 1;
                VectorXd x = m.colPivHouseholderQr().solve(b);
                VectorXd b_prime = m*x;
                bool error = false;
                for (int i = 0; i < 4; i++) {
                    if ( abs(b_prime(i)-b(i)) > 0.1*abs(b(i)) ) {
                        error = true;
                    }
                }
                if (error) {
                    continue;
                }
                if (x(0) >= 0 && x(1) >= 0 && x(2) >= 0 && x(3) > 0) {
                    std::array<double,3> p;
                    for (int i = 0; i < 3; i++) {
                        p[i] = x(0)*triangle[0][i]+x(1)*triangle[1][i]+x(2)*triangle[2][i];
                    }
                    if (Polyhedron::inside_triangle(triangle,p)) {
                        Polyhedron::IsInsideIntersection intersection;
                        intersection.gamma = x(3);
                        intersection.point = p;
                        intersection.face_index = face_index;
                        output.push_back(intersection);
                        break;
                    }
                }
            }
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
};

class Camera {
  public:
    std::array<double,3> origin = {0,0,0};
    std::array<double,3> x_vector = {1,0,0};
    std::array<double,3> y_vector = {0,1,0};
    std::array<double,3> focal = {0,0,-1};
    double zoom = 1;
    std::array<double,2> project(std::array<double,3> point) {
        MatrixXd m(3,3);
        VectorXd b(3);
        for (int i = 0; i < 3; i++) {
            m(i,0) = x_vector[i];
            m(i,1) = y_vector[i];
            m(i,2) = point[i]-focal[i];
            b(i) = point[i]-origin[i];
        }
        VectorXd x = m.lu().solve(b);
        std::array<double,2> output = {x(0)*zoom,x(1)*zoom};
        return output;
    }
    void move(std::array<double,3> displacement) {
        origin[0] += displacement[0];
        origin[1] += displacement[1];
        origin[2] += displacement[2];
        focal[0] += displacement[0];
        focal[1] += displacement[1];
        focal[2] += displacement[2];
    }
    void rotate(std::array<double,3> angles) {
        focal[0] -= origin[0];
        focal[1] -= origin[1];
        focal[2] -= origin[2];
        vector<std::array<double,3>> camera_items = Polyhedron::rotate({focal, x_vector, y_vector}, angles);
        focal = camera_items[0];
        focal[0] += origin[0];
        focal[1] += origin[1];
        focal[2] += origin[2];
        x_vector = camera_items[1];
        y_vector = camera_items[2];
    }
    std::array<double,3> forward_vector() {
        std::array<double,3> vec = {origin[0]-focal[0],origin[1]-focal[1],origin[2]-focal[2]};
        double mag = sqrt(Polyhedron::dot(vec,vec));
        return {vec[0]/mag, vec[1]/mag, vec[2]/mag};
    }
};
class Box {
  public:
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
    Box(std::array<double,2> x, std::array<double,2> y, std::array<double,2> z) {
        x_min = x[0];
        x_max = x[1];
        y_min = y[0];
        y_max = y[1];
        z_min = z[0];
        z_max = z[1];
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
        MatrixXd m(3,2);
        VectorXd b(3);
        for (int i = 0; i < 3; i++) {
            m(i,0) = p1[i]-p2[i];
            m(i,1) = p4[i]-p3[i];
            b(i) = p4[i]-p2[i];
        }
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        for (int i=0; i < 3; i++) {
            if ( abs(b_prime(i)-b(i)) > 0.1*abs(b(i)) ) {
                return NULL;
            }
        }
        if (x(0) >= 0 && x(0) <= 1 && x(1) >= 0 && x(1) <= 1) {
            std::array<double,3>* output = new std::array<double,3>;
            for (int i = 0; i < 3; i++) {
               (*output)[i] = x(0)*p1[i]+(1-x(0))*p2[i];
            }
            return output;
        }
        return NULL;
    }
    bool is_inside(std::array<double,3> point) const {
        return point[0] >= x_min && point[0] <= x_max && point[1] >= y_min && point[1] <= y_max && point[2] >= z_min && point[2] <= z_max;
    }
    bool intersect(set<std::array<double,3>> edge) const {
        set<std::array<double,3>>::iterator it = edge.begin();
        std::array<double,3> p1 = *(it);
        it++;
        std::array<double,3> p2 = *(it);
        if (p1[0] >= x_min && p1[0] <= x_max && p1[1] >= y_min && p1[1] <= y_max && p1[2] >= z_min && p1[2] <= z_max) {
            return true;
        }
        if (p2[0] >= x_min && p2[0] <= x_max && p2[1] >= y_min && p2[1] <= y_max && p2[2] >= z_min && p2[2] <= z_max) {
            return true;
        }
        std::array<std::array<double,3>,2> segment_shadow = {{{p1[0],0,p1[2]},{p2[0],0,p2[2]}}};
        set<vector<std::array<double,3>>> circuit_shadow = {{{x_min,0,z_min},{x_max,0,z_min},{x_max,0,z_max},{x_min,0,z_max}}};
        vector<Polyhedron::CircuitIntersection> intersections = Polyhedron::circuit_intersect(segment_shadow, circuit_shadow);
        for (const Polyhedron::CircuitIntersection& intersection : intersections) {
            double y1 = p1[1];
            double y2 = y1+(p2[1]-p1[1])*intersection.alpha;
            if (y2 >= y_min && y2 <= y_max) {
                return true;
            }
        }
        return false;
    }
    static bool inside_polyhedron(Polyhedron poly, std::array<double,3> point) {
        int intersections = 0;
        for (int face_index = 0; face_index < poly.faces.size(); face_index++) {
            set<vector<std::array<double,3>>> circuits = Polyhedron::make_clockwise(poly.circuits(face_index));
            vector<std::array<double,3>>* circuit = Polyhedron::circuit_cut(circuits);
            vector<std::array<std::array<double,3>,3>> triangles = Polyhedron::triangulate(*circuit);
            for (const std::array<std::array<double,3>,3>& triangle : triangles) {
                if (Polyhedron::inside_triangle(triangle, point)) {
                    delete circuit;
                    return true;
                }
            }
            if ((*circuit)[0][2] == (*circuit)[1][2] && (*circuit)[0][2] == (*circuit)[2][2] && (*circuit)[0][2] > point[2]) {
                std::array<double,3> projection = {point[0],point[1],(*circuit)[0][2]};
                bool projection_inside_circuit = false;
                for (const std::array<std::array<double,3>,3>& triangle : triangles) {
                    if (Polyhedron::inside_triangle(triangle, projection)) {
                        intersections++;
                        break;
                    }
                }
            }
            delete circuit;
        }
        return intersections % 2 == 1;
    }
    static set<vector<std::array<double,3>>> del_circuit_helper(set<set<std::array<double,3>>>path_edges, set<std::array<double,3>> start, set<std::array<double,3>> previous, set<std::array<double,3>> current, vector<std::array<double,3>> path, set<set<std::array<double,3>>>& seen) {
        map<std::array<double,3>,set<set<std::array<double,3>>>> edge_lookup;
        for (const set<std::array<double,3>>& edge : path_edges) {
            for (const std::array<double,3>& point : edge) {
                edge_lookup[Polyhedron::round_point(point)].insert(edge);
            }
        }
        set<vector<std::array<double,3>>> output;
        if (current == start) {
            output.insert(path);
            return output;
        }
        if (seen.find(current) != seen.end()) {
            return output;
        }
        seen.insert(current);
        set<std::array<double,3>> difference;
        set_difference(current.begin(),current.end(),previous.begin(),previous.end(),inserter(difference,difference.begin()));
        std::array<double,3> point = *(difference.begin());
        if (find(path.begin(), path.end(), point) != path.end()) {
            return output;
        }
        path.push_back(point);
        previous = current;
        set<set<std::array<double,3>>> temp(edge_lookup[Polyhedron::round_point(point)].begin(),edge_lookup[Polyhedron::round_point(point)].end());
        temp.erase(previous);
        for (const set<std::array<double,3>>& x : temp) {
            set<std::array<double,3>> point_union;
            for (const std::array<double,3>& y : path) {
                point_union.insert(Polyhedron::round_point(y));
            }
            set<std::array<double,3>> edge = Polyhedron::round_edge(x);
            point_union.insert(edge.begin(),edge.end());
            if (Polyhedron::coplanar(point_union)) {
                current = x;
                set<vector<std::array<double,3>>> intermediate = Box::del_circuit_helper(path_edges, start, previous, current, path, seen);
                output.insert(intermediate.begin(), intermediate.end());
            }
        }
        vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
        for (int i = 0; i < output_list.size(); i++) {
            for (int j = 0; j < output_list.size(); j++) {
                set<std::array<double,3>> difference2;
                set_difference(output_list[i].begin(),output_list[i].end(),output_list[j].begin(),output_list[j].end(),inserter(difference2,difference2.begin()));
                if (j > i && output_list[i].size() == output_list[j].size() && !(difference2.size()) && output.find(output_list[j]) != output.end()) {
                    output.erase(output_list[j]);
                }
            }
        }
        return output;
    }
    static set<vector<std::array<double,3>>> del_circuit_helper(set<set<std::array<double,3>>>path_edges) {
        map<std::array<double,3>,set<set<std::array<double,3>>>> edge_lookup;
        for (const set<std::array<double,3>>& edge : path_edges) {
            for (const std::array<double,3>& point : edge) {
                edge_lookup[Polyhedron::round_point(point)].insert(edge);
            }
        }
        set<set<std::array<double,3>>> seen;
        set<vector<std::array<double,3>>> output;
        for (const set<std::array<double,3>>& edge : path_edges) {
            vector<std::array<double,3>> path;
            set<std::array<double,3>> start = edge;
            std::array<double,3> point = *(start.begin());
            path.push_back(point);
            set<set<std::array<double,3>>> temp(edge_lookup[Polyhedron::round_point(point)].begin(),edge_lookup[Polyhedron::round_point(point)].end());
            temp.erase(start);
            for (const set<std::array<double,3>>& x : temp) {
                set<std::array<double,3>> current = x;
                set<vector<std::array<double,3>>> intermediate = Box::del_circuit_helper(path_edges, start, start, current, path, seen);
                output.insert(intermediate.begin(),intermediate.end());
            }
        }
        vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
        for (int i = 0; i < output_list.size(); i++) {
            for (int j = 0; j < output_list.size(); j++) {
                set<std::array<double,3>> difference2;
                set_difference(output_list[i].begin(),output_list[i].end(),output_list[j].begin(),output_list[j].end(),inserter(difference2,difference2.begin()));
                if (j > i && output_list[i].size() == output_list[j].size() && !(difference2.size()) && output.find(output_list[j]) != output.end()) {
                    output.erase(output_list[j]);
                }
            }
        }
        return output;
    }
    class PointDistanceComparator {
      public:
        std::array<double,3> point = {0,0,0};
        PointDistanceComparator(std::array<double,3> point) {
            this->point = point;
        }
        bool operator()(std::array<double,3> a, std::array<double,3> b) const {
            return Polyhedron::distance(point,a) < Polyhedron::distance(point,b);
        }
    };
    struct DelPreprocessOutput {
        set<set<std::array<double,3>>> new_edges;
        set<set<std::array<double,3>>> path_edges;
        int box_index;
    };
    DelPreprocessOutput del_preprocess(Polyhedron poly, int face_index) const {
        set<int> face = poly.faces[face_index];
        set<vector<std::array<double,3>>> circuits = Polyhedron::make_clockwise(poly.circuits(face_index));
        vector<std::array<double,3>>* circuit = Polyhedron::circuit_cut(circuits);
        DelPreprocessOutput output;
        output.box_index = -1;
        bool all = true;
        for (int i = 0; i < circuit->size(); i++) {
            if ((*circuit)[i][0] < x_min || (*circuit)[i][0] > x_max || (*circuit)[i][1] < y_min || (*circuit)[i][1] > y_max || (*circuit)[i][2] < z_min || (*circuit)[i][2] > z_max) {
                all = false;
                break;
            }
        }
        if (all) {
            bool double_break = false;
            vector<set<std::array<double,3>>> edges_grounded;
            for (const set<int>& edge : poly.edges) {
                set<std::array<double,3>> edge_grounded;
                for (const int& index : edge) {
                    edge_grounded.insert(poly.verts[index]);
                }
                edges_grounded.push_back(edge_grounded);
            }
            for (const std::array<double,3>& p1 : (*circuit)) {
                if (find(poly.verts.begin(), poly.verts.end(), p1) != poly.verts.end()) {
                    for (const set<std::array<double,3>>& edge : edges_grounded) {
                        int intersection = 0;
                        std::array<double,3> p2;
                        if (edge.find(p1) != edge.end()) {
                            p2 = *(edge.begin());
                            if (p2 == p1) {
                                p2 = *(++edge.begin());
                            }
                        }
                        if (find(circuit->begin(), circuit->end(), p2) == circuit->end()) {
                            int dim;
                            for (int i = 0; i < 3; i++) {
                                if (p1[i] != p2[i]) {
                                    dim = i;
                                    break;
                                }
                            }
                            std::array<double,3> vec = {p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
                            std::array<double,2> min_max;
                            switch (dim) {
                                case 0:
                                    min_max = {x_min, x_max};
                                    break;
                                case 1:
                                    min_max = {y_min, y_max};
                                    break;
                                case 2:
                                    min_max = {z_min, z_max};
                                    break;
                            }
                            for (const double& w : min_max) {
                                std::array<std::array<std::array<double,3>,3>,2> triangles;
                                switch (dim) {
                                    case 0:
                                        triangles = {{std::array<std::array<double, 3>, 3>{{std::array<double, 3>{w,y_min,z_min},std::array<double, 3>{w,y_max,z_min},std::array<double, 3>{w,y_max,z_max}}},std::array<std::array<double, 3>, 3>{{std::array<double, 3>{w,y_min,z_min},std::array<double, 3>{w,y_min,z_max},std::array<double, 3>{w,y_max,z_max} }}}};
                                        break;
                                    case 1:
                                        triangles = {{std::array<std::array<double, 3>, 3>{{std::array<double, 3>{x_min,w,z_min},std::array<double, 3>{x_max,w,z_min},std::array<double, 3>{x_max,w,z_max}}},std::array<std::array<double, 3>, 3>{{std::array<double, 3>{x_min,w,z_min},std::array<double, 3>{x_min,w,z_max},std::array<double, 3>{x_max,w,z_max} }}}};
                                        break;
                                    case 2:
                                        triangles = {{std::array<std::array<double, 3>, 3>{{std::array<double, 3>{x_min,y_min,w},std::array<double, 3>{x_max,y_min,w},std::array<double, 3>{x_max,y_max,w}}},std::array<std::array<double, 3>, 3>{{std::array<double, 3>{x_min,y_min,w},std::array<double, 3>{x_min,y_max,w},std::array<double, 3>{x_max,y_max,w} }}}};
                                        break;
                                }
                                if ((vec[dim] > 0 && w > p1[dim]) || (vec[dim] < 0 && w < p1[dim])) {
                                    std::array<double,3> projection = {p1[0],p1[2],p1[3]};
                                    projection[dim] = w;
                                    for (const std::array<std::array<double,3>,3>& triangle : triangles) {
                                        if (Polyhedron::inside_triangle(triangle,projection)) {
                                            intersection++;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if (intersection % 2 == 0) {
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
                return output; 
            }
        }
        vector<std::array<double,3>> box;
        if ((*circuit)[0][0] == (*circuit)[1][0] && (*circuit)[0][0] == (*circuit)[2][0] && (*circuit)[0][0] >= x_min && (*circuit)[0][0] <= x_max) {
            box = {{(*circuit)[0][0],y_min,z_min},{(*circuit)[0][0],y_max,z_min},{(*circuit)[0][0],y_max,z_max},{(*circuit)[0][0],y_min,z_max}};
        } else if ((*circuit)[0][1] == (*circuit)[1][1] && (*circuit)[0][1] == (*circuit)[2][1] && (*circuit)[0][1] >= y_min && (*circuit)[0][1] <= y_max) {
            box = {{x_min,(*circuit)[0][1],z_min},{x_max,(*circuit)[0][1],z_min},{x_max,(*circuit)[0][1],z_max},{x_min,(*circuit)[0][1],z_max}};
        } else if ((*circuit)[0][2] == (*circuit)[1][2] && (*circuit)[0][2] == (*circuit)[2][2] && (*circuit)[0][2] >= z_min && (*circuit)[0][2] <= z_max) {
            box = {{x_min,y_min,(*circuit)[0][2]},{x_max,y_min,(*circuit)[0][2]},{x_max,y_max,(*circuit)[0][2]},{x_min,y_max,(*circuit)[0][2]}};
        }
        vector<std::array<double,3>> intersections;
        vector<bool*> intersections_inside;
        for (const int& edge_index : poly.faces[face_index]) {
            for (const int& index : poly.edges[edge_index]) {
                cout << "[" << poly.verts[index][0] << "," << poly.verts[index][1] << "," << poly.verts[index][2] << "] ";
            }
            cout << endl;
        }
        for (int i = 0; i < circuit->size(); i++) {
            bool any = false;
            std::array<double,3> p1 = (*circuit)[(i-1+circuit->size())%circuit->size()];
            std::array<double,3> p2 = (*circuit)[i];
            for (const vector<std::array<double,3>>& y : poly.circuits(face_index)) {
                if (find(y.begin(),y.end(),p1) != y.end() && find(y.begin(),y.end(), p2) != y.end()) {
                    any = true;
                    break;
                }
            }
            if (!any) {
                continue;
            }
            vector<std::array<double,3>> last_intersections;
            for (int j = 0; j < box.size(); j++) {
                if (!Polyhedron::colinear((set<std::array<double,3>>){p1,p2,box[(j-1+box.size())%box.size()],box[j]})) {
                    std::array<double,3>* intersection = Box::intersect_segments({p1, p2}, {box[(j-1+box.size())%box.size()], box[j]});
                    if (intersection != NULL) {
                        last_intersections.push_back(Polyhedron::round_point(*intersection));
                    }
                    delete intersection;
                } else {
                    if (Polyhedron::point_on_segment({box[(j-1+box.size())%box.size()], box[j]}, p1) && Polyhedron::point_on_segment({box[(j-1+box.size())%box.size()], box[j]}, p2)) {
                        if (find(box.begin(),box.end(),p1) == box.end()) {
                            last_intersections.push_back(p1);
                        }
                        if (find(box.begin(),box.end(),p2) == box.end()) {
                            last_intersections.push_back(p2);
                        }
                    }
                }
            }
            cout << "p1 [" << p1[0] << "," << p1[1] << "," << p1[2] << "] ";
            cout << "p2 [" << p2[0] << "," << p2[1] << "," << p2[2] << "] " << endl;
            cout << "last intersections: ";
            for (const std::array<double,3>& intersection : last_intersections) {
                cout << "[" << intersection[0] << "," << intersection[1] << "," << intersection[2] << "] ";
            }
            cout << endl; 
            switch (last_intersections.size()) {
                case 0:
                    if (!this->is_inside(p1) && !this->is_inside(p2)) {
                        output.new_edges.insert(Polyhedron::round_edge({p1,p2}));
                    }
                    break;
                case 1:
                    if (x_min > p2[0] || p2[0] > x_max || y_min > p2[1] || p2[1] > y_max || z_min > p2[2] || p2[2] > z_max) {
                        set<std::array<double,3>> edge = {p2,last_intersections[0]};
                        if (edge.size() == 2) {
                            output.new_edges.insert(Polyhedron::round_edge(edge));
                        }
                        vector<std::array<double,3>>::iterator it = find(intersections.begin(),intersections.end(), last_intersections[0]);
                        if (it != intersections.end()) {
                            int index = std::distance(intersections.begin(),it);
                            if (intersections_inside[index] != NULL && *(intersections_inside[index]) == true) {
                                intersections.erase(intersections.begin()+index);
                                delete intersections_inside[index];
                                intersections_inside.erase(intersections_inside.begin()+index);
                            } else if (intersections_inside[index] == NULL){
                                intersections_inside[index] = new bool(false);
                            }
                        } else {
                            intersections.push_back(last_intersections[0]);
                            intersections_inside.push_back(new bool(false));
                        }
                    } else if (x_min > p1[0] || p1[0] > x_max || y_min > p1[1] || p1[1] > y_max || z_min > p1[2] || p1[2] > z_max) {
                        set<std::array<double,3>> edge = {p1,last_intersections[0]};
                        if (edge.size() == 2) {
                            output.new_edges.insert(Polyhedron::round_edge(edge));
                        }
                        vector<std::array<double,3>>::iterator it = find(intersections.begin(),intersections.end(), last_intersections[0]);
                        if (it != intersections.end()) {
                            int index = std::distance(intersections.begin(),it);
                            if (intersections_inside[index] != NULL && *(intersections_inside[index]) == false) {
                                intersections.erase(intersections.begin()+index);
                                delete intersections_inside[index];
                                intersections_inside.erase(intersections_inside.begin()+index);
                            } else if (intersections_inside[index] == NULL) {
                                intersections_inside[index] = new bool(true);
                            }
                        } else {
                            intersections.push_back(last_intersections[0]);
                            intersections_inside.push_back(new bool(true));
                        }
                    }
                    break;
                case 2:
                    std::array<bool*,2> last_intersections_inside;
                    PointDistanceComparator comp(p1);
                    sort(last_intersections.begin(),last_intersections.end(),comp);
                    if ((x_min > p1[0] || p1[0] > x_max || y_min > p1[1] || p1[1] > y_max || z_min > p1[2] || p1[2] > z_max) && (x_min > p2[0] || p2[0] > x_max || y_min > p2[1] || p2[1] > y_max || z_min > p2[2] || p2[2] > z_max)) {
                        last_intersections_inside[1] = new bool(false);
                        last_intersections_inside[0] = new bool(true);
                        set<std::array<double,3>> edge = {p1,last_intersections[0]};
                        if (edge.size() == 2) {
                            output.new_edges.insert(Polyhedron::round_edge(edge));
                        }
                        edge = {p2,last_intersections[1]};
                        if (edge.size() == 2) {
                            output.new_edges.insert(Polyhedron::round_edge(edge));
                        }
                    } else if (x_min > p1[0] || p1[0] > x_max || y_min > p1[1] || p1[1] > y_max || z_min > p1[2] || p1[2] > z_max) {
                        last_intersections_inside[1] = NULL;
                        last_intersections_inside[0] = new bool(true);
                        set<std::array<double,3>> edge = {p1,last_intersections[0]};
                        if (edge.size() == 2) {
                            output.new_edges.insert(Polyhedron::round_edge(edge));
                        }
                    } else if (x_min > p2[0] || p2[0] > x_max || y_min > p2[1] || p2[1] > y_max || z_min > p2[2] || p2[2] > z_max) {
                        last_intersections_inside[1] = new bool(false);
                        last_intersections_inside[0] = NULL;
                        set<std::array<double,3>> edge = {p2,last_intersections[1]};
                        if (edge.size() == 2) {
                            output.new_edges.insert(Polyhedron::round_edge(edge));
                        } 
                    } else {
                        last_intersections_inside[1] = NULL;
                        last_intersections_inside[0] = NULL;
                    }
                    for (int j = 0; j < last_intersections.size(); j++) {
                        vector<std::array<double,3>>::iterator it = find(intersections.begin(),intersections.end(), last_intersections[j]);
                        if (it != intersections.end()) {
                            int index = std::distance(intersections.begin(),it);
                            if (last_intersections_inside[j] != NULL) {
                                if (intersections_inside[index] != NULL && *(intersections_inside[index]) != *(last_intersections_inside[j])) {
                                    intersections.erase(intersections.begin()+index);
                                    delete intersections_inside[index];
                                    intersections_inside.erase(intersections_inside.begin()+index);
                                } else if (intersections_inside[index] == NULL) {
                                    intersections_inside[index] = last_intersections_inside[j];
                                }
                            }
                        } else {
                            intersections.push_back(last_intersections[j]);
                            intersections_inside.push_back(last_intersections_inside[j]);
                        }
                    }
                    //for (bool* inside : last_intersections_inside) {
                    //    delete inside;
                    //}
            }
        }
        set<set<std::array<double,3>>> f;
        for (const int& edge_index : face) {
            set<std::array<double,3>> edge;
            for (const int& index : poly.edges[edge_index]) {
                edge.insert(poly.verts[index]);
            }
            f.insert(edge);
        }
        for (int box_index = 0; box_index < box.size(); box_index++) {
            set<std::array<double,3>> edge1 = {box[(box_index-1+box.size())%box.size()],box[box_index]};
            for (const set<std::array<double,3>>& edge2 : f) {
                bool no_break = true;
                for (std::array<double,3> point : edge1) {
                    if (!Polyhedron::point_on_segment(Polyhedron::round_edge(edge2), Polyhedron::round_point(point))) {
                        no_break = false;
                        break;
                    }
                }
                if (no_break) {
                    set<std::array<double,3>>::iterator it = edge2.begin();
                    std::array<double,3> p1 = *(it);
                    std::array<double,3> p2 = *(++it);
                    std::array<std::array<double,3>,2> edge1_array;
                    it = edge1.begin();
                    edge1_array[0] = *(it);
                    edge1_array[1] = *(++it);
                    PointDistanceComparator comp = PointDistanceComparator(p1);
                    sort(edge1_array.begin(), edge1_array.end(), comp);
                    set<std::array<double,3>> edge3 = {p1, edge1_array[0]};
                    if (Polyhedron::distance(p1, edge1_array[0])) {
                        output.new_edges.insert(Polyhedron::round_edge(edge3));
                    }
                    comp = PointDistanceComparator(p2);
                    sort(edge1_array.begin(), edge1_array.end(), comp);
                    edge3 = {p2, edge1_array[0]};
                    if (Polyhedron::distance(p2, edge1_array[0])) {
                        output.new_edges.insert(Polyhedron::round_edge(edge3));
                    }
                    break;
                }
            }
        }
        std::array<std::array<double,3>,4> box_old;
        for (int i = 0; i < box.size(); i++) {
            box_old[i] = box[i];
        }
        if (intersections.size() > 1) {
            if ((*circuit)[0][0] == (*circuit)[1][0] && (*circuit)[0][0] == (*circuit)[2][0]) {
                if (Polyhedron::round_float((*circuit)[0][0]) == Polyhedron::round_float(x_min)) {
                    output.box_index = 0;
                }
                if (Polyhedron::round_float((*circuit)[0][0]) == Polyhedron::round_float(x_max)) {
                    output.box_index = 1;
                }
            } else if ((*circuit)[0][1] == (*circuit)[1][1] && (*circuit)[0][1] == (*circuit)[2][1]) {
                if (Polyhedron::round_float((*circuit)[0][1]) == Polyhedron::round_float(y_min)) {
                    output.box_index = 2;
                }
                if (Polyhedron::round_float((*circuit)[0][1]) == Polyhedron::round_float(y_max)) {
                    output.box_index = 3;
                }
            } else if ((*circuit)[0][2] == (*circuit)[1][2] && (*circuit)[0][2] == (*circuit)[2][2]) {
                if (Polyhedron::round_float((*circuit)[0][2]) == Polyhedron::round_float(z_min)) {
                    output.box_index = 4;
                }
                if (Polyhedron::round_float((*circuit)[0][2]) == Polyhedron::round_float(z_max)) {
                    output.box_index = 5;
                }
            }
            vector<std::array<double,3>> intersections_rounded;
            for (const std::array<double,3>& intersection : intersections) {
                intersections_rounded.push_back(Polyhedron::round_point(intersection));
                for (int j = 0; j < box.size(); j++) {
                    if (Polyhedron::point_on_segment(Polyhedron::round_edge({box[(j-1+box.size())%box.size()],box[j]}), Polyhedron::round_point(intersection)) && Polyhedron::round_point(box[(j-1+box.size())%box.size()]) != Polyhedron::round_point(intersection) && Polyhedron::round_point(box[j]) != Polyhedron::round_point(intersection)) {
                        box.insert(box.begin() + j, intersection);
                        break;
                    }
                }
            }
            cout << "box: ";
            for (int i = 0; i < box.size(); i++) {
                cout << "[" << box[i][0] << "," << box[i][1] << "," << box[i][2] << "] ";
            }
            cout << endl;
            cout << "intersections: ";
            for (int i = 0; i < intersections.size(); i++) {
                cout << "[" << intersections[i][0] << "," << intersections[i][1] << "," << intersections[i][2] << "] ";
            }
            cout << endl;
            vector<set<set<std::array<double,3>>>> edge_components;
            for (const set<std::array<double,3>>& edge : output.new_edges) {
                edge_components.push_back({Polyhedron::round_edge(edge)});
            }
            while (true) {
                bool double_break = false;
                for (int i = 0; i < edge_components.size(); i++) {
                    for (int j = i+1; j < edge_components.size(); j++) {
                        bool any = false;
                        for (const set<std::array<double,3>>& edge1 : edge_components[i]) {
                            for (const set<std::array<double,3>>& edge2 : edge_components[j]) {
                                set<std::array<double,3>> edge_intersection;
                                set_intersection(edge1.begin(),edge1.end(),edge2.begin(),edge2.end(),inserter(edge_intersection,edge_intersection.begin()));
                                if (edge_intersection.size()) {
                                    any = true;
                                    break;
                                }
                            }
                            if (any) {
                                break;
                            }
                        }
                        if (any) {
                            edge_components[i].insert(edge_components[j].begin(),edge_components[j].end());
                            edge_components.erase(edge_components.begin()+j);
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
            vector<set<std::array<double,3>>> point_components;
            for (const set<set<std::array<double,3>>>& component : edge_components) {
                set<std::array<double,3>> point_component;
                for (const set<std::array<double,3>>& edge : component) {
                    for (const std::array<double,3>& point : edge) {
                        point_component.insert(point);
                    }
                }
                point_components.push_back(point_component);
            }
            map<std::array<double,3>,set<std::array<double,3>>> component_lookup;
            for (const set<std::array<double,3>>& component : point_components) {
                for (const std::array<double,3>& point : component) {
                    component_lookup[point] = component;
                }
            }
            for (int i = 0; i < intersections.size(); i++) {
                vector<std::array<double,3>>::iterator it = find(box.begin(),box.end(),intersections[i]);
                int index;
                vector<std::array<double,3>> path1;
                if (it != box.end()) {
                    index = std::distance(box.begin(),it);
                    path1.push_back(intersections[i]);
                } else {
                    it = find(box.begin(),box.end(),Polyhedron::round_point(intersections[i]));
                    index = std::distance(box.begin(),it);
                    path1.push_back(Polyhedron::round_point(intersections[i]));
                }
                bool* inside1 = intersections_inside[i];
                bool* inside2;
                index++;
                index %= box.size();
                while (true) {
                    path1.push_back(box[index]);
                    it = find(intersections_rounded.begin(),intersections_rounded.end(), Polyhedron::round_point(box[index]));
                    if (it != intersections_rounded.end()) {
                        inside2 = intersections_inside[std::distance(intersections_rounded.begin(),it)];
                        break;
                    }
                    index++;
                    index %= box.size();
                }
                it = find(box.begin(),box.end(),intersections[i]);
                vector<std::array<double,3>> path2;
                if (it != box.end()) {
                    index = std::distance(box.begin(),it);
                    path2.push_back(intersections[i]);
                } else {
                    it = find(box.begin(),box.end(),Polyhedron::round_point(intersections[i]));
                    index = std::distance(box.begin(),it);
                    path2.push_back(Polyhedron::round_point(intersections[i]));
                }
                bool* inside3 = intersections_inside[i];
                bool* inside4;
                index--;
                index += box.size();
                index %= box.size();
                while (true) {
                    path2.push_back(box[index]);
                    it = find(intersections_rounded.begin(),intersections_rounded.end(), Polyhedron::round_point(box[index]));
                    if (it != intersections_rounded.end()) {
                        inside4 = intersections_inside[std::distance(intersections_rounded.begin(),it)];
                        break;
                    }
                    index--;
                    index += box.size();
                    index %= box.size();
                }
                cout << (inside1 != NULL) << (inside2 != NULL) << endl;
                cout << (inside3 != NULL) << (inside4 != NULL) << endl;
                bool is_border_path1 = inside1 != NULL && inside2 != NULL && !(*inside1) && (*inside2);
                bool is_border_path2 = inside3 != NULL && inside4 != NULL && !(*inside3) && (*inside4);
                set<set<std::array<double,3>>> components;
                for (const std::array<double,3>& x : path1) {
                    if (component_lookup.find(Polyhedron::round_point(x)) != component_lookup.end()) {
                        components.insert(component_lookup[Polyhedron::round_point(x)]);
                    }
                }
                if (components.size() > 1) {
                    is_border_path1 = false;
                }
                components.clear();
                for (const std::array<double,3>& x : path2) {
                    if (component_lookup.find(Polyhedron::round_point(x)) != component_lookup.end()) {
                        components.insert(component_lookup[Polyhedron::round_point(x)]);
                    }
                }
                if (components.size() > 1) {
                    is_border_path2 = false;
                }
                if (path1[path1.size()-1] == path2[path2.size()-1] && is_border_path1 && is_border_path2) {
                    double distance1 = 0;
                    for (int j = 1; j < path1.size(); j++) {
                        distance1 += Polyhedron::distance(path1[j-1],path1[j]);
                    }
                    double distance2 = 0;
                    for (int j = 1; j < path2.size(); j++) {
                        distance2 += Polyhedron::distance(path2[j-1],path2[j]);
                    }
                    if (distance2 > distance1) {
                        is_border_path1 = false;
                    } else {
                        is_border_path2 = false;
                    }
                }
                std::array<vector<std::array<double,3>>,2> paths = {path1,path2};
                std::array<bool,2> is_border_path = {is_border_path1, is_border_path2};
                for (int path_index = 0; path_index < 2; path_index++) {
                    cout << "path " << path_index << ": ";
                    for (int k = 0; k < paths[path_index].size(); k++) {
                        cout << "[" << paths[path_index][k][0] << "," << paths[path_index][k][1] << "," << paths[path_index][k][2] << "] ";
                    }
                    cout << is_border_path[path_index] << endl;
                    if (is_border_path[path_index]) {
                        bool no_break = true;
                        for (int k = 0; k < paths[path_index].size()-1; k++) {
                            std::array<double,3> point;
                            for (int j = 0; j < 3; j++) {
                                point[j] = (paths[path_index][k][j]+paths[path_index][k+1][j])/2;
                            }
                            if (!Box::inside_polyhedron(poly,point)) {
                                no_break = false;
                                break;
                            }
                        }
                        if (no_break) {
                            for (int j = 0; j < paths[path_index].size()-1; j++) {
                                set<std::array<double,3>> edge1 = {paths[path_index][j],paths[path_index][j+1]};
                                std::array<std::array<double,3>,2> edge1_array = {paths[path_index][j],paths[path_index][j+1]};
                                no_break = true;
                                for (const set<std::array<double,3>>& edge2 : output.new_edges) {
                                    bool no_break2 = true;
                                    for (const std::array<double,3> point : edge1) {
                                        if (!Polyhedron::point_on_segment(Polyhedron::round_edge(edge2), Polyhedron::round_point(point))) {
                                            no_break2 = false;
                                            break;
                                        }
                                    }
                                    if (no_break2) {
                                        set<std::array<double,3>>::iterator it2 = edge2.begin();
                                        std::array<double,3> p1 = *(it2);
                                        std::array<double,3> p2 = *(++it2);
                                        PointDistanceComparator comp = PointDistanceComparator(p1);
                                        sort(edge1_array.begin(),edge1_array.end(),comp);
                                        if (Polyhedron::distance(p1,edge1_array[0])) {
                                            output.new_edges.insert(Polyhedron::round_edge({p1,edge1_array[0]}));
                                        }
                                        comp = PointDistanceComparator(p2);
                                        sort(edge1_array.begin(),edge1_array.end(),comp);
                                        if (Polyhedron::distance(p2,edge1_array[0])) {
                                            output.new_edges.insert(Polyhedron::round_edge({p2,edge1_array[0]}));
                                        }
                                        no_break = false;
                                        break;
                                    }
                                }
                                if (no_break) {
                                    if (edge1.size() == 2) {
                                        output.new_edges.insert(Polyhedron::round_edge(edge1));
                                        output.path_edges.insert(Polyhedron::round_edge(edge1));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            bool all = true;
            for (const std::array<double,3>& box_point : box_old) {
                vector<std::array<double,3>>::iterator it = find(intersections.begin(),intersections.end(),box_point);
                if (it == intersections.end() || intersections_inside[std::distance(intersections.begin(),it)] != NULL) {
                    all = false;
                    break;
                }
            }
            set<std::array<double,3>> difference;
            set_difference(circuit->begin(), circuit->end(), box_old.begin(), box_old.end(), inserter(difference, difference.begin()));
            if (all && !(circuit->size() == 4 && !difference.size())) {
                for (int k = 0; k < intersections.size(); k++) {
                    set<std::array<double,3>> edge1 = {intersections[(k-1+intersections.size())%intersections.size()],intersections[k]};
                    std::array<std::array<double,3>,2> edge1_array = {intersections[(k-1+intersections.size())%intersections.size()],intersections[k]};
                    std::array<double,3> p1 = edge1_array[0];
                    std::array<double,3> p2 = edge1_array[1];
                    set<set<std::array<double,3>>> edge_union(output.new_edges.begin(),output.new_edges.end());
                    for (const int& edge_index : face) {
                        set<std::array<double,3>> edge;
                        for (const int& index : poly.edges[edge_index]) {
                            edge.insert(poly.verts[index]);
                        }
                        edge_union.insert(edge);
                    }
                    bool any = false;
                    std::array<double,3> p3 = {(p1[0]+p2[0])/2,(p1[1]+p2[1])/2,(p1[2]+p2[2])/2};
                    for (set<std::array<double,3>> edge2 : edge_union) {
                        if (Polyhedron::point_on_segment(Polyhedron::round_edge(edge2), Polyhedron::round_point(p3))) {
                            any = true;
                            break;
                        } 
                    }
                    if (!any) {
                        continue;
                    }
                    bool no_break = true;
                    for (const set<std::array<double,3>>& edge2 : edge_union) {
                        bool no_break2 = true;
                        for(const std::array<double,3>& point : edge1) {
                            if (!Polyhedron::point_on_segment(Polyhedron::round_edge(edge2), Polyhedron::round_point(point))) {
                                no_break2 = false;
                                break;
                            }
                        }
                        if (no_break2) {
                            if (output.new_edges.find(edge2) != output.new_edges.end()) {
                                set<std::array<double,3>>::iterator it = edge2.begin();
                                p1 = *(it);
                                p2 = *(++it);
                                output.new_edges.erase(Polyhedron::round_edge(edge2));
                                PointDistanceComparator comp = PointDistanceComparator(p1);
                                sort(edge1_array.begin(),edge1_array.end(),comp);
                                if (Polyhedron::distance(p1,edge1_array[0])) {
                                    output.new_edges.insert({p1,edge1_array[0]});
                                }
                                comp = PointDistanceComparator(p2);
                                sort(edge1_array.begin(),edge1_array.end(),comp);
                                if (Polyhedron::distance(p2,edge1_array[0])) {
                                    output.new_edges.insert({p2,edge1_array[0]});
                                }
                            } else {
                                output.new_edges.insert(Polyhedron::round_edge(edge1));
                                output.path_edges.insert(Polyhedron::round_edge(edge1));
                            }
                            no_break = false;
                            break;
                        }
                    }
                    if (no_break) {
                        output.new_edges.insert(Polyhedron::round_edge(edge1));
                        output.path_edges.insert(Polyhedron::round_edge(edge1));
                    }
                }
            }
        } else if (!intersections.size()) {
            bool any = false;
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*circuit)) {
                if (box.size() && Polyhedron::inside_triangle(triangle,box[0])) {
                    any = true;
                    break;
                }
            }
            if (any) {
                for (int j = 0; j < box.size(); j++) {
                    output.new_edges.insert(Polyhedron::round_edge({box[(j-1+box.size())%box.size()],box[j]}));
                }
                if ((*circuit)[0][0] == (*circuit)[1][0] && (*circuit)[0][0] == (*circuit)[2][0]) {
                    if (Polyhedron::round_float((*circuit)[0][0]) == Polyhedron::round_float(x_min)) {
                        output.box_index = 0;
                    }
                    if (Polyhedron::round_float((*circuit)[0][0]) == Polyhedron::round_float(x_max)) {
                        output.box_index = 1;
                    }
                } else if ((*circuit)[0][1] == (*circuit)[1][1] && (*circuit)[0][1] == (*circuit)[2][1]) {
                    if (Polyhedron::round_float((*circuit)[0][1]) == Polyhedron::round_float(y_min)) {
                        output.box_index = 2;
                    }
                    if (Polyhedron::round_float((*circuit)[0][1]) == Polyhedron::round_float(y_max)) {
                        output.box_index = 3;
                    }
                } else if ((*circuit)[0][2] == (*circuit)[1][2] && (*circuit)[0][2] == (*circuit)[2][2]) {
                    if (Polyhedron::round_float((*circuit)[0][2]) == Polyhedron::round_float(z_min)) {
                        output.box_index = 4;
                    }
                    if (Polyhedron::round_float((*circuit)[0][2]) == Polyhedron::round_float(z_max)) {
                        output.box_index = 5;
                    }
                }
            }
        }
        
        vector<std::array<double,3>>* new_circuit = Polyhedron::circuit_cut(Polyhedron::make_clockwise(del_circuit_helper(output.new_edges)));
        std::array<double,3> center = {0,0,0};
        for (const std::array<double,3>& x : box_old) {
            center[0] += x[0]/box_old.size();
            center[1] += x[1]/box_old.size();
            center[2] += x[2]/box_old.size();
        }
        center = Polyhedron::round_point(center);
        if (center[0] != x_min && center[0] != x_max && center[1] != y_min && center[1] != y_max && center[2] != z_min && center[2] != z_max) {
            bool any = false;
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*new_circuit)) {
                if (Polyhedron::inside_triangle(triangle, center)) {
                    any = true;
                }
            }
            if (any) {
                cout << "howdy " << center[0] << "," << center[1] << "," << center[2] <<endl;
                for (int i = 0; i < box_old.size(); i++) {
                    bool no_break = true;
                    std::array<std::array<double,3>,2> edge1_array = {Polyhedron::round_point(box_old[(i-1+box_old.size())%box_old.size()]), Polyhedron::round_point(box_old[i])};
                    for (const set<std::array<double,3>>& edge2 : output.new_edges) {
                        if (!Polyhedron::point_on_segment(Polyhedron::round_edge(edge2), Polyhedron::round_point(box_old[(i-1+box_old.size())%box_old.size()]))) {
                            continue;
                        } 
                        if (!Polyhedron::point_on_segment(Polyhedron::round_edge(edge2), Polyhedron::round_point(box_old[i]))) {
                            continue;
                        }
                        set<std::array<double,3>>::iterator it = edge2.begin();
                        std::array<double,3> p1 = *(it);
                        std::array<double,3> p2 = *(++it);
                        PointDistanceComparator comp = PointDistanceComparator(p1);
                        sort(edge1_array.begin(),edge1_array.end(),comp);
                        if (Polyhedron::distance(p1,edge1_array[0])) {
                            output.new_edges.insert(Polyhedron::round_edge({p1,edge1_array[0]}));
                        }
                        comp = PointDistanceComparator(p2);
                        sort(edge1_array.begin(),edge1_array.end(),comp);
                        if (Polyhedron::distance(p2,edge1_array[0])) {
                            output.new_edges.insert(Polyhedron::round_edge({p2,edge1_array[0]}));
                        }
                        output.new_edges.erase(edge2);
                        no_break = false;
                        break;
                    }
                    if (no_break) {
                        output.new_edges.insert(Polyhedron::round_edge({box_old[(i-1+box_old.size())%box_old.size()],box_old[i]}));   
                    }
                }
            }
        }
        for (bool* inside : intersections_inside) {
            delete inside;
        }
        delete circuit;
        return output;
    }
    Polyhedron del(Polyhedron poly) const {
        std::array<double,2> x_min_max = {x_min, x_max};
        std::array<double,2> y_min_max = {y_min, y_max};
        std::array<double,2> z_min_max = {z_min, z_max};
        vector<set<std::array<double,3>>> box_faces;
        for (const double& x : x_min_max) {
            set<std::array<double,3>> box_face;
            for (const double& y : y_min_max) {
                for (const double& z : z_min_max) {
                    box_face.insert({x,y,z});
                }
            }
            box_faces.push_back(box_face);
        }
        for (const double& y : y_min_max) {
            set<std::array<double,3>> box_face;
            for (const double& x : x_min_max) {
                for (const double& z : z_min_max) {
                    box_face.insert({x,y,z});
                }
            }
            box_faces.push_back(box_face);
        }
        for (const double& z : z_min_max) {
            set<std::array<double,3>> box_face;
            for (const double& x : x_min_max) {
                for (const double& y : y_min_max) {
                    box_face.insert({x,y,z});
                }
            }
            box_faces.push_back(box_face);
        }
        std::array<double,3> center = {(x_min+x_max)/2,(y_min+y_max)/2,(z_min+z_max)/2};
        cout << (Box::inside_polyhedron(poly, center)) << endl;
        if (!Box::inside_polyhedron(poly, center)) {
            bool triple_break = false;
            for (const set<int>& edge : poly.edges) {
                set<std::array<double,3>> edge_grounded;
                for (const int& index : edge) {
                    edge_grounded.insert(poly.verts[index]);
                }
                if (this->intersect(edge_grounded)) {
                    bool any1 = false;
                    for (const double& x1 : x_min_max) {
                        bool all = true;
                        for (const double& y1 : y_min_max) {
                            for (const double& z1 : z_min_max) {
                                for (const double& x2 : x_min_max) {
                                    for (const double& y2 : y_min_max) {
                                        for (const double& z2 : z_min_max) {
                                            if ((int)(x1 != x2)+(int)(y1 != y2)+(int)(z1 != z2) == 1) {
                                                for (const std::array<double,3>& point : edge_grounded) {
                                                    set<std::array<double,3>> temp_set = {{x1,y1,z1},{x2,y2,z2},point};
                                                    if (!Polyhedron::colinear(temp_set)) {
                                                        all = false;                    
                                                        break;
                                                    }
                                                }
                                                if(all) {
                                                    any1 = true;
                                                }
                                            }
                                            if (all) {
                                                break;
                                            }
                                        }
                                        if (all) {
                                            break;
                                        }
                                    }
                                    if (all) {
                                        break;
                                    }
                                }
                                if (all) {
                                    break;
                                }
                            }
                            if (all) {
                                break;
                            }
                        }
                        if (all) {
                            break;
                        }
                    }
                    if (!any1) {
                        bool any2 = false;
                        for (const set<std::array<double,3>>& box_face : box_faces) {
                            set<std::array<double,3>> temp_set(edge_grounded.begin(),edge_grounded.end());
                            temp_set.insert(box_face.begin(),box_face.end());
                            if (Polyhedron::coplanar(temp_set)) {
                                any2 = true;
                                break;
                            }
                        }
                        if (!any2) {
                            for (const std::array<double,3>& p1 : edge_grounded) {
                                if (x_min <= p1[0] && p1[0] <= x_max && y_min <= p1[1] && p1[1] <= y_max && z_min <= p1[2] && p1[2] <= z_max) {
                                    set<std::array<double,3>> difference(edge_grounded.begin(),edge_grounded.end());
                                    difference.erase(p1);
                                    std::array<double,3> p2 = *(difference.begin());
                                    std::array<double,3> vec;
                                    for (int i = 0; i < 3; i++) {
                                        vec[i] = p2[i]-p1[i];
                                    }
                                    for (int i = 0; i < 20; i++) {
                                        if (vec[0]/pow(2,i) == 0 && vec[1]/pow(2,i) == 0 && vec[2]/pow(2,i) == 0) {
                                            break;
                                        }
                                        std::array<double,3> p3;
                                        for (int j = 0; j < 3; j++) {
                                            p3[j] = p1[j]+vec[j]/pow(2,i); 
                                        }
                                        if (x_min <= p3[0] && p3[0] <= x_max && y_min <= p3[1] && p3[1] <= y_max && z_min <= p3[2] && p3[2] <= z_max) {
                                            triple_break = true;
                                            break;
                                        }
                                    }
                                }
                                if (triple_break) {
                                    break;
                                }
                            }
                        }
                    }
                    if (triple_break) {
                        break;
                    }    
                }
            }
            if (!triple_break) {
                return poly;
            }
        }
        vector<set<set<std::array<double,3>>>> faces;
        set<set<std::array<double,3>>> path_edges;
        std::array<int,6> box_map = {-1,-1,-1,-1,-1,-1};
        for (int i = 0; i < poly.faces.size(); i++) {
            DelPreprocessOutput output = this->del_preprocess(poly, i);
            faces.push_back(output.new_edges);
            path_edges.insert(output.path_edges.begin(),output.path_edges.end());
            if (output.box_index > -1) {
                box_map[output.box_index] = i;
            }
            cout << "new_edges " << output.new_edges.size() << endl;
            for (const set<std::array<double,3>>& edge : output.new_edges) {
                for (const std::array<double,3>& point : edge) {
                    cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
                }
                cout << endl;
            }
            /*
            cout << endl;
            for (const set<std::array<double,3>>& edge1 : output.new_edges) {
                for (const set<std::array<double,3>>& edge2 : output.new_edges) {
                    for (const std::array<double,3>& point : edge1) {
                        cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
                    }
                    cout << ";";
                    for (const std::array<double,3>& point : edge2) {
                        cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
                    }
                    cout << (edge1 == edge2) << endl;
                }
            }
            cout << endl;
            set<std::array<double,3>> point_set;
            for (const set<std::array<double,3>>& edge : output.new_edges) {
                for (const std::array<double,3>& point : edge) {
                    point_set.insert(point);
                }
            } 
            for (const std::array<double,3>& point : point_set) {
                cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
            }
            cout << endl;
            for (const std::array<double,3>& p1 : point_set) {
                for (const std::array<double,3>& p2 : point_set) {
                    cout << "[" << p1[0] << "," << p1[1] << "," << p1[2] << "] ";
                    cout << "[" << p2[0] << "," << p2[1] << "," << p2[2] << "] ";
                    cout << (p1 == p2) << " " << (Polyhedron::round_point(p1) == Polyhedron::round_point(p2)) << " " << Polyhedron::distance(p1,p2) << endl;
                }
            }*/
        }
        
        vector<int> old_face_indices;
        for (int i = 0; i < faces.size(); i++) {
            if (faces[i].size()) {
                old_face_indices.push_back(i);
            }
        }
        vector<set<int>> old_faces;
        for (const int& face_index : old_face_indices) {
            old_faces.push_back(poly.faces[face_index]);
        }
        bool all = true;
        for (const int& index : box_map) {
            if (index > -1) {
                all = false;
            }
        }
        if (all) {
            for (int i = 0; i < 6; i++) {
                box_map[i] = -2;
            }
        }
        for (int i = 0; i < 6; i++) {
            std::array<std::array<double,3>,4> box;
            switch (i) {
                case 0:
                    box = {{ {{x_min, y_min, z_min}},{{x_min, y_max, z_min}},{{x_min, y_max, z_max}},{{x_min, y_min, z_max}} }};
                    break;
                case 1:
                    box = {{ {{x_max, y_min, z_min}},{{x_max, y_max, z_min}},{{x_max, y_max, z_max}},{{x_max, y_min, z_max}} }};
                    break;
                case 2:
                    box = {{ {{x_min, y_min, z_min}},{{x_max, y_min, z_min}},{{x_max, y_min, z_max}},{{x_min, y_min, z_max}} }};
                    break;
                case 3:
                    box = {{ {{x_min, y_max, z_min}},{{x_max, y_max, z_min}},{{x_max, y_max, z_max}},{{x_min, y_max, z_max}} }};
                    break;
                case 4:
                    box = {{ {{x_min, y_min, z_min}},{{x_max, y_min, z_min}},{{x_max, y_max, z_min}},{{x_min, y_max, z_min}} }};
                    break;
                case 5:
                    box = {{ {{x_min, y_min, z_max}},{{x_max, y_min, z_max}},{{x_max, y_max, z_max}},{{x_min, y_max, z_max}} }};
                    break;
            }
            if (box_map[i] == -1) {
                set<set<std::array<double,3>>> new_face;
                for (int j = 0; j < box.size(); j++) {
                    new_face.insert(Polyhedron::round_edge({box[(j-1+box.size())%box.size()],box[j]}));
                }
                faces.push_back(new_face);
            }
            else if (box_map[i] != -2) {
                bool any1 = false;
                for (const std::array<double,3>& y : box) {
                    bool any2 = false;
                    vector<std::array<double,3>>* exterior_circuit = Polyhedron::find_exterior_circuit(poly.circuits(box_map[i]));
                    for (std::array<std::array<double,3>,3> triangle : Polyhedron::triangulate(*exterior_circuit)) {
                        for (const std::array<double,3>& x : box) {
                            if (Polyhedron::inside_triangle(triangle,x)) {
                                any2 = true;
                                break;
                            }
                        }
                        if (any2) {
                            break;
                        }
                    }
                    delete exterior_circuit;
                    if (!any2) {
                        any1 = true;
                        break; 
                    }
                }
                if (any1) {
                    all = true;
                    for (const std::array<double,3>& y : box) {
                        any1 = false;
                        if (!inside_polyhedron(poly,y)) {
                            all = false;
                            break;
                        }
                        for (int face_index = 0; face_index < poly.faces.size(); face_index++) {
                            for (std::array<std::array<double,3>,3> triangle : Polyhedron::triangulate(*Polyhedron::circuit_cut(poly.circuits(face_index)))) {
                                if (Polyhedron::inside_triangle(triangle,y)) {
                                    any1 = true;
                                    break;
                                }
                            }
                            if (any1) {
                                break;
                            }
                        }
                        if (any1) {
                            all = false;
                            break;
                        }
                    }
                    if (all) {
                        for (int j = 0; j < box.size(); j++) {
                            set<std::array<double,3>> edge1 = {box[(j-1+box.size())%box.size()],box[j]};
                            if (find(old_face_indices.begin(), old_face_indices.end(), box_map[i]) != old_face_indices.end()) {
                                bool no_break = true;
                                for (const set<std::array<double,3>>& edge2 : faces[box_map[i]]) {
                                    if (Polyhedron::point_on_segment(edge2,box[(j-1+box.size())%box.size()]) && Polyhedron::point_on_segment(edge2, box[j])) {
                                        set<std::array<double,3>>::iterator it = edge2.begin();
                                        std::array<double,3> p1 = *(it);
                                        std::array<double,3> p2 = *(++it);
                                        faces[box_map[i]].erase(edge2);
                                        Box::PointDistanceComparator comp = Box::PointDistanceComparator(p1);
                                        std::array<std::array<double,3>,2> edge1_array = {box[(j-1+box.size())%box.size()],box[j]};
                                        sort(edge1_array.begin(),edge1_array.end(), comp);
                                        if (Polyhedron::distance(p1,edge1_array[0]) > 0) {
                                            set<std::array<double,3>> edge3 = {p1, edge1_array[0]};
                                            faces[box_map[i]].insert(Polyhedron::round_edge(edge3));
                                        }
                                        comp = Box::PointDistanceComparator(p2);
                                        sort(edge1_array.begin(),edge1_array.end(), comp);
                                        if (Polyhedron::distance(p2,edge1_array[0]) > 0) {
                                            set<std::array<double,3>> edge3 = {p2, edge1_array[0]};
                                            faces[box_map[i]].insert(Polyhedron::round_edge(edge3));
                                        }
                                        no_break = false;
                                        break;
                                    }
                                }
                                if (no_break) {
                                    faces[box_map[i]].insert(edge1);
                                }
                            }
                        }
                    }
                }
            }
        }
        set<vector<std::array<double,3>>> circuits = Box::del_circuit_helper(path_edges);
        for (const vector<std::array<double,3>>& circuit : circuits) {
            bool no_break = true;
            for (int face_index = 0; face_index < faces.size(); face_index++) {
                if (face_index >= old_face_indices.size()) {
                    break;
                }
                set<std::array<double,3>> point_set;
                for (const std::array<double,3>& point : circuit) {
                    point_set.insert(Polyhedron::round_point(point));
                }
                for (const set<std::array<double,3>>& edge : faces[face_index]) {
                    for (const std::array<double,3>& point : edge) {
                        point_set.insert(Polyhedron::round_point(point));
                    }
                }
                if (Polyhedron::coplanar(point_set)) {
                    vector<std::array<double,3>>* exterior_circuit = Polyhedron::find_exterior_circuit(poly.circuits(old_face_indices[face_index]));
                    for (int i = 0; i < circuit.size(); i++) {
                        if (find(faces[face_index].begin(), faces[face_index].end(), set<std::array<double,3>>{circuit[(i-1+circuit.size())%circuit.size()],circuit[i]}) != faces[face_index].end()) {
                            bool any = false;
                            for (int j = 0; j < exterior_circuit->size(); j++) {
                                if (Polyhedron::point_on_segment({(*exterior_circuit)[(j-1+exterior_circuit->size())%exterior_circuit->size()],(*exterior_circuit)[j]},circuit[(i-1+circuit.size())%circuit.size()]) && Polyhedron::point_on_segment({(*exterior_circuit)[(j-1+exterior_circuit->size())%exterior_circuit->size()]},circuit[i])) {
                                    any = true;
                                    break;
                                }
                            }
                            if (!any) {
                                faces[face_index].erase(Polyhedron::round_edge({circuit[(i-1)%circuit.size()],circuit[i]}));
                            }
                        } else {
                            faces[face_index].insert(Polyhedron::round_edge({circuit[(i-1)%circuit.size()],circuit[i]}));
                        }
                    }
                    delete exterior_circuit;
                    no_break = false;
                    break;
                }
            }
            if (no_break) {
                set<set<std::array<double,3>>> new_face;
                for (int i = 0; i < circuit.size(); i++) {
                    new_face.insert(Polyhedron::round_edge({circuit[(i-1+circuit.size())%circuit.size()],circuit[i]}));
                }
                faces.push_back(new_face);
            }
        }
        while (true) {
            bool triple_break = false;
            for (set<set<std::array<double,3>>>& face : faces) {
                for (const set<std::array<double,3>>& edge1 : face) {
                    for (const set<std::array<double,3>>& edge2 : face) {
                        set<std::array<double,3>> edge_intersection;
                        set_intersection(edge1.begin(), edge1.end(), edge2.begin(), edge2.end(), inserter(edge_intersection, edge_intersection.begin()));
                        set<std::array<double,3>> edge_union;
                        set_union(edge1.begin(), edge1.end(), edge2.begin(), edge2.end(), inserter(edge_union, edge_union.begin()));
                        if (edge_intersection.size() == 1 && Polyhedron::colinear(edge_union)) {
                            map<std::array<double,3>,std::array<double,3>> point_map;
                            for (const std::array<double,3>& point : edge1) {
                                point_map[Polyhedron::round_point(point)] = point;
                            }
                            for (const std::array<double,3>& point : edge2) {
                                point_map[Polyhedron::round_point(point)] = point;
                            }
                            set<std::array<double,3>> edge1_rounded = Polyhedron::round_edge(edge1);
                            set<std::array<double,3>> edge2_rounded = Polyhedron::round_edge(edge2);
                            set<std::array<double,3>> edge_rounded_sym_diff;
                            set_symmetric_difference(edge1_rounded.begin(),edge1_rounded.end(), edge2_rounded.begin(), edge2_rounded.end(), inserter(edge_rounded_sym_diff, edge_rounded_sym_diff.begin()));
                            if (edge_rounded_sym_diff.size() == 2) {
                                set<std::array<double,3>> new_edge;
                                for (const std::array<double,3>& point : edge_rounded_sym_diff) {
                                    new_edge.insert(point_map[point]);
                                }
                                face.insert(new_edge);
                            } else {
                                int max_distance = 0;
                                set<std::array<double,3>> max_edge;
                                for (const std::array<double,3>& x : edge1) {
                                    for (const std::array<double,3>& y : edge2) {
                                        if (Polyhedron::distance(x,y) > max_distance) {
                                            max_distance = Polyhedron::distance(x,y);
                                            max_edge = {x,y};
                                        }
                                    }
                                }
                                face.insert(Polyhedron::round_edge(max_edge));
                            }
                            face.erase(Polyhedron::round_edge(edge1));
                            face.erase(Polyhedron::round_edge(edge2));
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
            if (!triple_break) {
                break;
            }
        }
        int size = faces.size();
        
        for (int face_index = 0; face_index < size; face_index++) {
            circuits = Box::del_circuit_helper(faces[face_index]);
            for (const vector<std::array<double,3>>& circuit : circuits) {
                cout << "del circuit ";
                for (const std::array<double,3>& point : circuit) {
                    cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
                }
                cout << endl;
            }
            vector<std::array<double, 3>>* exterior_circuit = Polyhedron::find_exterior_circuit(Polyhedron::make_clockwise(circuits));
            cout << "exterior circuit is null " << (exterior_circuit==NULL) << endl;
            if (circuits.size() > 1 && exterior_circuit == NULL) {
                vector<vector<std::array<double,3>>> circuits_vector(circuits.begin(),circuits.end());
                for (int i = 1; i < circuits_vector.size(); i++) {
                    set<set<std::array<double,3>>> face;
                    for (int j = 0; j < circuits_vector[i].size(); j++) {
                        face.insert({circuits_vector[i][(j-1+circuits_vector[i].size())%circuits_vector[i].size()], circuits_vector[i][j]});
                        faces[face_index].erase({circuits_vector[i][(j-1+circuits_vector[i].size())%circuits_vector[i].size()], circuits_vector[i][j]});
                    }
                    faces.push_back(face);
                }
            }
        }
        for (int face_index = 0; face_index < faces.size();) {
            if (!faces[face_index].size()) {
                faces.erase(faces.begin() + face_index);
            } else {
                face_index++;
            }
        }
        Polyhedron new_poly;
        set<std::array<double,3>> verts_set;
        set<set<std::array<double,3>>> edges_set;
        for (const set<set<std::array<double,3>>>& face : faces) {
            for (const set<std::array<double,3>>& edge : face) {
                for (const std::array<double,3>& point : edge) {
                    verts_set.insert(point);
                }
                edges_set.insert(edge);
            }
        }
        new_poly.verts.insert(new_poly.verts.end(),verts_set.begin(),verts_set.end());
        vector<set<std::array<double,3>>> edges_vector(edges_set.begin(),edges_set.end());
        for (const set<set<std::array<double,3>>>& face : faces) {
            set<int> face_ungrounded;
            for (const set<std::array<double,3>>& edge : face) {
                vector<set<std::array<double,3>>>::iterator it = find(edges_vector.begin(),edges_vector.end(),edge);
                face_ungrounded.insert(std::distance(edges_vector.begin(), it));
            }
            new_poly.faces.push_back(face_ungrounded);
        }
        for (const set<std::array<double,3>>& edge : edges_vector) {
            set<int> edge_ungrounded;
            for (const std::array<double,3>& point : edge) {
                vector<std::array<double,3>>::iterator it = find(new_poly.verts.begin(),new_poly.verts.end(),point);
                edge_ungrounded.insert(std::distance(new_poly.verts.begin(),it));
            }
            new_poly.edges.push_back(edge_ungrounded);
        }
        return new_poly;
    }
    static set<vector<std::array<double,3>>> add_circuit_helper(set<set<std::array<double,3>>>face, set<std::array<double,3>> start, set<std::array<double,3>> previous, set<std::array<double,3>> current, vector<std::array<double,3>> path) {
        map<std::array<double,3>,set<set<std::array<double,3>>>> edge_lookup;
        for (const set<std::array<double,3>>& edge : face) {
            for (const std::array<double,3>& point : edge) {
                edge_lookup[Polyhedron::round_point(point)].insert(edge);
            }
        }
        set<vector<std::array<double,3>>> output;
        if (current == start) {
            output.insert(path);
            return output;
        }
        set<std::array<double,3>> difference;
        set_difference(current.begin(),current.end(),previous.begin(),previous.end(),inserter(difference,difference.begin()));
        std::array<double,3> point = *(difference.begin());
        if (find(path.begin(), path.end(), point) != path.end()) {
            return output;
        }
        path.push_back(point);
        previous = current;
        set<set<std::array<double,3>>> temp(edge_lookup[Polyhedron::round_point(point)].begin(),edge_lookup[Polyhedron::round_point(point)].end());
        temp.erase(previous);
        for (const set<std::array<double,3>>& x : temp) {
            set<std::array<double,3>> point_union;
            for (const std::array<double,3>& y : path) {
                point_union.insert(Polyhedron::round_point(y));
            }
            set<std::array<double,3>> edge = Polyhedron::round_edge(x);
            point_union.insert(edge.begin(),edge.end());
            if (Polyhedron::coplanar(point_union)) {
                current = x;
                set<vector<std::array<double,3>>> intermediate = Box::add_circuit_helper(face, start, previous, current, path);
                output.insert(intermediate.begin(), intermediate.end());
            }
        }
        vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
        for (int i = 0; i < output_list.size(); i++) {
            for (int j = 0; j < output_list.size(); j++) {
                set<std::array<double,3>> difference2;
                set_difference(output_list[i].begin(),output_list[i].end(),output_list[j].begin(),output_list[j].end(),inserter(difference2,difference2.begin()));
                if (output_list[j].size() > output_list[i].size() && !difference2.size() && output.find(output_list[j]) != output.end()) {
                    output.erase(output_list[j]);
                }
                if (output_list[j].size() == output_list[i].size() && !difference2.size() && j > i && output.find(output_list[j]) != output.end()) {
                    output.erase(output_list[j]);
                }
                if (j > i && !(difference2.size()) && output.find(output_list[j]) != output.end()) {
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
    static set<vector<std::array<double,3>>> add_circuit_helper(set<set<std::array<double,3>>> face) {
        map<std::array<double,3>,set<set<std::array<double,3>>>> edge_lookup;
        for (const set<std::array<double,3>>& edge : face) {
            for (const std::array<double,3>& point : edge) {
                edge_lookup[Polyhedron::round_point(point)].insert(edge);
            }
        }
        set<vector<std::array<double,3>>> output;
        for (const set<std::array<double,3>>& edge : face) {
            vector<std::array<double,3>> path;
            set<std::array<double,3>> start = edge;
            std::array<double,3> point = *(start.begin());
            path.push_back(point);
            set<set<std::array<double,3>>> temp(edge_lookup[Polyhedron::round_point(point)].begin(),edge_lookup[Polyhedron::round_point(point)].end());
            temp.erase(start);
            for (const set<std::array<double,3>>& x : temp) {
                set<std::array<double,3>> current = x;
                set<vector<std::array<double,3>>> intermediate = Box::add_circuit_helper(face, start, start, current, path);
                output.insert(intermediate.begin(),intermediate.end());
            }
        }
        vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
        for (int i = 0; i < output_list.size(); i++) {
            for (int j = 0; j < output_list.size(); j++) {
                set<std::array<double,3>> difference2;
                set_difference(output_list[i].begin(),output_list[i].end(),output_list[j].begin(),output_list[j].end(),inserter(difference2,difference2.begin()));
                if (output_list[j].size() > output_list[i].size() && !difference2.size() && output.find(output_list[j]) != output.end()) {
                    output.erase(output_list[j]);
                }
                if (output_list[j].size() == output_list[i].size() && !difference2.size() && j > i && output.find(output_list[j]) != output.end()) {
                    output.erase(output_list[j]);
                }
                if (j > i && !(difference2.size()) && output.find(output_list[j]) != output.end()) {
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
    Polyhedron add(Polyhedron poly) const {
        vector<set<std::array<double,3>>> edges;
        for (const set<int>& edge : poly.edges) {
            set<std::array<double,3>> edge_grounded;
            for (const int& index : edge) {
                edge_grounded.insert(poly.verts[index]);
            }
            edges.push_back(edge_grounded);
        }
        vector<set<set<std::array<double,3>>>> faces;
        for (const set<int>& face : poly.faces) {
            set<set<std::array<double,3>>> face_grounded;
            for (const int& index : face) {
                face_grounded.insert(edges[index]);
            }
            faces.push_back(face_grounded);
        }
        vector<set<std::array<double,3>>> new_edges;
        set<std::array<double,3>> new_edge;
        std::array<set<set<std::array<double,3>>>,6> new_faces;
        new_edge = {{x_min,y_min,z_min},{x_max,y_min,z_min}};
        new_edges.push_back(new_edge);
        new_faces[2].insert(new_edge);
        new_faces[4].insert(new_edge);
        new_edge = {{x_min,y_max,z_min},{x_max,y_max,z_min}};
        new_edges.push_back(new_edge);
        new_faces[3].insert(new_edge);
        new_faces[4].insert(new_edge);
        new_edge = {{x_min,y_min,z_max},{x_max,y_min,z_max}};
        new_edges.push_back(new_edge);
        new_faces[2].insert(new_edge);
        new_faces[5].insert(new_edge);
        new_edge = {{x_min,y_max,z_max},{x_max,y_max,z_max}};
        new_edges.push_back(new_edge);
        new_faces[3].insert(new_edge);
        new_faces[5].insert(new_edge);
        new_edge = {{x_min,y_min,z_min},{x_min,y_max,z_min}};
        new_edges.push_back(new_edge);
        new_faces[0].insert(new_edge);
        new_faces[4].insert(new_edge);
        new_edge = {{x_max,y_min,z_min},{x_max,y_max,z_min}};
        new_edges.push_back(new_edge);
        new_faces[1].insert(new_edge);
        new_faces[4].insert(new_edge);
        new_edge = {{x_min,y_min,z_max},{x_min,y_max,z_max}};
        new_edges.push_back(new_edge);
        new_faces[0].insert(new_edge);
        new_faces[5].insert(new_edge);
        new_edge = {{x_max,y_min,z_max},{x_max,y_max,z_max}};
        new_edges.push_back(new_edge);
        new_faces[1].insert(new_edge);
        new_faces[5].insert(new_edge);
        new_edge = {{x_min,y_min,z_min},{x_min,y_min,z_max}};
        new_edges.push_back(new_edge);
        new_faces[0].insert(new_edge);
        new_faces[2].insert(new_edge);
        new_edge = {{x_max,y_min,z_min},{x_max,y_min,z_max}};
        new_edges.push_back(new_edge);
        new_faces[1].insert(new_edge);
        new_faces[2].insert(new_edge);
        new_edge = {{x_min,y_max,z_min},{x_min,y_max,z_max}};
        new_edges.push_back(new_edge);
        new_faces[0].insert(new_edge);
        new_faces[3].insert(new_edge);
        new_edge = {{x_max,y_max,z_min},{x_max,y_max,z_max}};
        new_edges.push_back(new_edge);
        new_faces[1].insert(new_edge);
        new_faces[3].insert(new_edge);
        map<int,set<set<set<std::array<double,3>>>>> mapping;
        for (int face_index1 = 0; face_index1 < faces.size(); face_index1++) {
            vector<std::array<double,3>>* circuit = Polyhedron::circuit_cut(Polyhedron::make_clockwise(poly.circuits(face_index1)));
            bool any = false;
            for (const set<std::array<double,3>>& edge : faces[face_index1]) {
                if (this->intersect(edge)) {
                    any = true;
                    break;
                }
            }
            if (any) {
                for (int face_index2 = 0; face_index2 < new_faces.size(); face_index2++) {
                    any = false;
                    for (int i = 0; i < circuit->size(); i++) {
                        for (const set<std::array<double,3>>& y : new_faces[face_index2]) {
                            if (find(poly.verts.begin(),poly.verts.end(), (*circuit)[(i-1+circuit->size())%circuit->size()]) != poly.verts.end() && find(poly.verts.begin(),poly.verts.end(), (*circuit)[i]) != poly.verts.end()) {
                                bool any2 = false;
                                for (const set<std::array<double,3>>& z : edges) {
                                    if (z.find((*circuit)[(i-1+circuit->size())%circuit->size()]) != z.end() && z.find((*circuit)[i]) != z.end()) {
                                        any2 = true;
                                        break;
                                    }
                                }
                                if (Box::intersect_segments({(*circuit)[(i-1+circuit->size())%circuit->size()],(*circuit)[i]},y) != NULL && any2) {
                                    any = true;
                                    break;
                                }
                            }
                        }
                        if (any) {
                            break; 
                        }
                    }
                    set<std::array<double,3>> point_union;
                    for (const set<std::array<double,3>>& edge : faces[face_index1]) {
                        for (const std::array<double,3>& point : edge) {
                            point_union.insert(Polyhedron::round_point(point));
                        }
                    }
                    for (const set<std::array<double,3>>& edge : new_faces[face_index2]) {
                        for (const std::array<double,3>& point : edge) {
                            point_union.insert(Polyhedron::round_point(point));
                        }
                    }
                    if (Polyhedron::coplanar(point_union) && any) {
                        set<set<std::array<double,3>>> face_union(new_faces[face_index2/2].begin(),(new_faces[face_index2/2].end()));
                        face_union.insert(new_faces[face_index2/2+1].begin(),(new_faces[face_index2/2+1].end()));
                        set<set<std::array<double,3>>> difference;
                        set_difference(new_edges.begin(),new_edges.end(),face_union.begin(),face_union.end(), inserter(difference, difference.begin()));
                        bool double_break = false;
                        bool triple_break = false;
                        for (const std::array<double,3>& p1 : (*circuit)) {
                            if (find(poly.verts.begin(), poly.verts.end(), p1) != poly.verts.end()) {
                                set<std::array<double,3>> point_set;
                                for (const set<std::array<double,3>>& edge : faces[face_index1]) {
                                    for (const std::array<double,3>& point : edge) {
                                        point_set.insert(point);
                                    }
                                }
                                for (const set<std::array<double,3>>& edge1 : edges) {
                                    if (edge1.find(p1) != edge1.end()) {
                                        std::array<double,3> p2;
                                        for (const std::array<double,3>& p : edge1) {
                                            if (p != p1) {
                                                p2 = p;
                                            }
                                        }
                                        if (point_set.find(p2) == point_set.end()) {
                                            for (const set<std::array<double,3>>& edge2 : difference) {
                                                set<std::array<double,3>>::iterator it = edge2.begin();
                                                std::array<double,3> p3 = *(it);
                                                it++;
                                                std::array<double,3> p4 = *(it);
                                                if (Polyhedron::point_on_segment(edge2, p1) && Polyhedron::point_on_segment(edge2, p2) || (Polyhedron::point_on_segment(edge1, p3) && Polyhedron::point_on_segment(edge1, p4))) {
                                                    triple_break = true;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if (triple_break) {
                                        break;
                                    }
                                }
                                for (const set<std::array<double,3>>& edge1 : edges) {
                                    std::array<double,3> p2;
                                    for (const std::array<double,3>& p : edge1) {
                                        if (p != p1) {
                                            p2 = p;
                                        }
                                    }
                                    if (point_set.find(p2) == point_set.end() && this->intersect(edge1)) {
                                        double_break = true;
                                        break;
                                    }
                                }
                            }
                            if (double_break || triple_break) {
                                break;
                            }
                        }
                        if (double_break || triple_break) {
                            mapping[face_index2].insert(faces[face_index1]);
                        }
                    }
                }
            }
            for (int face_index2 = 0; face_index2 < new_faces.size(); face_index2++) {
                set<std::array<double,3>> point_union;
                for (const set<std::array<double,3>>& edge : faces[face_index1]) {
                    for (const std::array<double,3>& point : edge) {
                        point_union.insert(Polyhedron::round_point(point));
                    }
                }
                for (const set<std::array<double,3>>& edge : new_faces[face_index2]) {
                    for (const std::array<double,3>& point : edge) {
                        point_union.insert(Polyhedron::round_point(point));
                    }
                }
                any = false;
                for (int i = 0; i < circuit->size(); i++) {
                    for (const set<std::array<double,3>>& y : new_faces[face_index2]) {
                        if (find(poly.verts.begin(),poly.verts.end(), (*circuit)[(i-1+circuit->size())%circuit->size()]) != poly.verts.end() && find(poly.verts.begin(),poly.verts.end(), (*circuit)[i]) != poly.verts.end()) {
                            bool any2 = false;
                            for (const set<std::array<double,3>>& z : edges) {
                                if (z.find((*circuit)[(i-1+circuit->size())%circuit->size()]) != z.end() && z.find((*circuit)[i]) != z.end()) {
                                    any2 = true;
                                    break;
                                }
                            }
                            if (Box::intersect_segments({(*circuit)[(i-1+circuit->size())%circuit->size()],(*circuit)[i]},y) != NULL && any2) {
                                any = true;
                                break;
                            }
                        }
                    }
                }
                set<std::array<double,3>> point_set;
                for (const set<std::array<double,3>>& edge : new_faces[face_index2]) {
                    for (const std::array<double,3>& point : edge) {
                        point_set.insert(point);
                    }
                }
                bool all1 = true;
                for (const std::array<double,3>& x : point_set) {
                    bool any2 = false;
                    for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*circuit)) {
                        if (Polyhedron::inside_triangle(triangle, x)) {
                            any2 = true;
                            break;
                        }
                    }
                    if (!any2) {
                        all1 = false;
                        break;
                    }
                }
                bool all2 = true; 
                for (const std::array<double,3>& x : *circuit) {
                    if (x[0] < x_min || x[0] > x_max || x[1] < y_min || x[1] > y_max || x[2] < z_min || x[2] > z_max) {
                        all2 = false;
                        break;
                    }
                }
                if (Polyhedron::coplanar(point_union) && (any || all1 || all2)) {
                    mapping[face_index2].insert(faces[face_index1]);
                }
            }
        }
        cout << "mapping" << endl;
        for (const pair<int,set<set<set<std::array<double,3>>>>>& p : mapping) {
            cout << p.first << endl;
            for (const set<set<std::array<double,3>>>& face : p.second) {
                cout << "face" << endl;
                for (const set<std::array<double,3>>& edge : face) {
                    for (const std::array<double,3>& point : edge) {
                        cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
                    }
                    cout << endl;
                }
            }
        }
        set<set<std::array<double,3>>> erased;
        for (const pair<int,set<set<set<std::array<double,3>>>>>& p : mapping) {
            for (const set<set<std::array<double,3>>>& face2 : p.second) {
                for (const set<std::array<double,3>>& edge1 : face2) {
                    bool no_break1 = true;
                    for (const set<std::array<double,3>>& edge2 : new_faces[p.first]) {
                        bool no_break2 = true;
                        for (const std::array<double,3>& point : edge2) {
                            if (!Polyhedron::point_on_segment(Polyhedron::round_edge(edge1),Polyhedron::round_point(point))) {
                                no_break2 = false;
                                break;
                            }
                        }
                        if (no_break2) {
                            set<std::array<double,3>>::iterator it = edge1.begin();
                            std::array<double,3> p1 = *(it);
                            std::array<double,3> p2 = *(++it);
                            it = edge2.begin();
                            std::array<std::array<double,3>,2> edge2_array;
                            edge2_array[0] = *(it);
                            edge2_array[1] = *(++it);
                            cout << "ERASE [" << edge2_array[0][0] << "," << edge2_array[0][1] << "," << edge2_array[0][2] << "] [" << edge2_array[1][0] << "," << edge2_array[1][1] << "," << edge2_array[1][2] << "]" << endl;
                            erased.insert(edge2);
                            new_faces[p.first].erase(edge2);
                            PointDistanceComparator comp = PointDistanceComparator(p1);
                            sort(edge2_array.begin(),edge2_array.end(),comp);
                            if (Polyhedron::distance(p1,edge2_array[0])) {
                                new_faces[p.first].insert({p1,edge2_array[0]});
                            }
                            comp = PointDistanceComparator(p2);
                            sort(edge2_array.begin(),edge2_array.end(),comp);
                            if (Polyhedron::distance(p2,edge2_array[0])) {
                                new_faces[p.first].insert({p2,edge2_array[0]});
                            }
                            no_break1 = false;
                            break;
                        }
                        no_break2 = true;
                        for (const std::array<double,3>& point : edge1) {
                            if (!Polyhedron::point_on_segment(Polyhedron::round_edge(edge2),Polyhedron::round_point(point))) {
                                no_break2 = false;
                                break;
                            }
                        }
                        if (no_break2) {
                            set<std::array<double,3>>::iterator it = edge2.begin();
                            std::array<double,3> p1 = *(it);
                            it++;
                            std::array<double,3> p2 = *(it);
                            new_faces[p.first].erase(edge2);
                            it = edge1.begin();
                            std::array<std::array<double,3>,2> edge1_array;
                            edge1_array[0] = *(it);
                            edge1_array[1] = *(++it);
                            PointDistanceComparator comp = PointDistanceComparator(p1);
                            sort(edge1_array.begin(),edge1_array.end(),comp);
                            if (Polyhedron::distance(p1,edge1_array[0])) {
                                new_faces[p.first].insert({p1,edge1_array[0]});
                            }
                            comp = PointDistanceComparator(p2);
                            sort(edge1_array.begin(),edge1_array.end(),comp);
                            if (Polyhedron::distance(p2,edge1_array[0])) {
                                new_faces[p.first].insert({p2,edge1_array[0]});
                            }
                            no_break1 = false;
                            break;
                        }
                        set<std::array<double,3>> edge_union;
                        set_union(edge1.begin(),edge1.end(),edge2.begin(),edge2.end(),inserter(edge_union, edge_union.begin()));
                        if (Polyhedron::colinear(edge_union)) {
                            set<std::array<double,3>>::iterator it = edge1.begin();
                            std::array<double,3> p1 = *(it);
                            it++;
                            std::array<double,3> p2 = *(it);
                            it = edge2.begin();
                            std::array<double,3> p3 = *(it);
                            it++;
                            std::array<double,3> p4 = *(it);
                            std::array<std::array<double,3>,2> edge1_array;
                            edge1_array[0] = p1;
                            edge1_array[1] = p2;
                            std::array<std::array<double,3>,2> edge2_array;
                            edge2_array[0] = p3;
                            edge2_array[1] = p4;
                            set<std::array<double,3>> edge_intersection;
                            set_intersection(edge1.begin(),edge1.end(),edge2.begin(),edge2.end(),inserter(edge_intersection, edge_intersection.begin()));
                            cout << endl;
                            cout << "[" << p1[0] << "," << p1[1] << "," << p1[2] << "] ";
                            cout << "[" << p2[0] << "," << p2[1] << "," << p2[2] << "] " << endl;
                            cout << "[" << p3[0] << "," << p3[1] << "," << p3[2] << "] ";
                            cout << "[" << p4[0] << "," << p4[1] << "," << p4[2] << "] " << endl;
                            if (Polyhedron::point_on_segment(edge2, p1) and Polyhedron::point_on_segment(edge1, p3)) {
                                new_faces[p.first].erase(edge2);
                                PointDistanceComparator comp = PointDistanceComparator(p2);
                                sort(edge2_array.begin(),edge2_array.end(),comp);
                                std::array<double,3> p5 = edge2_array[0];
                                comp = PointDistanceComparator(p4);
                                sort(edge1_array.begin(),edge1_array.end(),comp);
                                std::array<double,3> p6 = edge1_array[0];
                                if (erased.find({p2,p5}) == erased.end()) {
                                    new_faces[p.first].insert({p2,p5});
                                }
                                if (erased.find({p4,p6}) == erased.end()) {
                                    new_faces[p.first].insert({p4,p6});
                                }
                                cout << "INSERT ";
                                cout << "[" << p2[0] << "," << p2[1] << "," << p2[2] << "] ";
                                cout << "[" << p5[0] << "," << p5[1] << "," << p5[2] << "] " << endl;
                                cout << "[" << p4[0] << "," << p4[1] << "," << p4[2] << "] ";
                                cout << "[" << p6[0] << "," << p6[1] << "," << p6[2] << "] " << endl;
                                no_break1 = false;
                                break;
                            }
                            if (Polyhedron::point_on_segment(edge2, p2) and Polyhedron::point_on_segment(edge1, p3)) {
                                new_faces[p.first].erase(edge2);
                                PointDistanceComparator comp = PointDistanceComparator(p1);
                                sort(edge2_array.begin(),edge2_array.end(),comp);
                                std::array<double,3> p5 = edge2_array[0];
                                comp = PointDistanceComparator(p4);
                                sort(edge1_array.begin(),edge1_array.end(),comp);
                                std::array<double,3> p6 = edge1_array[0];
                                if (erased.find({p1,p5}) == erased.end()) {
                                    new_faces[p.first].insert({p1,p5});
                                }
                                if (erased.find({p4,p6}) == erased.end()) {
                                    new_faces[p.first].insert({p4,p6});
                                }
                                cout << "INSERT ";
                                cout << "[" << p1[0] << "," << p1[1] << "," << p1[2] << "] ";
                                cout << "[" << p5[0] << "," << p5[1] << "," << p5[2] << "] " << endl;
                                cout << "[" << p4[0] << "," << p4[1] << "," << p4[2] << "] ";
                                cout << "[" << p6[0] << "," << p6[1] << "," << p6[2] << "] " << endl;
                                no_break1 = false;
                                break;
                            }
                            if (Polyhedron::point_on_segment(edge2, p1) and Polyhedron::point_on_segment(edge1, p4)) {
                                new_faces[p.first].erase(edge2);
                                PointDistanceComparator comp = PointDistanceComparator(p2);
                                sort(edge2_array.begin(),edge2_array.end(),comp);
                                std::array<double,3> p5 = edge2_array[0];
                                comp = PointDistanceComparator(p3);
                                sort(edge1_array.begin(),edge1_array.end(),comp);
                                std::array<double,3> p6 = edge1_array[0];
                                if (erased.find({p2,p5}) == erased.end()) {
                                    new_faces[p.first].insert({p2,p5});
                                }
                                if (erased.find({p3,p6}) == erased.end()) {
                                    new_faces[p.first].insert({p3,p6});
                                }
                                cout << "INSERT ";
                                cout << "[" << p2[0] << "," << p2[1] << "," << p2[2] << "] ";
                                cout << "[" << p5[0] << "," << p5[1] << "," << p5[2] << "] " << endl;
                                cout << "[" << p3[0] << "," << p3[1] << "," << p3[2] << "] ";
                                cout << "[" << p6[0] << "," << p6[1] << "," << p6[2] << "] " << endl;
                                no_break1 = false;
                                break;
                            }
                            if (Polyhedron::point_on_segment(edge2, p2) and Polyhedron::point_on_segment(edge1, p4)) {
                                new_faces[p.first].erase(edge2);
                                PointDistanceComparator comp = PointDistanceComparator(p1);
                                sort(edge2_array.begin(),edge2_array.end(),comp);
                                std::array<double,3> p5 = edge2_array[0];
                                comp = PointDistanceComparator(p3);
                                sort(edge1_array.begin(),edge1_array.end(),comp);
                                std::array<double,3> p6 = edge1_array[0];
                                if (erased.find({p1,p5}) == erased.end()) {
                                    new_faces[p.first].insert({p1,p5});
                                }
                                if (erased.find({p3,p6}) == erased.end()) {
                                    new_faces[p.first].insert({p3,p6});
                                }
                                cout << "INSERT ";
                                cout << "[" << p1[0] << "," << p1[1] << "," << p1[2] << "] ";
                                cout << "[" << p5[0] << "," << p5[1] << "," << p5[2] << "] " << endl;
                                cout << "[" << p3[0] << "," << p3[1] << "," << p3[2] << "] ";
                                cout << "[" << p6[0] << "," << p6[1] << "," << p6[2] << "] " << endl;
                                no_break1 = false;
                                break;
                            }
                        }
                    }
                    if (no_break1) {
                        cout << "INSERT [" << (*edge1.begin())[0] << "," << (*edge1.begin())[1] << "," << (*edge1.begin())[2] << "] [" << (*++edge1.begin())[0] << "," << (*++edge1.begin())[1] << "," << (*++edge1.begin())[2] << "]" << endl;
                        if (erased.find(edge1) == erased.end()) {
                            new_faces[p.first].insert(edge1);
                        }
                    }
                }
            }
            for (int face_index = 0; face_index < faces.size();) {
                if (p.second.find(faces[face_index]) != p.second.end()) {
                    faces.erase(faces.begin()+face_index);
                } else {
                    face_index++;
                }
            }
        }
        vector<set<set<std::array<double,3>>>> new_faces_vector;
        for (const set<set<std::array<double,3>>>& face : new_faces) {
            if (face.size()) {
                new_faces_vector.push_back(face);
            }
        }
        faces.insert(faces.end(),new_faces_vector.begin(),new_faces_vector.end());
        int size = faces.size();
        for (int face_index = 0; face_index < size; face_index++) {
            set<vector<std::array<double,3>>> circuits = Box::add_circuit_helper(faces[face_index]);
            for (const vector<std::array<double,3>>& circuit : circuits) {
                cout << "circuit " << endl;
                for (const std::array<double,3>& point : circuit) {
                    cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
                }
                cout << endl;
            }
            vector<std::array<double,3>>* exterior_circuit = Polyhedron::find_exterior_circuit(circuits);
	    cout << "null exterior " << (exterior_circuit == NULL) << endl;
            if (circuits.size() > 1 && exterior_circuit == NULL) {
                vector<vector<std::array<double,3>>> circuits_vector(circuits.begin(),circuits.end());
                for (int circuit_index = 1; circuit_index < circuits_vector.size(); circuit_index++) {
                    set<set<std::array<double,3>>> new_face;
                    for (int i = 0; i < circuits_vector[circuit_index].size(); i++) {
                        set<std::array<double,3>> edge = {circuits_vector[circuit_index][(i-1+circuits_vector[circuit_index].size())%circuits_vector[circuit_index].size()],circuits_vector[circuit_index][i]};
                        new_face.insert(edge);
                        faces[face_index].erase(edge);
                    }
                    faces.push_back(new_face);
                }
            }
            delete exterior_circuit;
        }
        while (true) {
            bool triple_break = false;
            for (set<set<std::array<double,3>>>& face : faces) {
                for (const set<std::array<double,3>>& edge1: face) {
                    for (const set<std::array<double,3>>& edge2: face) {
                        set<std::array<double,3>> edge_intersection;
                        set_intersection(edge1.begin(),edge1.end(),edge2.begin(),edge2.end(),inserter(edge_intersection, edge_intersection.begin()));
                        set<std::array<double,3>> edge_union;
                        set_union(edge1.begin(),edge1.end(),edge2.begin(),edge2.end(),inserter(edge_union, edge_union.begin()));
                        if (edge_intersection.size()==1 && Polyhedron::colinear(edge_union)) {
                            bool to_continue = false;
                            for (const set<std::array<double,3>>& edge3: face) {
                                if (edge3 != edge1 && edge3 != edge2 && edge3.find(*(edge_intersection.begin())) != edge3.end()) {
                                    to_continue = true;
                                    break;
                                }
                            }
                            if (to_continue) {
                                continue;
                            }
                            map<std::array<double,3>,std::array<double,3>> point_map;
                            for (const std::array<double,3>& point : edge1) {
                                point_map[Polyhedron::round_point(point)] = point;
                            }
                            for (const std::array<double,3>& point : edge2) {
                                point_map[Polyhedron::round_point(point)] = point;
                            }
                            set<std::array<double,3>> edge1_rounded = Polyhedron::round_edge(edge1);
                            set<std::array<double,3>> edge2_rounded = Polyhedron::round_edge(edge2);
                            set<std::array<double,3>> edge_rounded_sym_diff;
                            set_symmetric_difference(edge1_rounded.begin(),edge1_rounded.end(), edge2_rounded.begin(), edge2_rounded.end(), inserter(edge_rounded_sym_diff, edge_rounded_sym_diff.begin()));
                            if (edge_rounded_sym_diff.size() == 2) {
                                set<std::array<double,3>> new_edge;
                                for (const std::array<double,3>& point : edge_rounded_sym_diff) {
                                    new_edge.insert(point_map[point]);
                                }
                                face.insert(new_edge);
                            } else {
                                int max_distance = 0;
                                set<std::array<double,3>> max_edge;
                                for (const std::array<double,3>& x : edge1) {
                                    for (const std::array<double,3>& y : edge2) {
                                        if (Polyhedron::distance(x,y) > max_distance) {
                                            max_distance = Polyhedron::distance(x,y);
                                            max_edge = {x,y};
                                        }
                                    }
                                }
                                face.insert(max_edge);
                            }
                            face.erase(edge1);
                            face.erase(edge2);
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
            if (!triple_break) {
                break;
            }
        }
        for (const set<set<std::array<double,3>>>& face : faces) {
            cout << face.size() << " ";
        }
        cout << endl;
        for (int face_index = 0; face_index < faces.size();) {
            if (!faces[face_index].size()) {
                faces.erase(faces.begin()+face_index);
            } else {
                face_index++;
            }
        }
        Polyhedron new_poly;
        set<std::array<double,3>> verts_set;
        set<set<std::array<double,3>>> edges_set;
        for (const set<set<std::array<double,3>>>& face : faces) {
            for (const set<std::array<double,3>>& edge : face) {
                for (const std::array<double,3>& point : edge) {
                    verts_set.insert(point);
                }
                edges_set.insert(edge);
            }
        }
        new_poly.verts.insert(new_poly.verts.end(),verts_set.begin(),verts_set.end());
        vector<set<std::array<double,3>>> edges_vector(edges_set.begin(),edges_set.end());
        for (const set<set<std::array<double,3>>>& face : faces) {
            set<int> face_ungrounded;
            for (const set<std::array<double,3>>& edge : face) {
                vector<set<std::array<double,3>>>::iterator it = find(edges_vector.begin(),edges_vector.end(),edge);
                face_ungrounded.insert(std::distance(edges_vector.begin(), it));
            }
            new_poly.faces.push_back(face_ungrounded);
        }
        for (const set<std::array<double,3>>& edge : edges_vector) {
            set<int> edge_ungrounded;
            for (const std::array<double,3>& point : edge) {
                vector<std::array<double,3>>::iterator it = find(new_poly.verts.begin(),new_poly.verts.end(),point);
                edge_ungrounded.insert(std::distance(new_poly.verts.begin(),it));
            }
            new_poly.edges.push_back(edge_ungrounded);
        }
        return new_poly;
    }
};
Polyhedron get_cube(std::array<double,3>displacement, std::array<double,3>factors, std::array<double,3> angles) {
        vector<std::array<double,3>> points = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};
        std::array<double,3> center = {0.5,0.5,0.5};
        vector<std::array<double,3>> temp1 = Polyhedron::shift(points, {-center[0], -center[1], -center[2]});
        vector<std::array<double,3>> temp2 = Polyhedron::scale(temp1, factors);
        vector<std::array<double,3>> temp3 = Polyhedron::rotate(temp2, angles);
        vector<std::array<double,3>> temp4 = Polyhedron::shift(temp3, displacement);
        vector<std::set<int>> edges;
        for (int i = 0; i < 4; i++) {
            edges.push_back({i, (i+1)%4});
        }
        for (int i = 0; i < 4; i++) {
            edges.push_back({4+i, 4+(i+1)%4});
        }
        for (int i = 0; i < 4; i++) {
            edges.push_back({i, i+4});
        }
        Polyhedron polyhedron;
        polyhedron.verts = temp4;
        polyhedron.edges = edges;
        polyhedron.faces = {{0,1,2,3},{4,5,6,7},{8,4,9,0},{10,6,11,2},{1,10,5,9},{3,11,7,8}};
        return polyhedron;
}
Polyhedron get_cube(std::array<double,3>displacement, std::array<double,3>factors) {
    return get_cube(displacement, factors, {0,0,0});
}
// Screen dimension constants
const int SCREEN_WIDTH = 1280;//640;
const int SCREEN_HEIGHT = 960;//480;
Camera camera;

class Crosshair {
  public:
    std::array<double,2> pos = {0,0};
    void draw(SDL_Renderer* gRenderer) {
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
        SDL_RenderDrawLine(gRenderer, pos[0]+SCREEN_WIDTH/2-10, -pos[1]+SCREEN_HEIGHT/2, pos[0]+SCREEN_WIDTH/2+10,-pos[1]+SCREEN_HEIGHT/2);
        SDL_RenderDrawLine(gRenderer, pos[0]+SCREEN_WIDTH/2, -pos[1]+SCREEN_HEIGHT/2-10, pos[0]+SCREEN_WIDTH/2,-pos[1]+SCREEN_HEIGHT/2+10);
        }
};

int fill_polygon(SDL_Renderer* gRenderer, vector<std::array<double,2>> points) {
    if (points.size() == 0) {
        return 0;
    }
    if (points.size() == 1) {
        return SDL_RenderDrawPoint(gRenderer, round(points.front()[0]), round(points.front()[1]));
    }
    double y0 = points[0][1];
    double y1 = y0;
    for (const std::array<double,2>& point : points) {
        if (point[1] < y0) {
            y0 = point[1];
        }
        if (point[1] > y1) {
            y1 = point[1];
        }
    }
    for (int y = round(y0); y <= round(y1); y++) {
        int j = points.size()-1;
        vector<int> xs;
        for (int i= 0;  i < points.size();  j = i++) {
            if ((points[i][1] < y && y <= points[j][1]) || (points[j][1] < y && y <= points[i][1])) {
                // the x value of the vertex + the height delta between the vertex and the scanline *the reciprocal of the slope =
                // the x value of the intersection of the preceding edge and the scanline 
                xs.push_back(round(points[i][0] + (y-points[i][1]) * (points[i][0]-points[j][0]) / (points[i][1]-points[j][1]) ));
            }
            // sort the intersections, the new item bubbles up to its place
            for (int k = xs.size()-1;  k > 0 && xs[k-1] > xs[k];  --k) {
	            swap(xs[k-1], xs[k]);
	        }
        }
        
        for (int i = 0; i < xs.size(); i += 2) {
            SDL_RenderDrawLine(gRenderer, xs[i], y, xs[i+1], y);
        }
    }
    return 0;
}

namespace voxel_editor {
class Block {
  public:
    std::array<double,3> origin;
    std::array<double,3> size;
    Polyhedron block;
    Polyhedron poly;
    std::array<double,3> select = {0,0,0};
    int select_dimension = 0;
    double unit;
    std::array<double,3> select_size;
    set<double> meters = {1};
    set<vector<std::array<double,3>>> polygons;
    bool space_down = false;
    bool cylinder_tool = false;
    int cylinder_dimension = 0;
    set<string> cylinder_options;
    Polyhedron cylinder_poly;
    vector<Box> cylinder_boxes;
    Block (std::array<double,2> x, std::array<double,2> y, std::array<double,2> z, double unit) {
        origin = {x[0],y[0],z[0]};
        size = {x[1]-x[0], y[1]-y[0], z[1]-z[0]};
        block = get_cube({size[0]/2,origin[1]+size[1]/2,size[2]/2},{size[0],size[1],size[2]});
        this->unit = unit;
        select_size = {unit,unit,unit};
    }
    Block (double width, double height, double depth, double unit) : Block({0,width},{0,height},{0,depth}, unit) {
    }
    void flip() {
        cout << "flip" << endl;
        polygons = {};
        for (int face_index = 0; face_index < poly.faces.size(); face_index++) {
            cout << poly.faces[face_index].size() << endl;
            for (const int& edge_index : poly.faces[face_index]) {
                for(const int& index : poly.edges[edge_index]) {
                    cout << index << " [" << poly.verts[index][0] << "," << poly.verts[index][1] << "," << poly.verts[index][2] << "] ";
                }
                cout << endl;
            }
            set<vector<std::array<double,3>>> circuits = poly.circuits(face_index);
            for (const vector<std::array<double,3>>& circuit : circuits) {
                cout << "circuit " << endl;
                for (const std::array<double,3>& point : circuit) {
                    cout << "[" << point[0] << "," << point[1] << "," << point[2] << "] ";
                }
                cout << endl;
            }
            vector<std::array<double,3>>* points = Polyhedron::circuit_cut(Polyhedron::make_clockwise(poly.circuits(face_index)));
            for (const std::array<double,3> point : *points) {
                cout << " [" << point[0] << "," << point[1] << "," << point[2] << "] ";
            }
            cout << endl;
            polygons.insert(*points);
            delete points;
        }
    }
    void draw(SDL_Renderer* gRenderer) {
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
        for (const set<int>& edge : block.edges) {
            std::array<double,2> p1 = camera.project(block.verts[*(edge.begin())]);
            std::array<double,2> p2 = camera.project(block.verts[*(++edge.begin())]);
            p1[0] += SCREEN_WIDTH/2;
            p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
            p2[0] += SCREEN_WIDTH/2;
            p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
            SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
        }
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
        for (const vector<std::array<double,3>>& points : polygons) {
            vector<std::array<double,2>> points_2D;
            for (const std::array<double,3>& point : points) {
                points_2D.push_back(camera.project(point));
                points_2D.back()[0] += SCREEN_WIDTH/2;
                points_2D.back()[1] = points_2D.back()[1]*-1+SCREEN_HEIGHT/2; 
            }
            fill_polygon(gRenderer, points_2D);
        }
        SDL_SetRenderDrawColor( gRenderer, 0x80, 0x80, 0x80, 0xFF );
        std::array<double,3> min_select;
        std::array<double,3> max_select;
        for (int i = 0; i < 3; i++) {
            min_select[i] = min(select[i],select[i]+select_size[i]);
            max_select[i] = max(select[i],select[i]+select_size[i]);
        }
        /*double meter = 1;
        for (const double& x : meters) {
            meter *= x;
        }
        for (int d1 = 0; d1 < 3; d1++) {
            for (const int& mult : {1,-1}) {
                int d2;
                int d3;
                if (d1 == 0) {
                    d2 = 1;
                    d3 = 2;
                }
                if (d1 == 1) {
                    d2 = 0;
                    d3 = 2;
                }
                if (d1 == 2) {
                    d2 = 0;
                    d3 = 1;
                }
                vector<vector<std::array<double,3>>> points;
                std::array<double,3> select1;
                std::array<double,3> select2;
                if (mult == 1) {
                    select1 = max_select;
                    select2 = min_select;
                } else {
                    select1 = min_select;
                    select2 = max_select;
                }
                for (int i = 0; i < (int)abs(select_size[d2]/meter); i++) {
                    points.push_back(vector<std::array<double,3>>());
                    for (int j = 0; j < (int)abs(select_size[d3]/meter); j++) {
                        std::array<double,3> point;
                        point[d2] = select2[d2] + mult*meter*(i+.5);
                        point[d3] = select2[d3] + mult*meter*(j+.5);
                        point[d1] = mult*std::numeric_limits<double>::infinity();
                        points[i].push_back(point);
                    }
                }
                for (int face_index = 0; face_index < poly.faces.size(); face_index++) {
                    vector<std::array<double,3>>* circuit = Polyhedron::circuit_cut(Polyhedron::make_clockwise(poly.circuits(face_index)));
                    if ((*circuit)[0][d1] == (*circuit)[1][d1] && (*circuit)[1][d1] == (*circuit)[2][d1]) {
                        for (vector<std::array<double,3>>& point_row : points) {
                            for (std::array<double,3>& point : point_row) {
                                std::array<double,3> projection;
                                projection[d2] = point[d2];
                                projection[d3] = point[d3];
                                projection[d1] = (*circuit)[0][d1];
                                bool any = false;
                                for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*circuit)) {
                                    if (Polyhedron::inside_triangle(triangle,projection)) {
                                        any = true;
                                        break;
                                    }
                                }
                                if (mult*(*circuit)[0][d1] < mult*point[d1] && mult*(*circuit)[0][d1] >= mult*select1[d1] && any) {
                                    point[d1] = (*circuit)[0][d1];
                                }
                            }
                        }
                    }
                }
                set<std::array<int,2>> seen;
                vector<set<std::array<int,2>>> components;
                for (int i = 0; i < points.size(); i++) {
                    for (int j = 0; j < points[i].size(); j++) {
                        if (seen.find({i,j}) != seen.end()) {
                            continue;
                        }
                        if (points[i][j][d1] == mult*std::numeric_limits<double>::infinity()) {
                            continue;
                        }
                        set<std::array<int,2>> component = {{i,j}};
                        seen.insert({i,j});
                        queue<std::array<int,2>> q;
                        q.push({i-1,j});
                        q.push({i,j-1});
                        q.push({i+1,j});
                        q.push({i,j+1});
                        while (!q.empty()) {
                            if (seen.find(q.front()) != seen.end()) {
                                q.pop();
                                continue;
                            }
                            if (q.front()[0] >= 0 && q.front()[0] < points.size() && q.front()[1] >= 0 && q.front()[1] < points[q.front()[0]].size() && points[i][j][d1] == points[q.front()[0]][q.front()[1]][d1]) {
                                component.insert(q.front());
                                seen.insert(q.front());
                                q.push({q.front()[0]-1,q.front()[1]});
                                q.push({q.front()[0],q.front()[1]-1});
                                q.push({q.front()[0]+1,q.front()[1]});
                                q.push({q.front()[0],q.front()[1]+1});
                            }
                            q.pop();
                        }
                        components.push_back(component);
                    }
                }
                for (const set<std::array<int,2>>& component : components) {
                    set<std::array<double,3>> point_set;
                    for (int i = 0; i < points.size(); i++) {
                        for (int j = 0; j < points[i].size(); j++) {
                            if (component.find({i,j}) != component.end()) {
                                if (component.find({i-1,j}) == component.end() || component.find({i,j-1}) == component.end() || component.find({i-1,j-1}) == component.end()) {
                                    std::array<double,3> point = points[i][j];
                                    point[d2] = select2[d2] + mult*i*meter;
                                    point[d3] = select2[d3] + mult*j*meter;
                                    point_set.insert(point);
                                }
                                if (component.find({i-1,j}) == component.end() || component.find({i,j+1}) == component.end() || component.find({i-1,j+1}) == component.end()) {
                                    std::array<double,3> point = points[i][j];
                                    point[d2] = select2[d2] + mult*i*meter;
                                    point[d3] = select2[d3] + mult*(j+1)*meter;
                                    point_set.insert(point);
                                }
                                if (component.find({i+1,j}) == component.end() || component.find({i,j-1}) == component.end() || component.find({i+1,j-1}) == component.end()) {
                                    std::array<double,3> point = points[i][j];
                                    point[d2] = select2[d2] + mult*(i+1)*meter;
                                    point[d3] = select2[d3] + mult*j*meter;
                                    point_set.insert(point);
                                }
                                if (component.find({i+1,j}) == component.end() || component.find({i,j+1}) == component.end() || component.find({i+1,j+1}) == component.end()) {
                                    std::array<double,3> point = points[i][j];
                                    point[d2] = select2[d2] + mult*(i+1)*meter;
                                    point[d3] = select2[d3] + mult*(j+1)*meter;
                                    point_set.insert(point);
                                }
                            }
                        }
                    }
                    std::array<double,3> point = {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
                    for (const std::array<double,3>& x : point_set) {
                        if (std::array<double,2>{x[d2],x[d3]} < std::array<double,2>{point[d2],point[d3]}) {
                            point = x;
                        }
                    }
                    vector<std::array<double,3>> path;
                    map<std::array<double,3>,std::array<double,3>> rotated_point_mapping;
                    for (const std::array<double,3>& x : point_set) {
                        rotated_point_mapping[x] = x;
                    }
                    map<std::array<double,3>,std::array<double,3>> reverse_point_mapping;
                    for (const pair<std::array<double,3>,std::array<double,3>>& p : rotated_point_mapping) {
                        reverse_point_mapping[p.second] = p.first;
                    }
                    std::array<double,3> angle = {0,0,0};
                    while (!path.size() || point != path[0]) {
                        path.push_back(point);
                        vector<std::array<double,3>> keys;
                        for (const pair<std::array<double,3>,std::array<double,3>>& p : rotated_point_mapping) {       
                            keys.push_back(p.first);
                        }
                        map<std::array<double,3>,std::array<double,3>> temp;
                        vector<std::array<double,3>>temp_keys(keys.begin(),keys.end());
                        for (std::array<double,3>& x: temp_keys) {
                            x[0] -= reverse_point_mapping[path.back()][0];
                            x[1] -= reverse_point_mapping[path.back()][1];
                            x[2] -= reverse_point_mapping[path.back()][2];
                        }
                        if (path.size() > 1) {
                            angle[d1] += atan2(reverse_point_mapping[path.back()][d3],reverse_point_mapping[path.back()][d2]);
                            temp_keys = Polyhedron::rotate(temp_keys, angle);
                        }
                        for (int i = 0; i < keys.size(); i++) {
                            temp[temp_keys[i]] = rotated_point_mapping[keys[i]];
                        }
                        rotated_point_mapping = temp;
                        reverse_point_mapping.clear();
                        for (const pair<std::array<double,3>,std::array<double,3>>& p : rotated_point_mapping) {
                            reverse_point_mapping[p.second] = p.first;
                        }
                        vector<std::array<double,2>> gift_wrap_angle_distance;
                        vector<std::array<double,3>> gift_wrap;
                        for (const pair<std::array<double,3>,std::array<double,3>>& p : rotated_point_mapping) {
                            std::array<double,3> origin = {0,0,0};
                            if (p.first == origin) {
                                continue;
                            }
                            gift_wrap_angle_distance.push_back({atan2(p.first[d3],p.first[d2]),Polyhedron::distance(p.first,origin)});
                            gift_wrap.push_back(rotated_point_mapping[p.first]);
                        }
                        std::array<double,2> maxi_angle_distance = {-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
                        std::array<double,3> maxi;
                        for (int i = 0; i < gift_wrap_angle_distance.size(); i++) {
                            if (gift_wrap_angle_distance[i] > maxi_angle_distance) {
                                maxi_angle_distance = gift_wrap_angle_distance[i];
                                maxi = gift_wrap[i];
                            }
                        }
                        point = maxi;
                    }
                    for (int i = 0; i < path.size();) {
                        if (Polyhedron::colinear(set<std::array<double,3>>{path[(i-1+path.size())%path.size()],path[i],path[(i+1)%path.size()]})) {
                            path.erase(path.begin()+i);
                        } else {
                            i++;
                        }
                    }
                    vector<std::array<double,2>> path_2D;
                    for (const std::array<double,3>& x : path) {
                        path_2D.push_back(camera.project(x));
                        path_2D.back()[0] += SCREEN_WIDTH/2;
                        path_2D.back()[1] = path_2D.back()[1]*-1+SCREEN_HEIGHT/2;
                    }
                    fill_polygon(gRenderer, path_2D);
                }
            }
        }*/
        SDL_SetRenderDrawColor( gRenderer, 0x00, 0x00, 0x00, 0xFF );
        for (const set<int>& edge : poly.edges) {
            std::array<double,2> p1 = camera.project(poly.verts[*(edge.begin())]);
            std::array<double,2> p2 = camera.project(poly.verts[*(++edge.begin())]);
            p1[0] += SCREEN_WIDTH/2;
            p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
            p2[0] += SCREEN_WIDTH/2;
            p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
            SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
        }
        double width = size[0];
        double height = size[1];
        double depth = size[2];
        Polyhedron axes;
        if (select_dimension == 0) {
            axes = get_cube({origin[0]+width/2,origin[1]+abs(select_size[1])/2+min(select[1],select[1]+select_size[1]),origin[2]+depth/2}, {width,abs(select_size[1]),depth});
        } else if (select_dimension == 1) {
            axes = get_cube({origin[0]+width/2,origin[1]+height/2,origin[2]+abs(select_size[2])/2+min(select[2],select[2]+select_size[2])}, {width,height,abs(select_size[2])});
        } else if (select_dimension == 2) {
            axes = get_cube({origin[0]+abs(select_size[0])/2+min(select[0],select[0]+select_size[0]),origin[1]+height/2,origin[2]+depth/2}, {abs(select_size[0]),height,depth});
        }
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
        for (const set<int>& edge : axes.edges) {
            std::array<double,2> p1 = camera.project(axes.verts[*(edge.begin())]);
            std::array<double,2> p2 = camera.project(axes.verts[*(++edge.begin())]);
            p1[0] += SCREEN_WIDTH/2;
            p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
            p2[0] += SCREEN_WIDTH/2;
            p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
            SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
        }
        SDL_SetRenderDrawColor( gRenderer, 0x80, 0x80, 0x80, 0xFF );
        if (select[2] > origin[2]) {
            Polyhedron shadow = get_cube({select[0]+select_size[0]/2,origin[1]+select[1]+select_size[1]/2,(select[2]-origin[2])/2}, {select_size[0],select_size[1],select[2]-origin[2]});
            for (const set<int>& edge : shadow.edges) {
                std::array<double,2> p1 = camera.project(shadow.verts[*(edge.begin())]);
                std::array<double,2> p2 = camera.project(shadow.verts[*(++edge.begin())]);
                p1[0] += SCREEN_WIDTH/2;
                p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
                p2[0] += SCREEN_WIDTH/2;
                p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
                SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
            }
        }
        if (select[1] > origin[1]) {
            Polyhedron shadow = get_cube({select[0]+select_size[0]/2,(select[1]-origin[1])/2,select[2]+select_size[2]/2}, {select_size[0],select[1]-origin[1],select_size[2]});
            for (const set<int>& edge : shadow.edges) {
                std::array<double,2> p1 = camera.project(shadow.verts[*(edge.begin())]);
                std::array<double,2> p2 = camera.project(shadow.verts[*(++edge.begin())]);
                p1[0] += SCREEN_WIDTH/2;
                p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
                p2[0] += SCREEN_WIDTH/2;
                p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
                SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
            }
        }
        if (select[0] > origin[0]) {
            Polyhedron shadow = get_cube({(select[0]-origin[0])/2,select[1]+select_size[1]/2,select[2]+select_size[2]/2}, {select[0]-origin[0],select_size[1],select_size[2]});
            for (const set<int>& edge : shadow.edges) {
                std::array<double,2> p1 = camera.project(shadow.verts[*(edge.begin())]);
                std::array<double,2> p2 = camera.project(shadow.verts[*(++edge.begin())]);
                p1[0] += SCREEN_WIDTH/2;
                p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
                p2[0] += SCREEN_WIDTH/2;
                p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
                SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
            }
        }
        Polyhedron select_cube = get_cube({select[0]+select_size[0]/2,select[1]+select_size[1]/2,select[2]+select_size[2]/2}, select_size);
        SDL_SetRenderDrawColor( gRenderer, 0xA0, 0xA0, 0xA0, 0xFF );
        for (const set<int>& edge : select_cube.edges) {
            std::array<double,2> p1 = camera.project(select_cube.verts[*(edge.begin())]);
            std::array<double,2> p2 = camera.project(select_cube.verts[*(++edge.begin())]);
            p1[0] += SCREEN_WIDTH/2;
            p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
            p2[0] += SCREEN_WIDTH/2;
            p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
            SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
        }
        if (cylinder_tool) {
            if (cylinder_dimension == 0) {
                std::array<std::array<double,2>,2> left_edge = {{{std::numeric_limits<double>::infinity(),0},{std::numeric_limits<double>::infinity(),0}}};
                std::array<std::array<double,2>,2> right_edge = {{{-std::numeric_limits<double>::infinity(),0},{-std::numeric_limits<double>::infinity(),0}}};
                std::array<double,3> center = {select[0]+select_size[0]/2, select[1], select[2]+select_size[2]/2};
                int point_count = 100;
                vector<std::array<double,3>> ellipsoid;
                for(int i = 0; i < point_count; i++) {
                    ellipsoid.push_back({center[0]+cos(i*2*pi/point_count)*select_size[0]/2,center[1],center[2]+sin(i*2*pi/point_count)*select_size[2]/2});
                }
                vector<std::array<double,2>> ellipsoid_2D;
                for (const std::array<double,3>& x : ellipsoid) {
                    ellipsoid_2D.push_back(camera.project(x));
                    ellipsoid_2D.back()[0] += SCREEN_WIDTH/2;
                    ellipsoid_2D.back()[1] = ellipsoid_2D.back()[1]*-1+SCREEN_HEIGHT/2;
                    if (ellipsoid_2D.back()[0] < left_edge[0][0]) {
                        left_edge[0] = ellipsoid_2D.back();
                    }
                    if (ellipsoid_2D.back()[0] > right_edge[0][0]) {
                        right_edge[0] = ellipsoid_2D.back();
                    }
                }
                for (int i = 0; i < ellipsoid_2D.size(); i ++) {
                    SDL_RenderDrawLine(gRenderer, ellipsoid_2D[i][0], ellipsoid_2D[i][1], ellipsoid_2D[(i+1)%ellipsoid_2D.size()][0], ellipsoid_2D[(i+1)%ellipsoid_2D.size()][1]);
                }
                center = {select[0]+select_size[0]/2, select[1]+select_size[1], select[2]+select_size[2]/2};
                ellipsoid.clear();
                for(int i = 0; i < point_count; i++) {
                    ellipsoid.push_back({center[0]+cos(i*2*pi/point_count)*select_size[0]/2,center[1],center[2]+sin(i*2*pi/point_count)*select_size[2]/2});
                }
                ellipsoid_2D.clear();
                for (const std::array<double,3>& x : ellipsoid) {
                    ellipsoid_2D.push_back(camera.project(x));
                    ellipsoid_2D.back()[0] += SCREEN_WIDTH/2;
                    ellipsoid_2D.back()[1] = ellipsoid_2D.back()[1]*-1+SCREEN_HEIGHT/2;
                    if (ellipsoid_2D.back()[0] < left_edge[1][0]) {
                        left_edge[1] = ellipsoid_2D.back();
                    }
                    if (ellipsoid_2D.back()[0] > right_edge[1][0]) {
                        right_edge[1] = ellipsoid_2D.back();
                    }
                }
                for (int i = 0; i < ellipsoid_2D.size(); i ++) {
                    SDL_RenderDrawLine(gRenderer, ellipsoid_2D[i][0], ellipsoid_2D[i][1], ellipsoid_2D[(i+1)%ellipsoid_2D.size()][0], ellipsoid_2D[(i+1)%ellipsoid_2D.size()][1]);
                }
                SDL_RenderDrawLine(gRenderer, left_edge[0][0], left_edge[0][1], left_edge[1][0], left_edge[1][1]);
                SDL_RenderDrawLine(gRenderer, right_edge[0][0], right_edge[0][1], right_edge[1][0], right_edge[1][1]);
            }
            set<set<std::array<double,3>>> ellipsoid_shadow;
            std::array<double,3> center = {select[0]+select_size[0]/2, min_select[1], select[2]+select_size[2]/2};
            for(double i = min_select[0]; i < max_select[0]; i += unit) {
                for(double j = min_select[2]; j < max_select[2]; j += unit) {
                    bool all = true;
                    for (const double& k : {i,i+unit}) {
                        for (const double& l : {j,j+unit}) {
                            double angle = atan2(l-center[2], k-center[0]);
                            if (Polyhedron::distance(center,{k,center[1],l}) > sqrt(pow(select_size[0]/2*cos(angle),2)+pow(select_size[2]/2*sin(angle),2))) {
                                all = false;
                                break;
                            }
                        }
                        if (!all) {
                            break;
                        }
                    }
		    if (!all) {
		    	continue;
		    }
                    set<set<std::array<double,3>>> edges;
                    edges.insert({{i,center[1],j},{i+unit,center[1],j}});
                    edges.insert({{i+unit,center[1],j},{i+unit,center[1],j+unit}});
                    edges.insert({{i+unit,center[1],j+unit},{i,center[1],j+unit}});
                    edges.insert({{i,center[1],j+unit},{i,center[1],j}});
                    for (const set<std::array<double,3>>& edge : edges) {
                        if (ellipsoid_shadow.find(edge) == ellipsoid_shadow.end()) {
                            ellipsoid_shadow.insert(edge);
                        } else {
                            ellipsoid_shadow.erase(edge);
                        }
                    }
                }
            }
            map<double,set<set<std::array<double,3>>>> map;
            set<set<std::array<double,3>>> ellipsoid_top;
            for (const set<std::array<double,3>>& edge : ellipsoid_shadow) {
                set<std::array<double,3>>::iterator it = edge.begin();
                map[(*it)[0]].insert(edge);
                map[(*++it)[0]].insert(edge);
                it = edge.begin();
                set<std::array<double,3>> new_edge = {{(*it)[0],max_select[1],(*it)[2]},{(*++it)[0],max_select[1],(*it)[2]}};
                ellipsoid_top.insert(new_edge);
                //std::array<double,2> p1 = camera.project(*it);
                //std::array<double,2> p2 = camera.project(*(++it));
                //p1[0] += SCREEN_WIDTH/2;
                //p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
                //p2[0] += SCREEN_WIDTH/2;
                //p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
                //SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
            }
            vector<set<set<std::array<double,3>>>> faces;
            faces.push_back(ellipsoid_shadow);
            faces.push_back(ellipsoid_top);
            cylinder_boxes.clear();
            double i;
            for(i = min_select[0]; i <= max_select[0]; i += unit) {
                if (map[i].size() == 2) {
                    std::array<double,3> p1;
                    std::array<double,3> p4;
                    set<set<std::array<double,3>>>::iterator it1 = map[i].begin();
                    set<std::array<double,3>>::iterator it2 = it1->begin();
                    if ((*it2)[0] == i) {
                        p1 = (*it2);
                    } else {
                        p1 = *(++it2);
                    }
                    it2 = (++it1)->begin();
                    if ((*it2)[0] == i) {
                        p4 = (*it2);
                    } else {
                        p4 = *(++it2);
                    }
                    vector<std::array<double,3>> circuit = {p1,{p1[0],max_select[1],p1[2]},{p4[0],max_select[1],p4[2]},p4};
                    set<set<std::array<double,3>>> face;
                    for (int j = 0; j < circuit.size(); j++) {
                        face.insert({circuit[(j-1+circuit.size())%circuit.size()],circuit[j]});
                    }
                    faces.push_back(face);
                }
                if (i == max_select[0]) {
                    break;
                }
                if (map[i].size() && map[i+unit].size()) {
                    set<set<std::array<double,3>>> edge_intersect;
                    set_intersection(map[i].begin(),map[i].end(),map[i+unit].begin(),map[i+unit].end(),inserter(edge_intersect, edge_intersect.begin()));
                    set<set<std::array<double,3>>>::iterator it1 = edge_intersect.begin();
                    set<std::array<double,3>>::iterator it2 = (*it1).begin();
                    std::array<std::array<double,3>,2> edge1_array = {(*it2),*(++it2)};
                    it2 = (*++it1).begin();
                    std::array<std::array<double,3>,2> edge2_array = {(*it2),*(++it2)};
                    Box::PointDistanceComparator comp = Box::PointDistanceComparator(*it2);
                    sort(edge1_array.begin(),edge1_array.end(), comp);
                    sort(edge2_array.begin(),edge2_array.end(), comp);
                    std::array<double,2> x = {std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
                    std::array<double,2> y = {min_select[1],max_select[1]};
                    std::array<double,2> z = {std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
                    for (const std::array<double,3>& point : edge1_array) {
                        if (point[0] < x[0]) {
                            x[0] = point[0];
                        }
                        if (point[0] > x[1]) {
                            x[1] = point[0];
                        }
                        if (point[2] < z[0]) {
                            z[0] = point[2];
                        }
                        if (point[2] > z[1]) {
                            z[1] = point[2];
                        }
                    }
                    for (const std::array<double,3>& point : edge2_array) {
                        if (point[0] < x[0]) {
                            x[0] = point[0];
                        }
                        if (point[0] > x[1]) {
                            x[1] = point[0];
                        }
                        if (point[2] < z[0]) {
                            z[0] = point[2];
                        }
                        if (point[2] > z[1]) {
                            z[1] = point[2];
                        }
                    }
                    if (cylinder_boxes.size() && cylinder_boxes.back().z_min == z[0] && cylinder_boxes.back().z_max == z[1]) {
                        cylinder_boxes.back() = Box({cylinder_boxes.back().x_min,x[1]},y,z);
                    } else {
                        cylinder_boxes.push_back(Box(x,y,z));
                    }
                    set<set<std::array<double,3>>> face;
                    face.insert({edge1_array[0],edge1_array[1]});
                    face.insert({edge1_array[1],{edge1_array[1][0],max_select[1],edge1_array[1][2]}});
                    face.insert({{edge1_array[1][0],max_select[1],edge1_array[1][2]},{edge1_array[0][0],max_select[1],edge1_array[0][2]}});
                    face.insert({{edge1_array[0][0],max_select[1],edge1_array[0][2]},edge1_array[0]});
                    faces.push_back(face);
                    face.clear();
                    face.insert({edge2_array[0],edge2_array[1]});
                    face.insert({edge2_array[1],{edge2_array[1][0],max_select[1],edge2_array[1][2]}});
                    face.insert({{edge2_array[1][0],max_select[1],edge2_array[1][2]},{edge2_array[0][0],max_select[1],edge2_array[0][2]}});
                    face.insert({{edge2_array[0][0],max_select[1],edge2_array[0][2]},edge2_array[0]});
                    faces.push_back(face);
                }
            }
            for (int i = 0; i < faces.size(); i++) {
                for (int j = i+1; j < faces.size();) {
                    set<set<std::array<double,3>>> face_intersect;
                    set_intersection(faces[i].begin(),faces[i].end(),faces[j].begin(),faces[j].end(),inserter(face_intersect, face_intersect.begin()));
                    set<std::array<double,3>> point_set;
                    for (const set<std::array<double,3>>& edge : faces[i]) {
                        for (const std::array<double,3>& point : edge) {
                            point_set.insert(point);
                        }
                    } 
                    for (const set<std::array<double,3>>& edge : faces[j]) {
                        for (const std::array<double,3>& point : edge) {
                            point_set.insert(point);
                        }
                    } 
                    if (face_intersect.size() && Polyhedron::coplanar(point_set)) {
                        faces[i].insert(faces[j].begin(),faces[j].end());
                        for (const set<std::array<double,3>>& edge : face_intersect) {
                            faces[i].erase(edge);
                        }
                        faces.erase(faces.begin()+j);
                    } else {
                        j++;
                    }
                }
            }
            set<std::array<double,3>> verts_set;
            set<set<std::array<double,3>>> edges_set;
            for (const set<set<std::array<double,3>>>& face : faces) {
                for (const set<std::array<double,3>>& edge : face) {
                    for (const std::array<double,3>& point : edge) {
                        verts_set.insert(point);
                    }
                    edges_set.insert(edge);
                }
            }
            cylinder_poly.verts.clear();
            cylinder_poly.verts.insert(cylinder_poly.verts.end(),verts_set.begin(),verts_set.end());
            vector<set<std::array<double,3>>> edges_vector(edges_set.begin(),edges_set.end());
            cylinder_poly.faces.clear();
            for (const set<set<std::array<double,3>>>& face : faces) {
                set<int> face_ungrounded;
                for (const set<std::array<double,3>>& edge : face) {
                    vector<set<std::array<double,3>>>::iterator it = find(edges_vector.begin(),edges_vector.end(),edge);
                    face_ungrounded.insert(std::distance(edges_vector.begin(), it));
                }
                cylinder_poly.faces.push_back(face_ungrounded);
            }
            cylinder_poly.edges.clear();
            for (const set<std::array<double,3>>& edge : edges_vector) {
                set<int> edge_ungrounded;
                for (const std::array<double,3>& point : edge) {
                    vector<std::array<double,3>>::iterator it = find(cylinder_poly.verts.begin(),cylinder_poly.verts.end(),point);
                    edge_ungrounded.insert(std::distance(cylinder_poly.verts.begin(),it));
                }
                cylinder_poly.edges.push_back(edge_ungrounded);
            }
            for (const set<int>& edge : cylinder_poly.edges) {
                set<int>::iterator it = edge.begin();
                std::array<double,2> p1 = camera.project(cylinder_poly.verts[*(it)]);
                std::array<double,2> p2 = camera.project(cylinder_poly.verts[*(++it)]);
                p1[0] += SCREEN_WIDTH/2;
                p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
                p2[0] += SCREEN_WIDTH/2;
                p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
                SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
            }
        }
    }
    voxel_editor::Block reunit(double new_unit) {
        voxel_editor::Block b({origin[0],origin[0]+size[0]}, {origin[1],origin[1]+size[1]},{origin[2],origin[2]+size[2]}, new_unit);
        b.select = select;
        b.select[0] = floor(b.select[0]/b.unit)*b.unit;
        b.select[1] = floor(b.select[1]/b.unit)*b.unit;
        b.select[2] = floor(b.select[2]/b.unit)*b.unit;
	    b.select_dimension = select_dimension;
	    b.select_size = select_size;
	    b.space_down = space_down;
	    if (!space_down || select_dimension == 2) {
       	    b.select_size[0] = floor(b.select_size[0]/unit)*b.unit;
        }
	    if (!space_down || select_dimension == 0) {
        	b.select_size[1] = floor(b.select_size[1]/unit)*b.unit;
	    }
        if (!space_down || select_dimension == 1) {
            b.select_size[2] = floor(b.select_size[2]/unit)*b.unit;
        }
        b.poly = poly;
        b.meters = meters;
        b.polygons = polygons;
    	b.cylinder_tool = cylinder_tool;
    	b.cylinder_dimension = cylinder_dimension;
    	b.cylinder_options = cylinder_options;
        return b;
    }
    voxel_editor::Block resize(std::array<double,2>x, std::array<double,2>y, std::array<double,2>z) {
        voxel_editor::Block b(x, y, z, unit);
        b.select = select;
        b.poly = poly;
        b.meters = meters;
        b.polygons = polygons;
        return b;
    }
    bool select_by_void() {
        return true;
    }
    bool dump(string filename) {
        string header_insert;
        header_insert += "element block 1\n";
        header_insert += "property float x_min\n";
        header_insert += "property float x_max\n";
        header_insert += "property float y_min\n";
        header_insert += "property float y_max\n";
        header_insert += "property float z_min\n";
        header_insert += "property float z_max\n";
        header_insert += "property float unit\n";
        header_insert += "property list uchar float meters\n";
        string output = this->poly.dump();
        string after = "format ascii 1.0\n";
        int index = output.find(after)+after.size();
        output = output.substr(0,index) + header_insert + output.substr(index,output.size()-index);
        string body_insert;
        body_insert += to_string(origin[0]) + " ";
        body_insert += to_string(origin[0]+size[0]) + " ";
        body_insert += to_string(origin[1]) + " ";
        body_insert += to_string(origin[1]+size[1]) + " ";
        body_insert += to_string(origin[2]) + " ";
        body_insert += to_string(origin[2]+size[2]) + " ";
        body_insert += to_string(unit) + " ";
        body_insert += to_string(meters.size()) + " ";
        for (const double& meter : meters) {
            body_insert += to_string(meter) + " ";
        }
        body_insert += "\n";
        after = "end_header\n";
        index = output.find(after)+after.size();
        output = output.substr(0,index) + body_insert + output.substr(index,output.size()-index);
        ofstream output_file(filename);
        if (output_file.is_open()) {
            output_file << output;
            return true;
        }
        return false;
    }
    void loads(string text) {
        int element_index = -1;
        vector<string> elements_vector;
        map<string,vector<string>> elements;
        string poly_text;
        bool block_element = false;
        bool first_body_line = true;
        while (text.size()) {
            std::string line;
            if (text.find("\n") != std::string::npos) {
                line = text.substr(0,text.find("\n")+1);
                text = text.substr(text.find("\n")+1,text.size()-1-1);
            } else {
                line = text;
                text = "";
            }
            string trim_line = trim(line);
            vector<string> split_line;
            while (trim_line.size()) {
                string token;
                if (trim_line.find(" ") != std::string::npos) {
                    token = trim_line.substr(0, trim_line.find(" "));
                    trim_line = trim_line.substr(trim_line.find(" ")+1,trim_line.size()-trim_line.find(" ")-1);
                } else {
                    token = trim_line;
                    trim_line = "";
                }
                if (token.size()) {
                    split_line.push_back(token);
                }
            }
            if (split_line.size() >= 1 && split_line[0] == "end_header") {
                element_index = 0;
                poly_text += line;
            } else if (element_index == -1) {
                if (split_line.size() >= 3) {
                    if (split_line[0] == "element") {
                        elements_vector.push_back(split_line[1]);
                        int element_count = stoi(split_line[2]);
                        block_element = split_line[1] == "block";
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
                if (!block_element) {
                    poly_text += line;
                }
            } else {
                if (first_body_line) {
                    first_body_line = false;
                } else {
                    poly_text += line;
                    continue;
                }
                string element = elements_vector[element_index];
                if (split_line.size() >= elements[element].size()) {
                    int index = 0;
                    std::array<double,3> vert;
                    for (const string& property : elements[element]) {
                        string meter_list = "meters";
                        if (property == meter_list) {
                            int property_count = stoi(split_line[index]);
                            for (int i = 0; i < property_count; i++) {
                                meters.insert(stod(split_line[index+1+i]));
                            }
                            index += 1+property_count;
                        } else {
                            if (property == "x_min") {
                                origin[0] = stod(split_line[index]);
                            } else if (property == "x_max") {
                                size[0] = stod(split_line[index])-origin[0];
                            } else if (property == "y_min") {
                                origin[1] = stod(split_line[index]);
                            } else if (property == "y_max") {
                                size[1] = stod(split_line[index])-origin[1];
                            } else if (property == "z_min") {
                                origin[2] = stod(split_line[index]);
                            } else if (property == "z_max") {
                                size[2] = stod(split_line[index])-origin[2];
                            } else if (property == "unit") {
                                unit = stod(split_line[index]);
                            }
                            index++;
                        }
                    }
                }
            }
        }
        poly.loads(poly_text); 
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
};
}

bool init();

//Frees media and shuts down SDL
void close();

//Loads individual image
SDL_Surface* loadSurface( std::string path );

//The window we'll be rendering to
SDL_Window* gWindow = NULL;

//The surface contained by the window
SDL_Surface* gScreenSurface = NULL;

//The window renderer
SDL_Renderer* gRenderer = NULL;

bool init() {
    // Initialization flag
    bool success = true;

    // Initialize SDL
    if( SDL_Init( SDL_INIT_VIDEO ) < 0 )
    {
        printf( "SDL could not initialize! SDL Error: %s\n", SDL_GetError() );
        success = false;
    }
    else
    {
        // Create window
        gWindow = SDL_CreateWindow( "SDL Tutorial", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
        if( gWindow == NULL )
        {
            printf( "Window could not be created! SDL Error: %s\n", SDL_GetError() );
            success = false;
        }
        else
        {
            //Create renderer for window
			gRenderer = SDL_CreateRenderer( gWindow, -1, SDL_RENDERER_ACCELERATED );
			if( gRenderer == NULL )
			{
				printf( "Renderer could not be created! SDL Error: %s\n", SDL_GetError() );
				success = false;
			}
			else
			{
				//Initialize renderer color
				SDL_SetRenderDrawColor( gRenderer, 0x00, 0x00, 0x00, 0xFF );

            }

        }
    }

    return success;
}

int main( int argc, char* args[] ) {
    if( !init() ) {
        printf( "Failed to initialize!\n" );
    } else {
        camera.zoom = 100;
        camera.focal[0] = 0;
        camera.focal[1] = 0;
        camera.focal[2] = -100;
        Polyhedron cube = get_cube({0,0,-1},{1,1,1},{0,0,0});
        Crosshair crosshair;
        bool mouse_down = false;
        int mousedown_x;
        int mousedown_y;
        SDL_ShowCursor(SDL_DISABLE);
        bool del = false;
        bool meta_down = false;
        bool shift_down = true;
        voxel_editor::Block block(3,3,3,1);
        // Main loop flag
        bool quit = false;
        // Event handler
        SDL_Event e;
        //block.add_poly({0,0,0},{1,1,1},true);
        // While application is running
        while( !quit ) {
            auto start = std::chrono::high_resolution_clock::now();
            // Handle events on queue
            while( SDL_PollEvent( &e ) != 0 ) {
                 std::array<SDL_KeyCode,10> number_keys = {SDLK_1,SDLK_2,SDLK_3,SDLK_4,SDLK_5,SDLK_6,SDLK_7,SDLK_8,SDLK_9,SDLK_0};    
                //User requests quit
                if( e.type == SDL_QUIT ) {
                    quit = true;
                } else if( e.type == SDL_MOUSEMOTION) {
                    int x, y;
                    SDL_GetMouseState( &x, &y );
                    crosshair.pos[0] = x-SCREEN_WIDTH/2;
                    crosshair.pos[1] = -(y-SCREEN_HEIGHT/2);
                    if (mouse_down) {
                        double delta_x = x-mousedown_x;
                        double delta_y = y-mousedown_y;
                        camera.move({camera.x_vector[0]*-delta_x/camera.zoom,camera.x_vector[1]*-delta_x/camera.zoom,camera.x_vector[2]*-delta_x/camera.zoom});
                        camera.move({camera.y_vector[0]*delta_y/camera.zoom,camera.y_vector[1]*delta_y/camera.zoom,camera.y_vector[2]*delta_y/camera.zoom});
                        mousedown_x = x;
                        mousedown_y = y;
                    }
                } else if (e.type == SDL_MOUSEBUTTONDOWN) {
                    SDL_GetMouseState( &mousedown_x, &mousedown_y );
                    mouse_down = true;
                } else if (e.type == SDL_MOUSEBUTTONUP) {
                    mouse_down = false;
                } else if(e.type == SDL_KEYDOWN && find(number_keys.begin(),number_keys.end(),e.key.keysym.sym) != number_keys.end()) {
                    if (e.key.keysym.sym == SDLK_1) {
                        block = block.reunit(1);
                    } else {
                        int index = distance(number_keys.begin(), find(number_keys.begin(),number_keys.end(),e.key.keysym.sym));
                        if (!meta_down) {
                            block = block.reunit(block.unit/(index+1));
                            block.meters = Polyhedron::remeter(block.meters, block.unit);
                        } else {
                            block = block.reunit(block.unit*(index+1));
                        }
                    }
                } else if(e.type == SDL_KEYDOWN) {
                    int direction = 1;
                    switch( e.key.keysym.sym ) {
                        case SDLK_a:
                            camera.rotate({0,M_PI/180*10,0});
                            break;

                        case SDLK_d:
                            camera.rotate({0,-M_PI/180*10,0});
                            break;

                        case SDLK_w:
                            camera.rotate({M_PI/180*10,0,0});
                            break;

                        case SDLK_s:
                            if (meta_down) {
                                if (shift_down) {
                                    cout << "(Exporting...) Enter filename: ";
                                    string filename;
                                    cin >> filename;
                                    block.poly.exp(filename);
                                } else {
                                    cout << "(Saving...) Enter filename: ";
                                    string filename;
                                    cin >> filename;
                                    block.dump(filename);
                                }
                            } else {
                                camera.rotate({-M_PI/180*10,0,0});
                            }
                            break;

                        case SDLK_q:
                            camera.rotate({0,0,-M_PI/180*10});
                            break;

                        case SDLK_e:
                            camera.rotate({0,0,M_PI/180*10});
                            break;
                        case SDLK_o:
                            if (meta_down) {
                                cout << "(Loading...) Enter filename: ";
                                string filename;
                                cin >> filename;
                                block.load(filename);
                                cout << block.poly.verts.size() << " " << block.poly.edges.size() << " " << block.poly.faces.size() << endl;
                                block.flip();
                            }
                            break; 
                        case SDLK_n:
                            if (meta_down) {
                                cout << "(Creating New...) Enter the origin (x y z): ";
                                string line;
                                std::getline(cin, line);
                                std::array<double,3> origin;
                                try {
                                    int count = 0;
                                    while (line.size() && count < 3) {
                                        string token;
                                        if (line.find(" ") != std::string::npos) {
                                            token = line.substr(0, line.find(" "));
                                            line = line.substr(line.find(" ")+1,line.size()-line.find(" ")-1);
                                        } else {
                                            token = line;
                                            line = "";
                                        }
                                        origin[count] = stod(token);
                                        count++;
                                    }
                                } catch (const std::invalid_argument& e) {
                                    origin = {0,0,0};
                                }
                                cout << "(Creating New...) Enter the size: (x y z): ";
                                line.clear();
                                std::getline(cin, line);
                                std::array<double,3> size;
                                try {
                                    int count = 0;
                                    while (line.size() && count < 3) {
                                        string token;
                                        if (line.find(" ") != std::string::npos) {
                                            token = line.substr(0, line.find(" "));
                                            line = line.substr(line.find(" ")+1,line.size()-line.find(" ")-1);
                                        } else {
                                            token = line;
                                            line = "";
                                        }
                                        size[count] = stod(token);
                                        count++;
                                    }
                                } catch (const std::invalid_argument& e) {
                                    size = {3,3,3};
                                }
                                block = voxel_editor::Block({origin[0],origin[0]+size[0]},{origin[1],origin[1]+size[1]},{origin[2],origin[2]+size[2]},1);
                            }
                            break; 
                        case SDLK_p: {
                            cout << "(Padding...) Enter padding dimension (left, right, bottom, top, front, back) and unit multiple: ";
                            string line;
                            string padding_dimension;
                            int padding;
                            std::getline(cin, line);
                            try {
                                int count = 0;
                                while (line.size() && count < 2) {
                                    string token;
                                    if (line.find(" ") != std::string::npos) {
                                        token = line.substr(0, line.find(" "));
                                        line = line.substr(line.find(" ")+1,line.size()-line.find(" ")-1);
                                    } else {
                                        token = line;
                                        line = "";
                                    }
                                    if (count == 0) {
                                       padding_dimension = token;
                                    } else {
                                       padding = stoi(token);
                                    }
                                    count++;
                                }
                            } catch (const std::invalid_argument& e) {
                                break;
                            }
                            if (padding_dimension == "left") {
                                if (padding < 0) {
                                    Box box({block.origin[0],block.origin[0]-padding*block.unit},{block.origin[1],block.origin[1]+block.size[1]},{block.origin[2],block.origin[2]+block.size[2]});
                                    block.poly = box.del(block.poly);
                                }
                                block = block.resize({block.origin[0]-padding*block.unit,block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]+block.size[1]},{block.origin[2],block.origin[2]+block.size[2]});
                                if (block.select[0] > block.origin[0]) {
                                    block.select[0] = block.origin[0];
                                }
                            } else if (padding_dimension == "right") {
                                if (padding < 0) {
                                    Box box({block.origin[0]+block.size[0]+padding*block.unit,block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]+block.size[1]},{block.origin[2],block.origin[2]+block.size[2]});
                                    block.poly = box.del(block.poly);
                                }
                                block = block.resize({block.origin[0],block.origin[0]+block.size[0]+padding*block.unit},{block.origin[1],block.origin[1]+block.size[1]},{block.origin[2],block.origin[2]+block.size[2]});
                                if (block.select[0] > block.origin[0]+block.size[0]-block.unit) {
                                    block.select[0] = block.origin[0]+block.size[0]-block.unit;
                                }
                            } else if (padding_dimension == "bottom") {
                                if (padding < 0) {
                                    Box box({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]-padding*block.unit},{block.origin[2],block.origin[2]+block.size[2]});
                                    block.poly = box.del(block.poly);
                                }
                                block = block.resize({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1]-padding*block.unit,block.origin[1]+block.size[1]},{block.origin[2],block.origin[2]+block.size[2]});
                                if (block.select[1] > block.origin[1]) {
                                    block.select[1] = block.origin[1];
                                }
                            } else if (padding_dimension == "top") {
                                if (padding < 0) {
                                    Box box({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1]+block.size[1]+padding*block.unit,block.origin[1]+block.size[1]},{block.origin[2],block.origin[2]+block.size[2]});
                                    block.poly = box.del(block.poly);
                                }
                                block = block.resize({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]+block.size[1]+padding*block.unit},{block.origin[2],block.origin[2]+block.size[2]});
                                if (block.select[1] > block.origin[1]+block.size[1]-block.unit) {
                                    block.select[1] = block.origin[1]+block.size[1]-block.unit;
                                }
                            } else if (padding_dimension == "front") {
                                if (padding < 0) {
                                    Box box({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]+block.size[1]},{block.origin[2],block.origin[2]-padding*block.unit});
                                    block.poly = box.del(block.poly);
                                }
                                block = block.resize({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]+block.size[1]},{block.origin[2]-padding*block.unit,block.origin[2]+block.size[2]});
                                if (block.select[2] > block.origin[2]) {
                                    block.select[2] = block.origin[2];
                                }
                            } else if (padding_dimension == "back") {
                                if (padding < 0) {
                                    Box box({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]+block.size[1]},{block.origin[2]+block.size[2]+padding*block.unit,block.origin[2]+block.size[2]});
                                    block.poly = box.del(block.poly);
                                }
                                block = block.resize({block.origin[0],block.origin[0]+block.size[0]},{block.origin[1],block.origin[1]+block.size[1]+padding*block.unit},{block.origin[2],block.origin[2]+block.size[2]+padding*block.unit});
                                if (block.select[2] > block.origin[2]+block.size[2]-block.unit) {
                                    block.select[2] = block.origin[2]+block.size[2]-block.unit;
                                }
                            }
                            if (padding < 0) {
                                Polyhedron::RoundPointMeter rounder(block.meters);
                                block.poly.round_verts(rounder);
                                block.flip();
                            }
                            break;
                        } case SDLK_c:
                            block.cylinder_tool = !block.cylinder_tool;
                            cout << "cylinder_tool " << block.cylinder_tool << endl;
                            break;
                        case SDLK_LSHIFT:
                            if (!meta_down) {
                                block.select_dimension = (block.select_dimension+3-1)%3;
                            }
                            shift_down = true;   
                            break;
                        
                        case SDLK_RSHIFT:
                            block.select_dimension = (block.select_dimension+1)%3;
                            break;

                        case SDLK_z:
                            camera.zoom *= 1.1;
                            break;

                        case SDLK_x:
                            camera.zoom /= 1.1;
                            break;

                        case SDLK_SPACE:
                            block.space_down = true;
                            break;
                        case SDLK_DOWN:
                        case SDLK_LEFT:
                            direction = -1;
                            break;
                        case SDLK_BACKSPACE:
                            del = !del;  
                            break;
                        case SDLK_LGUI:
                        case SDLK_RGUI:
                            meta_down = true;
                            break;
                        default:
                            break;
                    }
                    int dir_mult1 = 1;
                    if (Polyhedron::dot(camera.forward_vector(),{0,0,1}) < Polyhedron::dot(camera.forward_vector(),{0,0,-1})) {
                        dir_mult1 = -1;
                    }
                    int dir_mult2 = 1;
                    if (Polyhedron::dot(camera.forward_vector(),{1,0,0}) < Polyhedron::dot(camera.forward_vector(),{-1,0,0})) {
                        dir_mult2 = -1;
                    }
                    int dir_mult3 = 1;
                    if ((dir_mult1 == -1 || dir_mult2 == -1) && (dir_mult1 != dir_mult2)) {
                        dir_mult3 = -1;
                    }
                    bool z_forward = true;
                    if (max(Polyhedron::dot(camera.forward_vector(),{0,0,1}),Polyhedron::dot(camera.forward_vector(),{0,0,-1})) < max(Polyhedron::dot(camera.forward_vector(),{1,0,0}),Polyhedron::dot(camera.forward_vector(),{-1,0,0}))) {
                        z_forward = false;
                    }
                    //cout << dir_mult1 << " " << dir_mult2 << " " << dir_mult3 << " " << z_forward << endl;
                    //cout << "del: " << del << endl; 
                    int index;
                    switch( e.key.keysym.sym ) {
                        case SDLK_UP:
                        case SDLK_DOWN:
                            index = 2;
                            if (!z_forward) {
                                index = 0;
                            }
                            if (!block.space_down) {
                                if (block.select_dimension == 0) {
                                    block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                    while (!block.select_by_void()) {
                                        block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                    }
                                } else if (block.select_dimension == 1 || block.select_dimension == 2) {
                                    block.select[1] += direction*block.unit;
                                    while (!block.select_by_void()) {
                                        block.select[1] += direction*block.unit;
                                    }
                                }
                                for (int i = 0; i < 3; i++) {
                                    if (block.select[i] < 0) {
                                        block.select[i] = 0;
                                    }
                                    if (block.select[i] > block.size[i]-block.unit) {
                                        block.select[i] = block.size[i]-block.unit;
                                    }
                                }
                            } else {
                                if (block.select_dimension == 0) {
                                    block.select_size[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward && dir_mult1 != dir_mult2));
                                    if (block.select_size[index] == 0) {
                                        block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                        if (block.select_size[index] == block.unit) {
                                            block.select[index] -= block.unit;
                                            block.select_size[index] += block.unit;
                                        } else {
                                            block.select[index] += block.unit;
                                            block.select_size[index] -= block.unit;
                                        }
                                    }
                                } else if (block.select_dimension == 1 || block.select_dimension == 2) {
                                    block.select_size[1] += direction*block.unit;
                                    if (block.select_size[1] == 0) {
                                        block.select_size[1] = direction*block.unit;
                                        if (block.select_size[1] == block.unit) {
                                            block.select[1] -= block.unit;
                                            block.select_size[1] += block.unit;
                                        } else {
                                            block.select[1] += block.unit;
                                            block.select_size[1] -= block.unit;
                                        }
                                    }
                                }
                                for (int i = 0; i < 3; i++) {
                                    if (block.select_size[i] == 0) {
                                        block.select_size[i] = direction*dir_mult1*dir_mult3*block.unit*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                        if (block.select_size[i] == block.unit) {
                                            block.select[i] -= block.unit;
                                        } else {
                                            block.select[i] += block.unit;
                                            block.select_size[i] -= block.unit;
					}
                                    }
                                    if (block.select[i]+block.select_size[i] > block.size[i]) {
                                        block.select_size[i] = block.size[i]-block.select[i];
                                        if (block.select_size[i] == 0) {
                                            block.select_size[i] -= block.unit;
                                        }
                                    }
                                    if (block.select[i]+block.select_size[i] < 0) {
                                        block.select_size[i] = -block.select[i];
                                        if (block.select_size[i] == 0) {
                                            block.select_size[i] += block.unit;
                                        }
                                    }
                                    if (block.select[i] < block.origin[i]) {
                                        block.select[i] += block.unit;
                                    }
                                    if (block.select[i] > block.origin[i]+block.size[i]) {
                                        block.select[i] -= block.unit;
                                    }
                                }
                            }
                            break;

                        case SDLK_RIGHT:
                        case SDLK_LEFT:
                            index = 0;
                            if (!z_forward) {
                                index = 2;
                            }
                            if (!block.space_down) {
                                if (block.select_dimension == 0) {
                                    block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                    while (!block.select_by_void()) {
                                        block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                    }
                                } else if (block.select_dimension == 1) {
                                    block.select[0] += direction*dir_mult2*block.unit;
                                    while (!block.select_by_void()) {
                                        block.select[0] += direction*dir_mult2*block.unit;
                                    }
                                } else if (block.select_dimension == 2) {
                                    block.select[2] -= direction*dir_mult2*block.unit;
                                    while (!block.select_by_void()) {
                                        block.select[2] -= direction*dir_mult2*block.unit;
                                    }
                                }
                                for (int i = 0; i <= 2; i += 2) {
                                    if (block.select[i] > block.size[i]-block.unit) {
                                        block.select[i] = block.size[i]-block.unit;
                                    }
                                    if (block.select[i] < block.origin[i]) {
                                        block.select[i] = block.origin[i];
                                    }
                                }
                            } else {
                                if (block.select_dimension == 0) {
                                    block.select_size[index] += direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                    if (block.select_size[index] == 0) {
                                        block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                        if (block.select_size[index] == block.unit) {
                                            block.select[index] -= block.unit;
                                            block.select_size[index] += block.unit;
                                        } else {
                                            block.select[index] += block.unit;
                                            block.select_size[index] -= block.unit;
                                        }
                                    }
                                } else if (block.select_dimension == 1) {
                                    block.select_size[0] += direction*dir_mult1*block.unit;
                                    if (block.select_size[0] == 0) {
                                        block.select_size[0] = direction*dir_mult1*block.unit;
                                        if (block.select_size[0] == block.unit) {
                                            block.select[0] -= block.unit;
                                            block.select_size[0] += block.unit;
                                        } else {
                                            block.select[0] += block.unit;
                                            block.select_size[0] -= block.unit;
                                        }
                                    }
                                } else if (block.select_dimension == 2) {
                                    block.select_size[2] -= direction*dir_mult2*block.unit;
                                    if (block.select_size[2] == 0) {
                                        block.select_size[2] = -direction*dir_mult2*block.unit;
                                        if (block.select_size[2] == block.unit) {
                                            block.select[2] -= block.unit;
                                            block.select_size[2] += block.unit;
                                        } else {
                                            block.select[2] += block.unit;
                                            block.select_size[2] -= block.unit;
                                        }
                                    }
                                }
                                for (int i = 0; i <= 2; i += 2) {
                                    if (block.select[i]+block.select_size[i] > block.size[i]) {
                                        block.select_size[i] = block.size[i]-block.select[i];
                                    }
                                    if (block.select[i]+block.select_size[i] < 0) {
                                        block.select_size[i] = -block.select[i];
                                        if (block.select_size[i] == 0) {
                                            block.select_size[i] += block.unit;
                                        }
                                    }
                                    if (block.select[i] < block.origin[i]) {
                                        block.select[i] += block.unit;
                                    }
                                    if (block.select[i] > block.origin[i]+block.size[i]) {
                                        block.select[i] -= block.unit;
                                    }
                                }
                            }
                            break;
                    }
                } else if( e.type == SDL_KEYUP ) {
                    switch( e.key.keysym.sym ) {
                        case SDLK_SPACE: {
                            block.space_down = false;
                            Polyhedron::RoundPointMeter rounder(block.meters);
                            if (block.cylinder_tool) {
                                for (const Box& box : block.cylinder_boxes) {
                                    block.poly = box.del(block.poly);
                                }
                                block.poly.round_verts(rounder);
                                if (!del) {
                                    for (const Box& box : block.cylinder_boxes) {
                                        block.poly = box.add(block.poly);
                                    }
                                }
                            } else {
                                double i = block.select[0];
                                double i_max = block.select[0]+block.select_size[0];
                                if (block.select_size[0] < 0) {
                                    swap(i,i_max);
                                }
                                double j = block.select[1];
                                double j_max = block.select[1]+block.select_size[1];
                                if (block.select_size[1] < 0) {
                                    swap(j,j_max);
                                }
                                double k = block.select[2];
                                double k_max = block.select[2]+block.select_size[2];
                                if (block.select_size[2] < 0) {
                                    swap(k,k_max);
                                }
                                i = min(block.select[0],block.select[0]+block.select_size[0]);
                                j = min(block.select[1],block.select[1]+block.select_size[1]);
                                k = min(block.select[2],block.select[2]+block.select_size[2]);
                                std::array<double,3> size = {abs(block.select_size[0]),abs(block.select_size[1]),abs(block.select_size[2])};
                                Box box({i,i+size[0]},{j,j+size[1]},{k,k+size[2]});
                                block.poly = box.del(block.poly);
                                if (!del) {
                                    block.poly = box.add(block.poly);
                                }
                            }
                            for (const set<int>& face : block.poly.faces) {
                                for (const int& edge_index : face) {
                                    for (const int& index : block.poly.edges[edge_index]) {
                                        cout << "[" << block.poly.verts[index][0] << "," << block.poly.verts[index][1] << "," << block.poly.verts[index][2] << "] ";
                                    }
                                    cout << endl;
                                }
                            }
                            block.poly.round_verts(rounder);
                            block.flip();
                            block.select[0] = block.select[0]+block.select_size[0]-((block.select_size[0] > 0) ? block.unit : 0);
                            block.select[1] = block.select[1]+block.select_size[1]-((block.select_size[1] > 0) ? block.unit : 0);
                            block.select[2] = block.select[2]+block.select_size[2]-((block.select_size[2] > 0) ? block.unit : 0);
                            block.select_size[0] = block.unit;
                            block.select_size[1] = block.unit;
                            block.select_size[2] = block.unit;
                            break;
                        } case SDLK_LSHIFT: {
                            shift_down = false;   
                            break;
                        } case SDLK_LGUI: {
                        } case SDLK_RGUI: {
                            meta_down = false;
                            break;
                        }
                    }
                }
            }
            SDL_SetRenderDrawColor( gRenderer, 0x00, 0x00, 0x00, 0xFF );
            SDL_RenderClear( gRenderer );
            crosshair.draw(gRenderer);
            SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
            for (set<int> edge : cube.edges) {
                std::array<double,3> p1 = cube.verts[*(edge.begin())];
                std::array<double,3> p2 = cube.verts[*(++edge.begin())];
                std::array<double,2> p1_2D = camera.project(p1);
                std::array<double,2> p2_2D = camera.project(p2);
                //cout << p1_2D[0] << " " << p1_2D[1] << " " << p2_2D[0] << " " << p2_2D[1] << endl;
                SDL_RenderDrawLine( gRenderer, p1_2D[0]+SCREEN_WIDTH/2, p1_2D[1]*-1+SCREEN_HEIGHT/2, p2_2D[0]+SCREEN_WIDTH/2, p2_2D[1]*-1+SCREEN_HEIGHT/2 );
            }
            block.draw(gRenderer);
            //Update screen
            SDL_RenderPresent( gRenderer );
            auto stop = std::chrono::high_resolution_clock::now();
            //cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms." << endl;
        }
    }
}
