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

using namespace Eigen;
using namespace std;

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
            std::array<double,3> temp1 = {p[0],p[1]*cos(angles[0])-p[2]*sin(angles[0]),p[2]*sin(angles[0])+p[2]*cos(angles[0])};
            std::array<double,3> temp2 = {temp1[0]*cos(angles[1])-temp1[2]*sin(angles[1]),temp1[1],temp1[0]*sin(angles[1])+temp1[2]*cos(angles[1])};
            std::array<double,3> temp3 = {temp2[0]*cos(angles[2])-temp2[1]*sin(angles[2]),temp2[0]*sin(angles[2])+temp2[1]*cos(angles[2]),temp2[2]};
            output.push_back(temp3);
        }
        return output;
    }
    static double dot(std::array<double,3> vector1, std::array<double,3> vector2) {
        double output = 0;
        for (int i = 0; i < vector1.size(); i++) {
            output += vector1[i]*vector2[i];
        }
        return output;
    }
    static std::array<double,3> cross3D(std::array<double,3> vector1, std::array<double,3> vector2) {
        std::array<double,3> output;
        for (int i = 0; i < 3; i++) {
            output[i] = vector1[(i+1)%3]*vector2[(i+2)%3]-vector1[(i+2)%3]*vector2[(i+1)%3];
        }
        return output
    }
    set<vector<std::array<double,3>>> circuits(int face_index, int start, int previous, int current, vector<std::array<double,3>> path, set<vector<std::array<double,3>>> old_circuits) {
        set<vector<std::array<double,3>>> output;
        if (current == start) {
            output.insert(path);
            return output;
        }
        set<std::array<double,3>> path_set(path.begin(),path.end());
        for (const vector<std::array<double,3>>& old_circuit : old_circuits) {
            set<std::array<double,3>> old_circuit_set(old_circuit.begin(),old_circuit.end());
            set<std::array<double,3>> difference = (path_set.begin(), path_set.end(), old_circuit_set.begin(), old_circuit_set.end(), std::inserter(difference, difference.begin()));
            if (!difference.size()) {
                return output; 
            }
        }
        set<int> face = faces[face_index];
        map<int,set<set<int>>> edge_lookup;
        for (const int& edge_index : face) {
            edge = edges[edge_index];
            for (const int& index : edge) {
                edge_lookup[index].insert(edge);
            }
        }
        set<int> difference = set_difference(edges[current].begin(), edges[current].end(), edges[previous].begin(), edges[previous].end(), std::inserter(difference, difference.begin()));
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
            if (Box::coplanar(union_set)) {
                vector<set<int>>::iterator it = find(edges.begin(),edges.end(),y);
                current = std::distance(edges.begin(), it);
                old_circuits.insert(output.begin(),output.end());
                intermediate = this.circuits(face_index, start, previous, current, path, old_circuits);
                output.insert(intermediate.begin(),intermediate.end());
            }
        }
        vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
        for (int i = 0; i < output_list.size(); i++) {
            for (int j = i+1; j < output_list.size(); j++) {
                set<std::array<double,3>> set_i(output_list[i].begin(),output_list[i].end());
                set<std::array<double,3>> set_j(output_list[j].begin(),output_list[j].end());
                set<std::array<double,3>> difference = set_difference(set_i.begin(), set_i.end(), set_j.begin(), set_j.end(), std::inserter(difference, difference.begin()));
                if (output_list[i] != output_list[j] && !difference.size()) {
                    vector<std::array<double,3>> y_r(output_list[j].begin(),output_list[j].end());
                    std::reverse(y_r.begin(),y_r.end());
                    vector<std::array<double,3>>::iterator it = find(output_list[j].begin(),output_list[j].end(), output_list[i][0])
                    int index_x0_in_y = std::distance(output_list[j].begin(), it);
                    it = find(y_r.begin(),y_r.end(), output_list[i][0]);
                    int index_x0_in_y_r = std::distance(y_r.begin(), it);
                    vector<std::array<double,3>> y_rearranged(output_list[j].begin()+index_x0_in_y,output_list[j].end());
                    y_rearranged.insert(output_list[j].begin(),output_list[j].begin()+index_x0_in_y);
                    vector<std::array<double,3>> y_rearranged_sliced(y_rearranged.begin(),y_rearranged.begin()+output_list[i].size());
                    vector<std::array<double,3>> y_r_rearranged(y_r.begin()+index_x0_in_y_r,y_r.end());
                    y_r_rearranged.insert(y_r.begin(),y_r.begin()+index_x0_in_y_r);
                    vector<std::array<double,3>> y_r_rearranged_sliced(y_r_rearranged.begin(),y_r_rearranged.begin()+output_list[i].size());
                    if (y_rearranged_sliced == output_list[i] || y_r_rearranged == output_list[i]) {
                        if (output.find(output_list[j]) != output.end()) {
                            output.erase(output_list[j])
                        }
                    }
                }
            }
        }
        return output;
    }
    set<vector<std::array<double,3>>> circuits(int face_index) {
        set<int> face = faces[face_index];
        map<int,set<set<int>>> edge_lookup;
        for (const int& edge_index : face) {
            edge = edges[edge_index];
            for (const int& index : edge) {
                edge_lookup[index].insert(edge);
            }
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
                intermediate = this.circuits(face_index, start, start, current, path, output);
                output.insert(intermediate.begin(),intermediate.end());
                for (const vector<std::array<double,3>>& circuit : intermediate) {
                    for (int i = 0; i < circuit.size(); i++) {
                        seen.insert({circuit[(i-1)%circuit.size()],circuit[i]});
                    }
                }
            }
        }
        vector<vector<std::array<double,3>>> output_list(output.begin(),output.end());
        for (int i = 0; i < output_list.size(); i++) {
            for (int j = i+1; j < output_list.size(); j++) {
                set<std::array<double,3>> set_i(output_list[i].begin(),output_list[i].end());
                set<std::array<double,3>> set_j(output_list[j].begin(),output_list[j].end());
                set<std::array<double,3>> difference = set_difference(set_i.begin(), set_i.end(), set_j.begin(), set_j.end(), std::inserter(difference, difference.begin()));
                if (output_list[i] != output_list[j] && !difference.size()) {
                    vector<std::array<double,3>> y_r(output_list[j].begin(),output_list[j].end());
                    std::reverse(y_r.begin(),y_r.end());
                    vector<std::array<double,3>>::iterator it = find(output_list[j].begin(),output_list[j].end(), output_list[i][0])
                    int index_x0_in_y = std::distance(output_list[j].begin(), it);
                    it = find(y_r.begin(),y_r.end(), output_list[i][0]);
                    int index_x0_in_y_r = std::distance(y_r.begin(), it);
                    vector<std::array<double,3>> y_rearranged(output_list[j].begin()+index_x0_in_y,output_list[j].end());
                    y_rearranged.insert(output_list[j].begin(),output_list[j].begin()+index_x0_in_y);
                    vector<std::array<double,3>> y_rearranged_sliced(y_rearranged.begin(),y_rearranged.begin()+output_list[i].size());
                    vector<std::array<double,3>> y_r_rearranged(y_r.begin()+index_x0_in_y_r,y_r.end());
                    y_r_rearranged.insert(y_r.begin(),y_r.begin()+index_x0_in_y_r);
                    vector<std::array<double,3>> y_r_rearranged_sliced(y_r_rearranged.begin(),y_r_rearranged.begin()+output_list[i].size());
                    if (y_rearranged_sliced == output_list[i] || y_r_rearranged == output_list[i]) {
                        if (output.find(output_list[j]) != output.end()) {
                            output.erase(output_list[j])
                        }
                    }
                }
            }
        }
        return output;
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
            b(i) = point[i] 
        }
        for (int i = 0; i < 3; i++) {
            m(3,i) = 1;
        }
        b(3) = 1;
        VectorXd x = m.colPivHouseholderQr().solve(b);
        VectorXd b_prime = m*x;
        for (int i=0; i < 3; i++) {
            if ( abs(b_prime(i)-b(i)) > 0.1*abs(b(i)) ) {
               return false;
            }
        }
        if (x(0) >= 0 && x(1) >= 0 && x(2) >= 0) {
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
            cross_product_dot_normal.push_back(Polyhedron::dot(Polyhedron::cross_product_triplet(circuit[(i-1)%circuit.size()],circuit[i],circuit[(i+1)%circuit.size()]),normal));
        }
        vector<double> signs;
        for (const double& x : cross_product_dot_normal) {
            signs.push_back(copysign(1,x));
        }
        vector<bool> output;
        for (int i = 0; i < signs.size(); i++) {
            output.push_back(signs[i]==signs[(i-1)%signs.size()] || signs[i]==signs[(i+1)%signs.size()]);
        }
        return output;
        
    }
    struct EarClip {
        std::array<std::array<double,3>,3> ear;
        vector<std::array<double,3>> remainder;
    }
    static EarClip clip_ear(vector<std::array<double,3>> circuit) {
        vector<bool> is_convex1 = Polyhedron::convex_angles(circuit);
        vector<bool> is_convex2;
        for (const bool& x : is_convex1) {
            is_convex2.push_back(!x);
        }
        std::array<vector<bool>,2> is_convexes = {is_convex1,is_convex2};
        for (const vector<bool>& is_convex : is_convexes) {
            for (int i = 0; i < circuit.size(); i++) {
                if (is_convex[i]) {
                    std::array<std::array<double,3>,3> triangle = {circuit[(i-1)%circuit.size()],circuit[i],circuit[(i+1)%circuit.size()]};
                    bool no_break = true;
                    for (int j = 0; j < circuit.size(); j++) {
                        if (j != i && j != (i+1)%circuit.size() && j != (i-1)%circuit.size()) {

                            if (Polyhedron.inside_triangle(triangle,circuit[j]) && circuit[j] != circuit[(i-1)%circuit.size()] && circuit[j] != circuit[i] && circuit[j] != circuit[(i+1)%circuit.size()]) {
                                no_break = false;
                                break
                            }
                        }
                    }
                    if (no_break) {
                        EarClip ear_clip;
                        ear_clip.ear = triangle;
                        for (int k = 0; k < circuit.size(); k++) {
                            if (k != i) {
                                ear_clip.remainder.push_back(circuit[k]);
                            }
                        }
                        return ear_clip;
                    }
                }
            }
        }
    }
    static vector<std::array<std::array<double,3>,3>> triangulate(vector<std::array<double,3>> circuit) {
        vector<std::array<std::array<double,3>,3>> output;
        vector<std::array<double,3>> remainder = circuit;
        while (remainder.size() >3) {
            EarClip ear_clip = Polyhedron::clip_ear(remainder);
            remainder = ear_clip.remainder;
            output.push_back(ear_clip.ear);
        }
        output.push_back({remainder[0],remainder[1],remainder[2]});
        return output;
    }
    static vector<std::array<double,3>>* find_exterior_circuit(set<vector<std::array<double,3>>> circuits) {
        vector<vector<std::array<double,3>>> circuits_list(circuits.begin(),circuits.end());
        for (int i = 0; i < circuits_list.size(); i++) {
            vector<std::array<std::array<double,3>,3>> triangulation = Polyhedron.triangulate(circuits_list[i])
            for (int j = 0; j < circuits_list.size(); j++) {
                if (j != i) {
                    bool double_break = false;
                    for (const std::array<double,3>& point : circuits_list[j]) {
                        bool any = false;
                        for (const std::array<std::array<double,3>,3>& z : triangulation) {
                            if (Polyhedron.inside_triangle(z,point)) {
                                any = true;
                                break;
                            }
                        }
                        if (!any) {
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
                vector<std::array<double,3>>* output = new vector<std::array<double,3>;
                output->insert(circuits_list[i].begin(),circuits_list[i].end());
                return output
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
                std::array<double,3> p3 = circuit[(i-1)%circuit.size()];
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
                for (int j=0; j < 3; j++) {
                    if ( abs(b_prime(j)-b(j)) > 0.1*abs(b(j)) ) {
                        continue;
                    }
                }
                Polyhedron::CircuitIntersection circuit_intersection;
                circuit_intersection.alpha = x(0);
                for (int j=0; j < 3; j++) {
                    circuit_intersection.edge[j] = x(0)*p2[j]+(1-x(0))*p1[j];
                }
                circuit_intersection.circuit = circuit;
                output.push_back(circuit_intersection)
            }
        }
        return output;
    }
    static bool intersections_comp(Polyhedron::CircuitIntersection a, Polyhedron::CircuitIntersection b) {
        return a.alpha < b.alpha;
    }
    static vector<std::<array,3>>* circuit_cut(set<vector<std::array<double,3>>> circuits) {
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
            for (const std::array<double,3>& x : interior) {
                bool double_break = false;
                for (const std::array<double,3>& y : output[output.size()-1]) {
                    std::array<std::array<double,3>,2> segment = {x,y};
                    set<vector<std::array<double,3>>> interior_circuits(output.begin(),output.end()-1);
                    vector<Polyhedron::CircuitIntersection> intersections = Polyhedron::circuit_intersect(segment, interior_circuits);
                    if (!intersections.size()) {
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
                if (intersection.circuit == output[output.size()-1]) {
                    first_exterior_intersection = intersection;
                    break;
                }
            }
            Polyhedron::CircuitIntersection last_interior_intersection;
            for (const Polyhedron::CircuitIntersection& intersection : intersections) {
                if (intersection.circuit == output[output.size()-1]) {
                    continue;
                }
                last_interior_intersection = intersection;
            }
            int i;
            for (i = 0; i < output[output.size()-1].size(); i++) {
                if (find(first_exterior_intersection.edge.begin(), first_exterior_intersection.edge.end(), output[output.size()-1][i]) != first_exterior_intersection.edge.end() && find(first_exterior_intersection.edge.begin(), first_exterior_intersection.edge.end(), output[output.size()-1][(i-1)%output[output.size()-1].size()]) != first_exterior_intersection.edge.end()) {
                    break;
                }
            }
            int j;
            for (j = 0; j < last_interior_intersection.circuit.size(); j++) {
                if (find(last_interior_intersection.edge.begin(), last_interior_intersection.edge.end(), last_interior_intersection.circuit[j]) != last_interior_intersection.edge.end() && find(last_interior_intersection.edge.begin(), last_interior_intersection.edge.end(), last_interior_intersection.circuit[(j-1)%last_interior_intersection.circuit.size()]) != last_interior_intersection.edge.end()) {
                    break;
                }
            }
            vector<std::array<double,3>> temp(output[output.size()-1].begin(),output[output.size()-1].begin()+i);
            temp.push_back(first_exterior_intersection.point);
            temp.push_back(last_interior_intersection.point);
            for(int k = j-1; k >= 0; k--) {
                temp.push_back(last_interior_intersection.circuit[k]);
            }
            for(int k = last_interior_intersection.circuit.size()-1; k >= j; k--) {
                temp.push_back(last_interior_intersection.circuit[k]);
            }
            temp.push_back(last_interior_intersection.point);
            temp.push_back(first_exterior_intersection.point);
            temp.insert(output[output.size()-1].begin()+i,output.end());
            output[output.size()-1] = temp;
            vector<vector<std::array<double,3>>>::iterator it = find(output.begin(),output.end(),last_interior_intersection.circuit);
            int index = std::distance(output.begin(),it);
            output.erase(output.begin()+index);
            for(i = 1; i < output[output.size()-1].size()+1;){
                if (output[output.size()-1][i%output[output.size()-1].size()] == output[output.size()-1][i-1]) {
                    output[output.size()-1].erase(output[output.size()-1].begin()+i%output[output.size()-1].size());
                }
            }
        }
        vector<std::array<double,3>>* output_pointer = new vector<std::array<double,3>>;
        output_pointer->insert(output[0].begin(),output[0].end());
        return output_pointer;
    }
    struct IsInsideIntersection {
        double gamma;
        std::array<double,3> point;
        int face_index;
    }
    bool is_inside(std::array<double,3> point) {
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(*Polyhedron::circuit_cut(this.circuits(face_index)))) {
                if (Polyhedron::inside_triangle(triangle, point)) {
                    return true;
                }
            }
        }
        std::array<double,3> vec = {(double)(rand()/RAND_MAX),(double)(rand()/RAND_MAX),(double)(rand()/RAND_MAX)};
        vector<Polyhedron::IsInsideIntersection> output;
        for (int face_index = 0; face_index < faces.size(); face_index++) {
            vector<std::array<double,3>> circuit = *Polyhedron::circuit_cut(this.circuits(face_index));
            for (const std::array<std::array<double,3>,3>& triangle : Polyhedron::triangulate(circuit)) {
                MatrixXd m(4,4);
                VectorXd b(4);
                for (int i = 0; i < 3; i++) {
                    m(i,0) = triangle[0][i];
                    m(i,1) = triangle[1][i];
                    m(i,2) = triangle[2][i];
                    m(i,3) = -vec[i]
                    b(i) = point[i] 
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
                    continue
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
    static bool colinear(vector<std::array<double,3>> points) {
        if (points.size() < 3) {
            return true;
        }
        MatrixXd m(3,1);
        for (int i=0; i < 3; i++) {
            m(i,0) = points[1][i]-points[0][i];
        }
        int c = 0;
        for (const std::array<double,3>& p : points) {
            if (c > 1) {
                VectorXd b(3);
                for (int i=0; i < 3; i++) {
                    b(i) = p[i]-points[0][i];
                }
                VectorXd x = m.colPivHouseholderQr().solve(b);
                VectorXd b_prime = m*x;
                for (int i=0; i < 3; i++) {
                    if ( abs(b_prime(i)-b(i)) > 0.1*abs(b(i)) ) {
                        return false;
                    }
                }
            }
            c++;
        }
        return true;
    }
    static bool colinear(set<std::array<double,3>> points) {
        vector<std::array<double,3>> points_vector(points.begin(),points.end());
        return Box::colinear(points_vector);
    }
    static bool coplanar(vector<std::array<double,3>> points) {
        if (points.size() < 4) {
            return true;
        }
        if (Box::colinear(points)) {
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
            for (int j=0; j < 3; j++) {
                if ( abs(b_prime(j)-b(j)) > 0.1*abs(b(j)) ) {
                    return false;
                }
            }
        }
        return true;
    }
    static bool coplanar(set<std::array<double,3>> points) {
        vector<std::array<double,3>> points_vector(points.begin(),points.end());
        return Box::coplanar(points_vector);
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
    static std::array<double,3>* intersect_segmeents(set<std::array<double,3>> edge1, set<std::array<double,3>> edge2) {
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
            m(i,1) = p4[i]-p3[i]
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
               *(output)[i] = x(0)*p1[i]+(1-x(0))*p2[i]
            }
            return output;
        }
    }
    bool intersect(set<std::array<double,3>> edge) {
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
        std::array<std::array<double,3>,2> segment_shadow = {{p1[0],0,p1[2]},{p2[0],0,p2[2]}};
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
            set<vector<std::array<double,3>>> circuits = poly.circuits(face_index);
            vector<std::array<double,3>> circuit = Polyhedron.circuit_cut(circuits);
            set<vector<std::array<double,3>>> interior_circuits = (circuits);
            vector<std::array<double,3>> exterior_circuit = Polyhedron::find_exterior_circuit(circuits);
            interior_circuits.erase(exterior_circuit);
            vector<std::array<std::array<double,3>,3>> triangles = Polyhedron::triangulate(circuit);
            for (const <std::array<std::array<double,3>,3>& triangle : triangles) {
                if (Polyhedron::inside_triangle(triangle, point)) {
                    return true;
                }
            }
            if (circuit[0][2] == circuit[1][2] && circuit[0][2] == circuit[2][2] && circuit[0][2] > point[2]) {
                for (const <std::array<std::array<double,3>,3>& triangle : triangles) {
                    if (Polyhedron::inside_triangle(triangle, point)) {
                        return true;
                    }
                }
                std::array<double,3> projection = {point[0],point[1],circuit[0][2]};
                bool projection_inside_circuit = false;
                for (const <std::array<std::array<double,3>,3>& triangle : triangles) {
                    if (Polyhedron::inside_triangle(triangle, projection)) {
                        project_inside_circuit = true;
                    }
                }
                bool projection_inside_interior_circuit = false;
                for (const vector<std::array<double,3>>& interior_circuit : interior_circuits) {
                    triangles = Polyhedron::triangulate(interior_circuit);
                    for (const <std::array<std::array<double,3>,3>& triangle : triangles) {
                        if (Polyhedron::inside_triangle(triangle, projection)) {
                            projection_inside_interior_circuit = true;
                        }
                    }
                }
                if (projection_inside_circuit && !projection_inside_interior_circuit) {
                    intersections++;
                }
            }
        }
        return intersections % 2 == 1;
    }
    Polyhedron del(Polyhedron poly) {
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
        if (!Box::inside_polyhedron(poly, center)) {
            bool triple_break = false;
            for (const set<int>& edge : poly.edges) {
                set<std::array<double,3>> edge_grounded;
                for (const int& index : edge) {
                    edge_grounded.insert(poly.verts[index]);
                }
                if (this.intersect(edge_grounded)) {
                    bool any1 = false;
                    for (const double& x1 : x_min_max) {
                        for (const double& y1 : y_min_max) {
                            for (const double& z1 : z_min_max) {
                                for (const double& x2 : x_min_max) {
                                    for (const double& y2 : y_min_max) {
                                        for (const double& z2 : z_min_max) {
                                            if ((int)(x1 != x2)+(int)(y1 != y2)+(int)(z1 != z2) == 1) {
                                                bool all = true;
                                                for (const std::array<double,3>& point : edge_grounded) {
                                                    set<std::array<double,3>> temp_set = {{x1,y1,z1},{x2,y2,z2},point};
                                                    if (!Box::colinear(temp_set)) {
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
                            if (Box::coplanar(temp_set)) {
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
            DelPreprocessOutput output = this.del_preprocess(poly, i);
            faces.push_back(output.new_edges);
            path_edges.insert(output.path_edges.begin(),output.path_edges.end());
            if (output.box_index > -1) {
                box_map[output.box_index] = i;
            }
        }
        vector<int> old_face_indices;
        for (int i = 0; i < faces.size(); i++) {
            if (faces[i].size()) {
                old_face_indices.push_back(i)
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
                    box = {{x_min, y_min, z_min},{x_min, y_max, z_min},{x_min, y_max, z_max},{x_min, y_min, z_max}};
                    break;
                case 1:
                    box = {{x_max, y_min, z_min},{x_max, y_max, z_min},{x_max, y_max, z_max},{x_max, y_min, z_max}};
                    break;
                case 2:
                    box = {{x_min, y_min, z_min},{x_max, y_min, z_min},{x_max, y_min, z_max},{x_min, y_min, z_max}};
                    break;
                case 3:
                    box = {{x_min, y_max, z_min},{x_max, y_max, z_min},{x_max, y_max, z_max},{x_min, y_max, z_max}};
                    break;
                case 4:
                    box = {{x_min, y_min, z_min},{x_max, y_min, z_min},{x_max, y_max, z_min},{x_min, y_max, z_min}}
                    break;
                case 5:
                    box = {{x_min, y_min, z_max},{x_max, y_min, z_max},{x_max, y_max, z_max},{x_min, y_max, z_max}}
                    break;
            }
            if (box_map[i] == -1) {
                set<set<std::array<double,3>>> new_face;
                for (int j = 0; j < box.size(); j++) {
                    new_face.insert
                }
                faces.push_back()
            }
            else if (box_map[i] != -2) {
                
            }
        }
    }
};
