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
#include "ortools/constraint_solver/routing.h"
#include "ortools/constraint_solver/routing_parameters.h"
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
using namespace operations_research;

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

    static vector<std::array<double,3>> solve_tsp(const set<std::array<double, 3>>& point_set, std::array<Polyhedron,2>& polys) {
        vector<std::array<double,3>> points(point_set.begin(),point_set.end());
        int n = points.size();
        vector<vector<int>> dist_matrix = create_distance_matrix(points, polys);

        // Fix: Use NodeIndex type for depot
        const RoutingIndexManager::NodeIndex depot{0};
        RoutingIndexManager manager(n, 1, depot);  // 1 vehicle, depot = NodeIndex{0}
        RoutingModel routing(manager);

        const int transit_callback_index = routing.RegisterTransitCallback(
            [&dist_matrix, &manager](int64_t from_index, int64_t to_index) -> int64_t {
                // Fix: IndexToNode returns NodeIndex, need to convert to int
                RoutingIndexManager::NodeIndex from_node = manager.IndexToNode(from_index);
                RoutingIndexManager::NodeIndex to_node = manager.IndexToNode(to_index);
                return dist_matrix[from_node.value()][to_node.value()];
            });

        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);

        routing.AddDimension(transit_callback_index, 0, 1'000'000, true, "Distance");

        RoutingSearchParameters search_params = DefaultRoutingSearchParameters();
        search_params.set_first_solution_strategy(FirstSolutionStrategy::SEQUENTIAL_CHEAPEST_INSERTION);

        const Assignment* solution = routing.SolveWithParameters(search_params);

        if (!solution) return {};

        std::vector<int> ordered_indices;
        int64_t index = routing.Start(0);
        while (!routing.IsEnd(index)) {
            // Fix: Convert NodeIndex to int using .value()
            RoutingIndexManager::NodeIndex node = manager.IndexToNode(index);
            ordered_indices.push_back(node.value());
            index = solution->Value(routing.NextVar(index));
        }
        
        vector<std::array<double,3>> output;
        for (const int& index : ordered_indices) {
            output.push_back(points[index]);
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
        m(5,3) = 1;
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
    Polyhedron edge_sets_to_poly(std::array<vector<set<set<std::array<double,3>>>>,2> edge_sets_per_poly, std::array<Polyhedron,2> polys) {
        Polyhedron poly;
        set<std::array<double,3>> points;
        for (const vector<set<set<std::array<double,3>>>>& edge_sets : edge_sets_per_poly) {
            for (const set<set<std::array<double,3>>>& edges : edge_sets) {
                for (const set<std::array<double,3>>& edge : edges) {
                    for (const std::array<double,3>& point : edge) {
                        points.insert(point);
                    }
                }
            }
        }
        map<std::array<double,3>,std::array<double,3>> point_map;
        for (const std::array<double,3> point : points) {
            if (point_map.find(Polyhedron::round_point(point)) != point_map.end()) {
                if (Polyhedron::distance(point, Polyhedron::round_point(point)) < distance(point_map[round_point(point)], Polyhedron::round_point(point))) {
                    point_map[Polyhedron::round_point(point)] = point;
                }
            } else {
                point_map[Polyhedron::round_point(point)] = point;
            }
        }
        points.clear();
        for (const pair<std::array<double,3>,std::array<double,3>>& p : point_map) {
            points.insert(p.second);
            poly.verts.push_back(p.second);
        }
        for (int poly_index = 0; poly_index < 2; poly_index++) {
            for (int face_index = 0; face_index < polys[poly_index].faces.size(); face_index++) {
                set<std::array<double,3>> cofacial_points;
                for (const std::array<double,3>& point : points) {
                    bool any = false;
                    for (const std::array<std::array<double,3>,3>& triangle: Polyhedron::triangulate(*Polyhedron::circuit_cut(polys[poly_index].circuits(face_index)))) {
                        if (Polyhedron::inside_triangle(triangle, point)) {
                            any = true;
                            break;
                        }
                    }
                    if (any) {
                        cofacial_points.insert(point);
                    }
                }
                if (cofacial_points.size() > 2) {
                    vector<std::array<double,3>> path = solve_tsp(cofacial_points, polys);
                    if (!path.size()) {
                        continue;
                    }
                    set<int> new_face;
                    for (int i = 0; i < path.size(); i++) {
                        set<int> edge = {(int)std::distance(poly.verts.begin(),find(poly.verts.begin(),poly.verts.end(),path[(i-1+path.size())%path.size()])), (int)std::distance(poly.verts.begin(), find(poly.verts.begin(),poly.verts.end(),path[i]))};
                        vector<set<int>>::iterator it = find(poly.edges.begin(),poly.edges.end(),edge);
                        if (it != poly.edges.end()) {
                            new_face.insert(std::distance(poly.edges.begin(),it));
                        } else {
                            poly.edges.push_back(edge);
                            new_face.insert(poly.edges.size()-1);
                        }
                        poly.faces.push_back(new_face);
                    }
                }
            }
        }
        for (int vert_index = 0; vert_index < poly.verts.size();) {
            bool all = true;
            for (const set<int>& edge : poly.edges) {
                if (edge.find(vert_index) != edge.end()) {
                    all = false;
                    break;
                }
            }
            if (all) {
                poly.verts.erase(poly.verts.begin() + vert_index);
                for(int edge_index = 0; edge_index < poly.edges.size(); edge_index++) {
                    set<int>::iterator it = poly.edges[edge_index].begin();
                    int index1 = *it;
                    it++;
                    int index2 = *it;
                    if (index1 > vert_index) {
                        index1--;
                    }
                    if (index2 > vert_index) {
                        index2--;
                    }
                    poly.edges[edge_index] = {index1, index2};
                }
            } else {
                vert_index++;
            }
        }
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
            for (const std::array<double,3>& point : poly.verts) {
                if (!other_poly.is_inside(point)) {
                    all = false;
                }
            }
            if (all) {
                edge_sets_per_poly[poly_index].push_back(set<set<std::array<double,3>>>());
                for (const set<int>& edge : poly.edges) {
                    set<int>::iterator it = edge.begin();
                    int p_i1 = *it;
                    it++;
                    int p_i2 = *it;
                    edge_sets_per_poly[poly_index].front().insert({poly.verts[p_i1], poly.verts[p_i2]});
                }
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
                if (other_poly.in_faces(vert).size()) {
                    continue;
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
                        vector<int> in_faces = other_poly.in_faces(poly.verts[v_i2]);
                        if (in_faces.size()) {
                            leaves.insert(poly.verts[v_i2]);
                            for (const int& face_index : in_faces){
                                face_lookup[poly.verts[v_i2]].insert(face_index);
                            }
                            q.push(v_i2);
                        } else if (intersects.size()) {
                            sort(intersects.begin(), intersects.end(), Polyhedron::compare_face_intersections);
                            cout << "alpha " << intersects[0].alpha << endl;
                            leaves.insert(intersects[0].point);
                            for (int i = 0; i < intersects.size(); i++) {
                                if (Polyhedron::round_float(intersects[i].alpha) == Polyhedron::round_float(intersects[0].alpha)) {
                                    face_lookup[intersects[i].point].insert(intersects[i].face_index);
                                } else {
                                    break;
                                }
                            }
                            if (root_in_poly) {
                                new_edges.insert({poly.verts[v_i1], intersects[0].point});
                                cout << "point 1 " << poly.verts[1][0] << "," << poly.verts[1][1] << "," << poly.verts[1][2] << endl; 
                            } else {
                                for (int i = 1; i < intersects.size(); i++) {
                                    if (round_float(intersects[i].alpha) > round_float(intersects[0].alpha)) {
                                        cout << intersects[i].alpha << " " << intersects[0].alpha << endl;
                                        new_edges.insert({intersects[0].point, intersects[i].point});
                                        break;
                                    }
                                }
                            }
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
};

int main(int argc, char* argv[]) {
    string filename1 = "poly1.ply";
    if (argc > 1) {
        filename1 = argv[0];
    }
    Polyhedron poly1;
    poly1.load(filename1);
    string filename2 = "poly2.ply";
    if (argc > 2) {
        filename1 = argv[1];
    }
    Polyhedron poly2;
    poly2.load("poly2.ply");
    string filename3 = "out.txt";
    if (argc > 3) {
        filename3 = argv[2];
    }
    auto start = std::chrono::high_resolution_clock::now();
    Polyhedron poly3 = poly1.intersect(poly2);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time: " << duration.count() / 1000000.0 << " s" << std::endl;
    for (const std::array<double,3> vert : poly3.verts) {
        cout << "[" << vert[0] << "," << vert[1] << "," << vert[2] << "]" << endl;
    }
    return 0;
}
