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
    static bool coplanar(vector<std::array<double,3>> points) {
        if (points.size() < 4) {
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
    void construct_faces() {
        map<std::array<double,3>,set<int>> edge_lookup;
        for (int i = 0; i < edges.size(); i++) {
            for (const int& index : edges[i]) {
                std::array<double,3> point = verts[index];
                edge_lookup[point].insert(i);
            }
        }
        set<vector<int>> cycles;
        vector<vector<int>> queue;
        for (int i = 0; i < edges.size(); i++) {
            queue.push_back({i});
        }
        while (queue.size()) {
            vector<int> path = queue.back();
            queue.pop_back();
            for (const int& index : edges[path.back()]) {
                if (path.size() == 1 or edges[path[path.size()-2]].find(index) == edges[path[path.size()-2]].end()) {
                    std::array<double,3> point = verts[index];
                    for (const int& edge : edge_lookup[point]) {
                        if (find(path.begin(),path.end(),edge) == path.end()) {
                            vector<int> new_path;
                            for (int i = 0; i < path.size(); i++) {
                                new_path.push_back(path[i]);
                            }
                            new_path.push_back(edge);
                            if (new_path.size() > 2) {
                                set<std::array<double,3>> point_set;
                                for (const int& edge_index : new_path) {
                                    for (const int& index : edges[edge_index]) {
                                        point_set.insert(verts[index]);
                                    }
                                }
                                vector<std::array<double,3>> points(point_set.begin(),point_set.end());

                                if (Polyhedron::coplanar(points)) {
                                    set<int> edge_intersection;
                                    set_intersection(edges[new_path.front()].begin(), edges[new_path.front()].end(), edges[new_path.back()].begin(),edges[new_path.back()].end(),inserter(edge_intersection, edge_intersection.begin()));
                                    if (edge_intersection.size()) {
                                        cycles.insert(new_path);
                                    } else {
                                        queue.push_back(new_path);
                                    }
                                } 
                            } else {
                                queue.push_back(new_path);
                            }
                        }
                    }
                }
            }
        }
        faces.clear();
        for (const vector<int>& cycle : cycles) {
            set<int> temp(cycle.begin(),cycle.end());
            faces.push_back(temp);
        }
    }
    vector<std::array<double,3>> circuit(int face_index) const {
        set<int> face = faces[face_index];
        vector<std::array<double,3>> output;
        if (!face.size()) {
            return output;
        }
        map<std::array<double,3>,set<set<int>>> edge_lookup;
        for (const int& edge_index : face) {
            for (const int& index : edges[edge_index]) {
                //cout << "index: " << index << " " << edge_index << " " << edges.size() << endl;
                std::array<double,3> point = verts[index];
                edge_lookup[point].insert(edges[edge_index]);
                //cout << '[' << point[0] << "," << point[1] << "," << point[2] << "] " << edge_index << endl;
                
            }
        }
        set<int> start = edges[*(face.begin())];
        set<int> previous = start;
        std::array<double,3> point = verts[*(start.begin())];
        //cout << '[' << point[0] << "," << point[1] << "," << point[2] << "]" << endl;
        output.push_back(point);
        set<set<int>> set_diff;
        set<set<int>> temp = {start};
        set_difference(edge_lookup[point].begin(),edge_lookup[point].end(),temp.begin(),temp.end(),inserter(set_diff, set_diff.begin()));
        set<int> current = *(set_diff.begin());
        while (current != start) {
            set<int> set_diff2;
            set_difference(current.begin(),current.end(),previous.begin(),previous.end(),inserter(set_diff2, set_diff2.begin()));
            point = verts[*(set_diff2.begin())];
            output.push_back(point);
            previous = current;
            set_diff.clear();
            set<set<int>> temp = {current};
            set_difference(edge_lookup[point].begin(),edge_lookup[point].end(),temp.begin(),temp.end(),inserter(set_diff, set_diff.begin()));
            current = *(set_diff.begin());
        }
        return output;
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
        polyhedron.construct_faces();
        return polyhedron;
}

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
    std::array<double,3> size;
    Polyhedron block;
    map<std::array<double,3>,set<std::array<double,3>>*> contig;
    map<set<std::array<double,3>>,map<set<std::array<double,3>>,int>*> segments;
    map<set<std::array<double,3>>,Polyhedron*> polys;
    std::array<double,3> select = {0,0,0};
    int select_dimension = 0;
    double unit;
    std::array<double,3> select_size;
    Block(double width, double height, double depth, double unit) {
        size = {width, height, depth};
        block = get_cube({width/2,height/2,depth/2},{width,height,depth});
        this->unit = unit;
        select_size = {unit,unit,unit};
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
        for (const pair<set<std::array<double,3>>,Polyhedron*>& p : polys) {
            for (int face_index = 0; face_index < p.second->faces.size(); face_index++) {
                vector<std::array<double,3>> points = p.second->circuit(face_index);
                vector<std::array<double,2>> points_2D;
                for (const std::array<double,3> point : points) {
                    points_2D.push_back(camera.project(point));
                    points_2D.back()[0] += SCREEN_WIDTH/2;
                    points_2D.back()[1] = points_2D.back()[1]*-1+SCREEN_HEIGHT/2; 
                }
                fill_polygon(gRenderer, points_2D);
            }
        }
        SDL_SetRenderDrawColor( gRenderer, 0x00, 0x00, 0x00, 0xFF );
        for (const pair<set<std::array<double,3>>,Polyhedron*>& p : polys) {
            for (const set<int>& edge : p.second->edges) {
                std::array<double,2> p1 = camera.project(p.second->verts[*(edge.begin())]);
                std::array<double,2> p2 = camera.project(p.second->verts[*(++edge.begin())]);
                p1[0] += SCREEN_WIDTH/2;
                p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
                p2[0] += SCREEN_WIDTH/2;
                p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
                SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
            }
        }
        double width = size[0];
        double height = size[1];
        double depth = size[2];
        Polyhedron axes;
        if (select_dimension == 0) {
            axes = get_cube({width/2,unit/2+select[1],depth/2}, {width,unit,depth});
        } else if (select_dimension == 1) {
            axes = get_cube({width/2,height/2,unit/2+select[2]}, {width,height,unit});
        } else if (select_dimension == 2) {
            axes = get_cube({unit/2+select[0],height/2,depth/2}, {unit,height,depth});
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
        Polyhedron select_cube = get_cube({select[0]+select_size[0]/2,select[1]+select_size[1]/2,select[2]+select_size[2]/2}, select_size);
        SDL_SetRenderDrawColor( gRenderer, 0x80, 0x80, 0x80, 0xFF );
        for (const set<int>& edge : select_cube.edges) {
            std::array<double,2> p1 = camera.project(select_cube.verts[*(edge.begin())]);
            std::array<double,2> p2 = camera.project(select_cube.verts[*(++edge.begin())]);
            p1[0] += SCREEN_WIDTH/2;
            p1[1] = p1[1]*-1+SCREEN_HEIGHT/2;
            p2[0] += SCREEN_WIDTH/2;
            p2[1] = p2[1]*-1+SCREEN_HEIGHT/2;
            SDL_RenderDrawLine(gRenderer, p1[0], p1[1], p2[0], p2[1]);
        }
        SDL_SetRenderDrawColor( gRenderer, 0x80, 0x80, 0x80, 0xFF );
        if (select[2] > 0 ) {
            Polyhedron shadow = get_cube({select[0]+select_size[0]/2,select[1]+select_size[1]/2,0+select[2]/2}, {select_size[0],select_size[1],select[2]});
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
        if (select[1] > 0 ) {
            Polyhedron shadow = get_cube({select[0]+select_size[0]/2,0+select[1]/2,select[2]+select_size[2]/2}, {select_size[0],select[1],select_size[2]});
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
        if (select[0] > 0 ) {
            Polyhedron shadow = get_cube({0+select[0]/2,select[1]+select_size[1]/2,select[2]+select_size[2]/2}, {select[0],select_size[1],select_size[2]});
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
    }
    void construct_poly(set<std::array<double,3>> cont) {
        map<set<std::array<double,3>>,int>seg = *segments[cont];
        map<std::array<double,3>,set<set<std::array<double,3>>>> segment_map;
        for (const pair<set<std::array<double,3>>,int>& p : seg) {
            if ((p.second+4) % 2 == 1) {
                for (const std::array<double,3> point : p.first) {
                    segment_map[point].insert(p.first);
                }
            }
        }
        for (const pair<std::array<double,3>,set<set<std::array<double,3>>>>& p : segment_map) {
            bool found = true;
            while (found) {
                found = false;
                set<std::array<double,3>> seg1;
                set<std::array<double,3>> seg2;
                set<std::array<double,3>> seg3;
                
                for (const set<std::array<double,3>>&segment1 : segment_map[p.first]) {
                    for (const set<std::array<double,3>>&segment2 : segment_map[p.first]) {
                        set<std::array<double,3>> point_set;
                        point_set.insert(segment1.begin(),segment1.end());
                        point_set.insert(segment2.begin(),segment2.end());
                        vector<std::array<double,3>> points(point_set.begin(),point_set.end());
                        if (segment1 != segment2 && Polyhedron::colinear(points)) {
                            vector<std::array<double,3>> temp;
                            set_symmetric_difference(segment1.begin(),segment1.end(),segment2.begin(),segment2.end(),inserter(temp, temp.begin()));
                            set<std::array<double,3>> segment3(temp.begin(),temp.end());
                            found = true;
                            seg1.insert(segment1.begin(),segment1.end());
                            seg2.insert(segment2.begin(),segment2.end());
                            seg3.insert(segment3.begin(),segment3.end());
                            break;
                        
                        }
                    }
                    if (found) {
                        break;
                    }
                }
                if (found) {
                    for (const std::array<double,3>& point : seg1 ) {
                        segment_map[point].erase(seg1);
                    }
                    for (const std::array<double,3>& point : seg2 ) {
                        segment_map[point].erase(seg2);
                    }
                    for (const std::array<double,3>& point : seg3 ) {
                        segment_map[point].insert(seg3);
                    }
                }
            }
        }
        set<set<std::array<double,3>>> edges;
        set<std::array<double,3>> to_remove;
        for (const pair<std::array<double,3>,set<set<std::array<double,3>>>>& p : segment_map) {
            if (p.second.size()) {
                edges.insert(p.second.begin(),p.second.end());
            } else {
                to_remove.insert(p.first);
            }
        }
        for (const std::array<double,3> x : to_remove) {
            segment_map.erase(x);
        }
        set<std::array<double,3>> verts;
        for (const set<std::array<double,3>>& edge : edges) {
            verts.insert(edge.begin(),edge.end());
        }
        map<set<std::array<double,3>>,set<set<set<std::array<double,3>>>>>face_map;
        for (const pair<std::array<double,3>,set<set<std::array<double,3>>>>& p : segment_map) {
            for (const set<std::array<double,3>>& edge1 : p.second) {
                for (const set<std::array<double,3>>& edge2 : p.second) {
                    set<std::array<double,3>> edge_intersection;
                    set_intersection(edge1.begin(),edge1.end(),edge2.begin(),edge2.end(),inserter(edge_intersection, edge_intersection.begin()));
                    if (edge1 != edge2 && edge_intersection.size()) {
                        face_map[edge1].insert({edge1,edge2});
                        face_map[edge2].insert({edge1,edge2});
                    }
                }
            }
        }
        for (const pair<set<std::array<double,3>>,set<set<set<std::array<double,3>>>>>& p : face_map) {
            bool found = true;
            while (found) {
                found = false;
                set<set<std::array<double,3>>> f1;
                set<set<std::array<double,3>>> f2;
                set<set<std::array<double,3>>> f3;
                for (const set<set<std::array<double,3>>>& face1 : face_map[p.first]) {
                    for (const set<set<std::array<double,3>>>& face2 : face_map[p.first]) {
                        set<set<std::array<double,3>>> face3;
                        face3.insert(face1.begin(),face1.end());
                        face3.insert(face2.begin(),face2.end());
                        vector<std::array<double,3>> points;
                        for (const set<std::array<double,3>>& edge : face3) {
                            for (const std::array<double,3> point : edge) {
                                points.push_back(point);
                            }
                        }
                        if (face1 != face2 && Polyhedron::coplanar(points)) {
                            found = true;
                            f1.insert(face1.begin(),face1.end());
                            f2.insert(face2.begin(),face2.end());
                            f3.insert(face3.begin(),face3.end());
                            break;
                        }
                    }
                    if (found) {
                        break;
                    }
                }
                if (found) {
                    for (const set<std::array<double,3>>& edge : f1) {
                        face_map[edge].erase(f1);
                    }
                    for (const set<std::array<double,3>>& edge : f2) {
                        face_map[edge].erase(f2);
                    }
                    for (const set<std::array<double,3>>& edge : f3) {
                        face_map[edge].insert(f3);
                    }
                }
            }
        }
        set<set<set<std::array<double,3>>>> faces;
        for (const pair<set<std::array<double,3>>,set<set<set<std::array<double,3>>>>> p : face_map) {
            faces.insert(p.second.begin(),p.second.end());
        }
        vector<std::array<double,3>> verts_vector(verts.begin(),verts.end());
        vector<set<std::array<double,3>>> edges_vector(edges.begin(),edges.end());
        Polyhedron* poly = new Polyhedron();
        
        for (const set<set<std::array<double,3>>>& face : faces) {
            set<int> temp;
            for (const set<std::array<double,3>>& edge : face) {
                vector<set<std::array<double,3>>>::iterator it = find(edges_vector.begin(),edges_vector.end(), edge);
                temp.insert(it - edges_vector.begin());
            }
            poly->faces.push_back(temp);
        }
        for (const set<std::array<double,3>>& edge : edges) {
            set<int> temp;
            for (const std::array<double,3>& vert : edge) {
                vector<std::array<double,3>>::iterator it = find(verts_vector.begin(),verts_vector.end(),vert);
                temp.insert(it - verts_vector.begin());
            }
            poly->edges.push_back(temp);
        }
        poly->verts = verts_vector;
        polys[cont] = poly; 
    }
    void add_poly(std::array<double,3>pos, std::array<double,3>size, bool reconstruct) {
        set<set<std::array<double,3>>> edited;
        for(double i = pos[0]; i < pos[0]+size[0]; i+= unit) {
            for(double j = pos[1]; j < pos[1]+size[1]; j+= unit) {
                for(double k = pos[2]; k < pos[2]+size[2]; k+= unit) {
                    contig[{i,j,k}] = new set<std::array<double,3>>;
                    contig[{i,j,k}]->insert({i,j,k});
                    const vector<std::array<double,3>> deltas = {{unit,0,0},{0,unit,0},{0,0,unit}};
                    const std::array<int,2> multipliers = {-1,1};
                    for (const std::array<double,3>& delta : deltas) {
                        for (const int& m : multipliers) {
                            std::array<double,3> other = {i+m*delta[0],j+m*delta[1],k+m*delta[2]};
                            if (contig.find(other) != contig.end() && contig[{i,j,k}] != contig[other]) {
                                //cout << i << " " << j << " " << k << " other: " << other[0] << " " << other[1] << " " << other[2] << endl;
                                set<std::array<double,3>> temp(contig[other]->begin(),contig[other]->end());
                                edited.insert(temp);
                                contig[other]->insert(contig[{i,j,k}]->begin(),contig[{i,j,k}]->end());
                                contig[{i,j,k}] = contig[other];
                                for (const std::array<double,3>& other : *contig[{i,j,k}]) {
                                    contig[other] = contig[{i,j,k}];
                                } 
                            }
                        }
                    }
                }
            }
        }
        set<std::array<double,3>>cont = *contig[pos];
        map<set<std::array<double,3>>,int>* seg = new map<set<std::array<double,3>>,int>;
        for (const set<std::array<double,3>>& e : edited) {
            if (polys.find(e) != polys.end()) {
                polys.erase(e);
            }
            if (segments.find(e) != segments.end()) {
                for(const pair<set<std::array<double,3>>,int>& p : *segments[e]) {
                    (*seg)[p.first] += p.second;
                }
            }
        }
        for (int i1 = 0; i1 < 3; i1++) {
            std::array<double,3> p = {pos[0],pos[1],pos[2]};
            std::array<double,3> start = {p[0],p[1],p[2]};
            while (p[i1] < start[i1]+size[i1]) { 
                std::array<double,3> end = {p[0], p[1], p[2]};
                end[i1] += unit;
                set<std::array<double,3>>segment = {p,end};
                (*seg)[segment] += 1;
                p[i1] += unit;
            }
            for (int i2 = 0; i2 < 3; i2++) {
                if (i1 < i2) {
                    std::array<double,3> p = {pos[0],pos[1],pos[2]};
                    p[i1] += size[i1];
                    std::array<double,3> start = {p[0],p[1],p[2]};
                    while (p[i2] < start[i2]+size[i2]) {
                        std::array<double,3> end = {p[0], p[1], p[2]};
                        end[i2] += unit;
                        set<std::array<double,3>>segment = {p,end};
                        (*seg)[segment] += 1;
                        p[i2] += unit;
                    }
                    p[0] = pos[0], p[1] = pos[1], p[2] = pos[2];
                    p[i2] += size[i2];
                    start[0] = p[0], start[1] = p[1], start[2] = p[2];
                    while (p[i1] < start[i1]+size[i1]) {
                        std::array<double,3> end = {p[0], p[1], p[2]};
                        end[i1] += unit;
                        set<std::array<double,3>>segment = {p,end};
                        (*seg)[segment] += 1;
                        p[i1] += unit;
                    }
                    for (int i3 = 0; i3 < 3; i3++) {
                        if (i1 != i3 && i2 != i3) {
                            std::array<double,3> p = {pos[0],pos[1],pos[2]};
                            p[i1] += size[i1];
                            p[i2] += size[i2];
                            start[0] = p[0], start[1] = p[1], start[2] = p[2];
                            while (p[i3] < start[i3]+size[i3]) {
                                std::array<double,3> end = {p[0], p[1], p[2]};
                                end[i3] += unit;
                                set<std::array<double,3>>segment = {p,end};
                                (*seg)[segment] += 1;
                                p[i3] += unit;
                            }
                        }
                    }
                }
            }
        }
        segments[cont] = seg;
        cout << segments.size() << endl;
        if (reconstruct) {
            this->construct_poly(cont);
        }
    }
    void del_poly(std::array<double,3> pos, std::array<double,3> size) {
        cout << "hello" << endl;
        set<std::array<double,3>>* cont = contig[pos];
        polys.erase(*cont);
        set<std::array<double,3>> removed;
        for (double i=pos[0]; i < pos[0]+size[0]; i += unit) {
            for (double j=pos[1]; j < pos[1]+size[1]; j += unit) {
                for (double k=pos[2]; k < pos[2]+size[2]; k += unit) {
                    contig.erase({i,j,k});
                    cont->erase({i,j,k});
                    removed.insert({i,j,k});
                }    
            }    
        }
        set<std::array<double,3>> visited;
        vector<std::array<double,3>> queue;
        bool triple_break = false;
        for (const std::array<double,3>& item : removed) {
            const vector<std::array<double,3>> deltas = {{unit,0,0},{0,unit,0},{0,0,unit}};
            const std::array<int,2> multipliers = {-1,1};
            for (const std::array<double,3>& delta : deltas) {
                for (const int& m : multipliers) {
                    std::array<double,3> other = {item[0]+m*delta[0],item[1]+m*delta[1],item[2]+m*delta[2]};
                    if (cont->find(other) != cont->end()) {
                        queue.push_back(other);
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
        while (queue.size()) {
            std::array<double,3> item = queue.back();
            queue.pop_back();
            visited.insert(item);
            const vector<std::array<double,3>> deltas = {{unit,0,0},{0,unit,0},{0,0,unit}};
            const std::array<int,2> multipliers = {-1,1};
            for (const std::array<double,3>& delta : deltas) {
                for (const int& m : multipliers) {
                    std::array<double,3> other = {item[0]+m*delta[0],item[1]+m*delta[1],item[2]+m*delta[2]};
                    if (cont->find(other) != cont->end() and visited.find(other) == visited.end()) {
                        queue.push_back(other);
                    }
                }
            }
        }
        set<std::array<double,3>> cont2(cont->begin(),cont->end());
        cont2.insert(removed.begin(),removed.end()); 
        map<set<std::array<double,3>>,int>* segs = segments[cont2];
        segments.erase(cont2);
        if (visited.size() == cont->size()) {
            const vector<std::array<double,3>> deltas = {{unit,0,0},{0,unit,0},{0,0,unit}};
            for (const std::array<double,3>& item : removed) {
                for (const std::array<double,3>& d1 : deltas) {
                    set<std::array<double,3>> seg = {item, {item[0]+d1[0],item[1]+d1[1],item[2]+d1[2]}};
                    (*segs)[seg] -= 1;
                    for (const std::array<double,3>& d2 : deltas) {
                        if (d1[0] < d2[0] || (d1[0] == d2[0] && d1[1] < d2[1]) ) {
                            cout << d1[0] << d1[1] << d1[2] << endl;
                            cout << d2[0] << d2[1] << d2[2] << endl;
                            set<std::array<double,3>> seg1 = {{item[0]+d1[0],item[1]+d1[1],item[2]+d1[2]},{item[0]+d1[0]+d2[0],item[1]+d1[1]+d2[1],item[2]+d1[2]+d2[2]}};
                            (*segs)[seg1] -= 1;
                            set<std::array<double,3>> seg2 = {{item[0]+d2[0],item[1]+d2[1],item[2]+d2[2]},{item[0]+d1[0]+d2[0],item[1]+d1[1]+d2[1],item[2]+d1[2]+d2[2]}};
                            (*segs)[seg2] -= 1;
                            for (const std::array<double,3>& d3 : deltas) {
                                if (d1 != d3 && d2 != d3) {
                                    //cout << d1[0] << d1[1] << d1[2] << endl;
                                    //cout << d2[0] << d2[1] << d2[2] << endl;
                                    //cout << d3[0] << d3[1] << d3[2] << endl;
                                    set<std::array<double,3>> seg = {{item[0]+d1[0]+d2[0],item[1]+d1[1]+d2[1],item[2]+d1[2]+d2[2]},{item[0]+d1[0]+d2[0]+d3[0],item[1]+d1[1]+d2[1]+d3[1],item[2]+d1[2]+d2[2]+d3[2]}};
                                    (*segs)[seg] -= 1;
                                }
                            }
                        }
                    }
                }
            }
            segments[*cont] = segs;
            this->construct_poly(*cont);
        } else {
            for (const std::array<double,3>& other : *cont) {
                contig.erase(other);
            }
            set<std::array<double,3>> edited;
            for (const std::array<double,3>& other : *cont) {
                this->add_poly(other,{unit,unit,unit},false);
                edited.insert(other);
            }
            set<set<std::array<double,3>>> temp;
            for (const std::array<double,3>& other : edited) {
                temp.insert(*contig[other]);
            }
            for (const set<std::array<double,3>>& c : temp) {
                this->construct_poly(c);
            }
        }
        cout << "goodbye" << endl; 
    }
    voxel_editor::Block subdivide(int divisor) {
        voxel_editor::Block b(size[0], size[1], size[2], unit/divisor);
        b.contig = contig;
        b.select = select;
        b.polys = polys;
        for (const pair<set<std::array<double,3>>,map<set<std::array<double,3>>,int>*>& p1 : segments) {
            map<set<std::array<double,3>>,int>* segs = new map<set<std::array<double,3>>,int>;
            for (const pair<set<std::array<double,3>>,int>& p2 : *p1.second) {
                std::array<double,3> point1 = *(p2.first.begin());
                std::array<double,3> point2 = *(++p2.first.begin());
                std::array<double,3> diff = {(point2[0]-point1[0])/divisor,(point2[1]-point1[1])/divisor,(point2[2]-point1[2])/divisor};
                while (point1 != point2) {
                    cout << point1[0] << " " << point1[1] << " " << point1[2] << " ";
                    cout << point2[0] << " " << point2[1] << " " << point2[2] << endl;
                    set<std::array<double,3>> seg = {point1, {point1[0]+diff[0], point1[1]+diff[1], point1[2]+diff[2]}};
                    (*segs)[seg] = p2.second;
                    point1[0] += diff[0];
                    point1[1] += diff[1];
                    point1[2] += diff[2];
                }
            }
            b.segments[p1.first] = segs;
        }
        set<set<std::array<double,3>>> contiguous;
        for (pair<std::array<double,3>,set<std::array<double,3>>*> p : contig) {
            cout << p.first[0] << " " << p.first[1] << " " << p.first[2] << endl;
            contiguous.insert(*(p.second));
        }
        for (const set<std::array<double,3>>& cont : contiguous) {
            for (const std::array<double,3>& point : cont) {
                for (double i = point[0]; i < point[0]+unit; i += b.unit) {
                    for (double j = point[1]; j < point[1]+unit; j += b.unit) {
                        for (double k = point[2]; k < point[2]+unit; k += b.unit) {
                            cout << i << " " << j << " " << k << endl;
                            if (b.contig.find({i,j,k}) == b.contig.end()) {
                                b.contig[point]->insert({i,j,k});
                                b.contig[{i,j,k}] = b.contig[point];
                            }
                        }
                    }
                }
            }
            if (b.segments.find(cont) != b.segments.end()) {
                set<std::array<double,3>>* new_cont = b.contig[*(cont.begin())];
                b.segments[*new_cont] = b.segments[cont];
                b.polys[*new_cont] = b.polys[cont];
                b.segments.erase(cont);
                b.polys.erase(cont);
            }
        }
        return b;
    }
    bool select_by_void() {
        if (contig.find(select) == contig.end()) {
            return true;
        }
        const vector<std::array<double,3>> deltas = {{unit,0,0},{0,unit,0},{0,0,unit}};
        const std::array<int,2> multipliers = {-1,1};
        for (const std::array<double,3>& delta : deltas) {
            for (const int& m : multipliers) {
                std::array<double,3> other = {select[0]+m*delta[0],select[1]+m*delta[1],select[2]+m*delta[2]};
                if (contig.find(other) == contig.end()) {
                    return true;
                }
            }
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
        bool space_down = false;
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
                } else if( e.type == SDL_KEYDOWN ) {
                    vector<std::array<double,3>> camera_items;
                    int direction = 1;
                    switch( e.key.keysym.sym ) {
                        case SDLK_a:
                        camera_items = Polyhedron::rotate({camera.origin, camera.focal, camera.x_vector, camera.y_vector}, {0,M_PI/180*10,0});
                        camera.origin = camera_items[0];
                        camera.focal = camera_items[1];
                        camera.x_vector = camera_items[2];
                        camera.y_vector = camera_items[3];
                        break;

                        case SDLK_d:
                        camera_items = Polyhedron::rotate({camera.origin, camera.focal, camera.x_vector, camera.y_vector}, {0,-M_PI/180*10,0});
                        camera.origin = camera_items[0];
                        camera.focal = camera_items[1];
                        camera.x_vector = camera_items[2];
                        camera.y_vector = camera_items[3];
                        break;

                        case SDLK_w:
                        camera_items = Polyhedron::rotate({camera.origin, camera.focal, camera.x_vector, camera.y_vector}, {M_PI/180*10,0,0});
                        camera.origin = camera_items[0];
                        camera.focal = camera_items[1];
                        camera.x_vector = camera_items[2];
                        camera.y_vector = camera_items[3];
                        break;

                        case SDLK_s:
                        camera_items = Polyhedron::rotate({camera.origin, camera.focal, camera.x_vector, camera.y_vector}, {-M_PI/180*10,0,0});
                        camera.origin = camera_items[0];
                        camera.focal = camera_items[1];
                        camera.x_vector = camera_items[2];
                        camera.y_vector = camera_items[3];
                        break;

                        case SDLK_q:
                        camera_items = Polyhedron::rotate({camera.origin, camera.focal, camera.x_vector, camera.y_vector}, {0,0,-M_PI/180*10});
                        camera.origin = camera_items[0];
                        camera.focal = camera_items[1];
                        camera.x_vector = camera_items[2];
                        camera.y_vector = camera_items[3];
                        break;

                        case SDLK_e:
                        camera_items = Polyhedron::rotate({camera.origin, camera.focal, camera.x_vector, camera.y_vector}, {0,0,M_PI/180*10});
                        camera.origin = camera_items[0];
                        camera.focal = camera_items[1];
                        camera.x_vector = camera_items[2];
                        camera.y_vector = camera_items[3];
                        break;
                        
                        case SDLK_LSHIFT:
                        block.select_dimension = (block.select_dimension+3-1)%3;
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
                            space_down = true;
                        break;
                        case SDLK_DOWN:
                        case SDLK_LEFT:
                        direction = -1;
                        break;
                        
                        case SDLK_2:
                        block = block.subdivide(2);
                        break;

                        case SDLK_3:
                        block = block.subdivide(3);
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
                    cout << dir_mult1 << " " << dir_mult2 << " " << dir_mult3 << " " << z_forward << endl;
                    int index;
                    switch( e.key.keysym.sym ) {
                        case SDLK_UP:
                        case SDLK_DOWN:
                        index = 2;
                        if (!z_forward) {
                            index = 0;
                        }
                        if (!space_down) {
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
                                block.select_size[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2));
                            } else if (block.select_dimension == 1 || block.select_dimension == 2) {
                                block.select_size[1] += direction*block.unit;
                            }
                            if (block.select_size[index] == 0) {
                                block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                            }
                            if (block.select_size[1] == 0) {
                                block.select_size[1] = direction*block.unit;
                            }
                            for (int i = 0; i < 3; i++) {
                                if (block.select[i]+block.select_size[i] > block.size[i]) {
                                    block.select_size[i] = block.size[i]-block.select[i];
                                }
                                if (block.select[i]+block.select_size[i] < 0) {
                                    block.select_size[i] = -block.select[i];
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
                        if (!space_down) {
                            if (block.select_dimension == 0) {
                                block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                while (!block.select_by_void()) {
                                    block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                }
                            } else if (block.select_dimension == 1) {
                                block.select[0] += direction*dir_mult1*block.unit;
                                while (!block.select_by_void()) {
                                    block.select[0] += direction*dir_mult1*block.unit;
                                }
                            } else if (block.select_dimension == 2) {
                                block.select[2] += direction*dir_mult2*block.unit;
                                while (!block.select_by_void()) {
                                    block.select[2] += direction*dir_mult2*block.unit;
                                }
                            }
                            for (int i = 0; i <= 2; i += 2) {
                                if (block.select[i] > block.size[i]-block.unit) {
                                    block.select[i] = block.size[i]-block.unit;
                                }
                                if (block.select[i] < 0) {
                                    block.select[i] = 0;
                                }
                            }
                        } else {
                            if (block.select_dimension == 0) {
                                block.select_size[index] += direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                if (block.select_size[index] == 0) {
                                    block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(2*(int)(z_forward)-1)*(1-2*(int)(z_forward && dir_mult1 != dir_mult2));
                                }
                            } else if (block.select_dimension == 1) {
                                block.select_size[0] += direction*dir_mult1*block.unit;
                                if (block.select_size[0] == 0) {
                                    block.select_size[0] = direction*dir_mult1*block.unit;
                                }
                            } else if (block.select_dimension == 2) {
                                block.select_size[2] -= direction*dir_mult2*block.unit;
                                if (block.select_size[2] == 0) {
                                    block.select_size[2] = -direction*dir_mult2*block.unit;
                                }
                            }
                            for (int i = 0; i <= 2; i += 2) {
                                if (block.select[i]+block.select_size[i] > block.size[i]) {
                                    block.select_size[i] = block.size[i]-block.select[i];
                                }
                                if (block.select[i]+block.select_size[i] < 0) {
                                    block.select_size[i] = -block.select[i];
                                }
                            }
                        }
                        break;
                    }
                } else if( e.type == SDL_KEYUP ) {
                    switch( e.key.keysym.sym ) {
                        case SDLK_SPACE:
                        space_down = false;
                        bool to_add = false;
                        double i = block.select[0];
                        double i_max = block.select[0]+block.select_size[0];
                        if (block.select_size[0] < 0) {
                            swap(i,i_max);
                        }
                        while (i < i_max) {
                            double j = block.select[1];
                            double j_max = block.select[1]+block.select_size[1];
                            if (block.select_size[1] < 0) {
                                swap(j,j_max);
                            }
                            while (j < j_max) {
                                double k = block.select[2];
                                double k_max = block.select[2]+block.select_size[2];
                                if (block.select_size[2] < 0) {
                                    swap(k,k_max);
                                }
                                while (k < k_max) {
                                    if (block.contig.find({i,j,k}) == block.contig.end()) {
                                        to_add = true;
                                    }
                                    k += block.unit;
                                }
                                j += block.unit;
                            }
                            i += block.unit;
                        }
                        i = min(block.select[0],block.select[0]+block.select_size[0]);
                        double j = min(block.select[1],block.select[1]+block.select_size[1]);
                        double k = min(block.select[2],block.select[2]+block.select_size[2]);
                        std::array<double,3> size = {abs(block.select_size[0]),abs(block.select_size[1]),abs(block.select_size[2])};
                        if (to_add) {
                            block.add_poly({i,j,k},size, true);
                        } else {
                            block.del_poly({i,j,k},size);
                        }
                        block.select[0] = block.select[0]+block.select_size[0]-block.unit;
                        block.select[1] = block.select[1]+block.select_size[1]-block.unit;
                        block.select[2] = block.select[2]+block.select_size[2]-block.unit;
                        block.select_size[0] = block.unit;
                        block.select_size[1] = block.unit;
                        block.select_size[2] = block.unit;
                        break;
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
