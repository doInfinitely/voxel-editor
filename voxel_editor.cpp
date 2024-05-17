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
};

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
        for (int i = 0; i < edges.size(); i++) {
            for (const int& index : edges[i]) {
                std::array<double,3> point = verts[index];
                edge_lookup[point].insert(edges[i]);
            }
        }
        set<int> start = edges[*(face.begin())];
        set<int> previous = start;
        std::array<double,3> point = verts[*(start.begin())];
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
Polyhedron get_cube(std::array<double,3>displacement, std::array<double,3>factors) {
    return get_cube(displacement, factors, {0,0,0});
}

// Screen dimension constants
const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
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
        std::vector<int> xs;
        for (int i= 0;  i < points.size();  j = i++) {
            if ((points[i][1] < y && y <= points[j][1]) || (points[j][1] < y && y <= points[i][1])) {
                // the x value of the vertex + the height delta between the vertex and the scanline *the reciprocal of the slope =
                // the x value of the intersection of the preceding edge and the scanline 
                xs.push_back(round(points[i][0] + (y-points[i][1]) * (points[i][0]-points[j][0]) / (points[i][1]-points[j][1]) ));
            }
            // sort the intersections, the new item bubbles up to its place
            for (int k = xs.size()-1;  k && xs[k-1] > xs[k];  --k) {
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
    map<std::array<double,3>,set<std::array<double,3>>> contig;
    map<set<std::array<double,3>>,set<std::array<double,3>>> segments;
    map<set<std::array<double,3>>,Polyhedron> polys;
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
        for (const pair<set<std::array<double,3>>,Polyhedron>& p : polys) {
            for (int face_index = 0; face_index < p.second.faces.size(); face_index++) {
                vector<std::array<double,3>> points = p.second.circuit(face_index);
                vector<std::array<double,2>> points_2D;
                for (const std::array<double,3> point : points) {
                    points_2D.push_back(camera.project(point));
                    points_2D.back()[0] += SCREEN_WIDTH/2;
                    points_2D.back()[1] = points_2D.back()[1]*-1+SCREEN_HEIGHT/2; 
                }
                fill_polygon(gRenderer, points_2D);
            }
            for (const set<int>& edge : p.second.edges) {
                std::array<double,2> p1 = camera.project(p.second.verts[*(edge.begin())]);
                std::array<double,2> p2 = camera.project(p.second.verts[*(++edge.begin())]);
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
        camera.focal[2] = -1000;
        Polyhedron cube = get_cube({0,0,-1},{1,1,1},{0,0,0});
        Crosshair crosshair;
        bool mouse_down = false;
        int mousedown_x;
        int mousedown_y;
        SDL_ShowCursor(SDL_DISABLE);
        voxel_editor::Block block(3,3,3,1);
        // Main loop flag
        bool quit = false;

        // Event handler
        SDL_Event e;

        // While application is running
        while( !quit ) {
            auto start = std::chrono::high_resolution_clock::now();
            // Handle events on queue
            while( SDL_PollEvent( &e ) != 0 ) {
                //User requests quit
                if( e.type == SDL_QUIT ) {
                    quit = true;
                }
                if( e.type == SDL_MOUSEMOTION) {
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
                    
                }
                if (e.type == SDL_MOUSEBUTTONDOWN) {
                    SDL_GetMouseState( &mousedown_x, &mousedown_y );
                    mouse_down = true;
                }
                if (e.type == SDL_MOUSEBUTTONUP) {
                    mouse_down = false;
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
            cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms." << endl;
        }
    }
}
