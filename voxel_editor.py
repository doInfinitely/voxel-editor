from sympy import nsolve, solve, symbols, Eq, diff
import sympy
from math import sin, cos, pi, sqrt, copysign
import random
import numpy as np
import multiprocessing as mp
import queue
from PIL import Image
from scipy.optimize import fsolve
from scipy.spatial import ConvexHull
import warnings
warnings.simplefilter("ignore")
import matplotlib.pyplot as plt
import operator

def shift(points, displacement=(0,0,0)):
    output = []
    for p in points:
        output.append(tuple(p[i]+displacement[i] for i in range(3)))
    return output

def scale(points, factors=(1,1,1)):
    output = []
    if type(factors) == tuple or type(factors) == list:
        for p in points:
            output.append(tuple(p[i]*factors[i] for i in range(3)))
    else:
        for p in points:
            output.append(tuple(p[i]*factors for i in range(3)))
    return output

def rotate(points, angles=(0,0,0)):
    output = []
    for p in points:
        p = (p[0],p[1]*cos(angles[0])-p[2]*sin(angles[0]),p[1]*sin(angles[0])+p[2]*cos(angles[0]))
        p = (p[0]*cos(angles[1])-p[2]*sin(angles[1]),p[1],p[0]*sin(angles[1])+p[2]*cos(angles[1]))
        p = (p[0]*cos(angles[2])-p[1]*sin(angles[2]),p[0]*sin(angles[2])+p[1]*cos(angles[2]),p[2])
        output.append(p)
    return output

def dot(vector1, vector2):
    return sum(x*vector2[i] for i,x in enumerate(vector1))

class Camera:
    def __init__(self):
        self.origin = (0,0,0)
        self.x_vector = (1,0,0)
        self.y_vector = (0,1,0)
        self.focal = (0,0,-1)
        self.zoom = 1
    def project(self, point):
        o = self.origin
        x = self.x_vector
        y = self.y_vector
        f = self.focal
        p = point
        a = np.array([[x[i], y[i], p[i]-f[i]] for i in range(3)], dtype=np.float64)
        b = np.array([p[i]-o[i] for i in range(3)], dtype=np.float64)
        try:
            x = np.linalg.solve(a, b)
        except np.linalg.LinAlgError:
            return None
        x = x.tolist()
        return x[0]*self.zoom, x[1]*self.zoom
    def move(self, displacement):
        self.origin = tuple(self.origin[i]+displacement[i] for i in range(3))
        self.focal = tuple(self.focal[i]+displacement[i] for i in range(3))
    def rotate(self, angles):
        p = tuple(self.focal[i]-self.origin[i] for i in range(3))
        p, self.x_vector, self.y_vector = rotate([p, self.x_vector, self.y_vector], angles)
        self.focal = tuple(p[i]+self.origin[i] for i in range(3))
    def look(self, point):
        vector1 = (self.origin[0]-point[0], self.origin[1]-point[1], self.origin[2]-point[2])
        mag1 = sqrt(dot(vector1,vector1))
        unit = (vector1[0]/mag1,vector1[1]/mag1,vector1[2]/mag1)
        vector2 = (self.focal[0]-self.origin[0], self.focal[1]-self.origin[1], self.focal[2]-self.origin[2])
        mag2 = sqrt(dot(vector2,vector2))
        focal = (unit[0]*mag2+self.origin[0],unit[1]*mag2+self.origin[1],unit[2]*mag2+self.origin[2])
        theta = symbols("theta1 theta2 theta3")
        point = symbols("point1 point2 point3")
        focal[1], self.focal[1]
        f = vector2
        p1 = (f[0],f[1]*sympy.cos(theta[0])-f[2]*sympy.sin(theta[0]),f[1]*sympy.sin(theta[0])+f[2]*sympy.cos(theta[0]))
        p2 = (p1[0]*sympy.cos(theta[1])-p1[2]*sympy.sin(theta[1]),p1[1],p1[0]*sympy.sin(theta[1])+p1[2]*sympy.cos(theta[1]))
        p3 = (p2[0]*sympy.cos(theta[2])-p2[1]*sympy.sin(theta[2]),p2[0]*sympy.sin(theta[2])+p2[1]*sympy.cos(theta[2]),p2[2])
        eqs = [Eq(focal[i],p3[i]+self.origin[i]) for i in range(3)]
        try:
            solutions = nsolve(eqs, theta, (0.1,0.1,0.1), dict=True)
            solution = solutions[0]
        except (ValueError, ZeroDivisionError):
            print("error")
            return
        self.x_vector, self.y_vector = rotate([self.x_vector, self.y_vector], tuple(solution[theta[i]] for i in range(3)))
        self.focal = focal
    def forward_vector(self):
        vec = tuple(self.origin[i]-self.focal[i] for i in range(3))
        mag = sqrt(dot(vec,vec))
        return tuple(vec[i]/mag for i in range(3))

class Polyhedron:
    def __init__(self):
        self.verts = []
        self.edges = []
        self.faces = []
    # return sequential points on the face
    def colinear(points):
        if len(points) < 3:
            return True
        for p in points[2:]:
            a = np.array([[points[1][i]-points[0][i]] for i in range(3)])
            b = np.array([p[i]-points[0][i] for i in range(3)])
            x, res, rank, s = np.linalg.lstsq(a, b)
            if not np.allclose(np.dot(a, x), b, rtol=0.1):
                return False
        return True
    def coplanar(points):
        if len(points) < 4:
            return True
        points = list(points)
        '''
        alpha, beta = symbols("alpha beta")
        exprs = [alpha*(points[1][i]-points[0][i])+beta*(points[2][i]-points[0][i]) for i in range(3)]
        for p in points[3:]:
            eqs = [Eq(expr,p[i]-points[0][i]) for i,expr in enumerate(exprs)]
            if not len(solve(eqs)):
                return False
        '''
        triple_break = False
        for i,x in enumerate(points):
            for j,y in enumerate(points):
                for k,z in enumerate(points):
                    if i != j and i != k and j != k and not Polyhedron.colinear([x,y,z]):
                        ind = [i,j,k]
                        triple_break = True
                        break
                if triple_break:
                    break
            if triple_break:
                break
        for index,p in enumerate(points):
            if index in ind:
                continue
            a = np.array([[points[ind[1]][i]-points[ind[0]][i],points[ind[2]][i]-points[ind[0]][i]] for i in range(3)])
            b = np.array([p[i]-points[ind[0]][i] for i in range(3)])
            x, res, rank, s = np.linalg.lstsq(a, b)
            #print(a, b, x, res)
            if not np.allclose(np.dot(a, x), b, rtol=0.1):
                return False
        return True
    def construct_faces(self):
        edge_lookup = dict()
        for edge_index, edge in enumerate(self.edges):
            for index in edge:
                point = self.verts[index]
                if point not in edge_lookup:
                    edge_lookup[point] = set()
                edge_lookup[point].add(edge_index)
        cycles = set()
        queue = [(index,) for index,edge in enumerate(self.edges)]
        while len(queue):
            #print([len(x) for x in queue], [len(x) for x in cycles])
            path = queue.pop()
            for index in self.edges[path[-1]]:
                if len(path) == 1 or index not in self.edges[path[-2]]:
                    point = self.verts[index]
                    for edge in edge_lookup[point]:
                        if edge not in path:
                            new_path = path + (edge,)
                            points = {self.verts[index] for edge_index in new_path for index in self.edges[edge_index]}
                            if len(new_path) > 2:
                                if Polyhedron.coplanar(points):
                                    if len(self.edges[new_path[0]] & self.edges[new_path[-1]]):
                                        cycles.add(frozenset(new_path))
                                    else:
                                        queue.append(new_path)
                            else:
                                queue.append(new_path)
        self.faces = list(cycles)
    def circuit(self, face_index):
        face = self.faces[face_index]
        if not len(face):
            return []
        circuit = []
        edge_lookup = dict()
        for edge_index in face:
            edge = self.edges[edge_index]
            for index in self.edges[edge_index]:
                point = self.verts[index]
                #print(point, edge)
                if point not in edge_lookup:
                    edge_lookup[point] = set()
                edge_lookup[point].add(edge)
        start = self.edges[list(face)[0]]
        previous = start
        point = self.verts[list(start)[0]]
        circuit.append(point)
        #print([[self.verts[index] for index in self.edges[edge_index]] for edge_index in face])
        current = list(edge_lookup[point] - set([start]))[0]
        #print([[self.verts[index] for index in edge] for edge in self.edges])
        while current != start:
            point = self.verts[list(current - previous)[0]]
            circuit.append(point)
            previous = current
            #print(point, [[self.verts[index] for index in edge] for edge in edge_lookup[point]])
            #print('circuit', circuit)
            current = list(edge_lookup[point] - set([current]))[0]
        return circuit
    # returns any edge that intersects the line segment between p1 and p2
    def intersect(self, p1, p2):
        output = []
        for edge in self.edges:
            p3, p4 = tuple(self.verts[index] for index in edge)
            alpha, beta = symbols("alpha beta")
            eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], beta*p3[i]+(1-beta)*p4[i]) for i in range(3)]
            solutions = solve(eqs, dict=True)
            if len(solutions):
                alpha, beta = solutions[0][alpha], solutions[0][beta]
                if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1:
                    output.append(edge)
        return output
    def inside_triangle(triangle, point):
        triangle = list(triangle)
        alpha, beta = symbols("alpha beta")
        #print(triangle)
        exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
        exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
        eqs = [Eq(x,point[i]) for i,x in enumerate(exprs)]
        solutions = solve(eqs)
        if len(solutions):
            alpha, beta = solutions[0][alpha], solutions[0][beta]
            if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1:
                return True
        return False
    def triangulation(circuit):
        output = []
        for i,x in enumerate(circuit):
            for j,y in enumerate(circuit):
                for k,z in enumerate(circuit):
                    points = set([x,y,z])
                    if len(points) != 3:
                        continue
                    '''
                    for w in circuit:
                        if w not in points and Polyhedron.inside_triangle(points, w):
                            continue
                    '''
                    output.append(points)
                    polygons = (circuit[k:]+circuit[:i+1], circuit[i:j+1], circuit[j:k+1])
                    for polygon in polygons:
                        if len(polygon) == 3:
                            output.append(set(polygon))
                        else:
                                output.extend(Polyhedron.triangulation(polygon))
                    return output
        return output
    def pinpoint(triangle, position, direction):
        alpha_range = [0,1]
        beta_range = [0,1]
        gamma = 0
        min_delta = 0.001
        prev_distance = float("inf")
        while True:
            distance = 0
            alpha_center = (alpha_range[0]+alpha_range[1])/2
            beta_center = (beta_range[0]+beta_range[1])/2
            delta = min_delta
            while True:
                exprs = [alpha_center*triangle[0][i]+(1-alpha_center)*triangle[1][i] for i in range(3)]
                exprs = [beta_center*x+(1-beta_center)*triangle[2][i] for i,x in enumerate(exprs)]
                exprs = [expr-(position[i]+direction[i]*(gamma+delta)) for i,expr in enumerate(exprs)]
                distance = sum(exprs)
                #print(prev_distance, distance, delta)
                if distance > prev_distance or distance < 0:
                    break
                prev_distance = distance
                delta *= 2
            delta /= 2
            if delta >= min_delta:
                gamma += delta
            #print(alpha_center, beta_center, gamma, delta, prev_distance)
            beta = beta_center
            alpha = (alpha_range[0]+alpha_center)/2
            exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
            exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
            exprs = [(expr-(position[i]+direction[i]*gamma))**2 for i,expr in enumerate(exprs)]
            distance1 = sqrt(sum(exprs))
            alpha = (alpha_range[1]+alpha_center)/2
            exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
            exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
            exprs = [(expr-(position[i]+direction[i]*gamma))**2 for i,expr in enumerate(exprs)]
            distance2 = sqrt(sum(exprs))
            if distance1 < distance2:
                alpha_range[1] = alpha_center
            else:
                alpha_range[0] = alpha_center
            alpha = alpha_center
            beta = (beta_range[0]+beta_center)/2
            exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
            exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
            exprs = [(expr-(position[i]+direction[i]*gamma))**2 for i,expr in enumerate(exprs)]
            distance1 = sqrt(sum(exprs))
            beta = (beta_range[1]+beta_center)/2
            exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
            exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
            exprs = [(expr-(position[i]+direction[i]*gamma))**2 for i,expr in enumerate(exprs)]
            distance2 = sqrt(sum(exprs))
            if distance1 < distance2:
                beta_range[1] = beta_center
            else:
                beta_range[0] = beta_center
            if delta < min_delta and alpha_range[1]-alpha_range[0] < min_delta and beta_range[1]-beta_range[0] < min_delta:
                break
        alpha, beta = (alpha_range[0]+alpha_range[1])/2, (beta_range[0]+beta_range[1])/2
        exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
        exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
        exprs = [(expr-(position[i]+direction[i]*gamma))**2 for i,expr in enumerate(exprs)]
        distance = sqrt(sum(exprs))
        return alpha, beta, gamma, distance

    # Cast a ray from position in direction, return where on the polyhedron it hits
    def ray(self, position, direction):
        output = []
        for face_index,face in enumerate(self.faces):
            circuit = self.circuit(face_index)
            min_distance, mini = float("inf"), None
            for triangle in Polyhedron.triangulation(circuit):
                triangle = list(triangle)
                alpha, beta, gamma = symbols("alpha beta gamma")
                exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
                exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
                eqs = [Eq(expr,position[i]+direction[i]*gamma) for i,expr in enumerate(exprs)]
                '''
                a,b,g = 0,0,0
                exprs = [(x-position[i]+direction[i]*gamma)**2 for i,x in enumerate(exprs)]
                for i in range(1000):
                    da = diff(sum(exprs), alpha).subs(alpha,a).subs(beta,b).subs(gamma,g)
                    db = diff(sum(exprs), beta).subs(alpha,a).subs(beta,b).subs(gamma,g)
                    dg = diff(sum(exprs), gamma).subs(alpha,a).subs(beta,b).subs(gamma,g)
                    a -= da*0.0001
                    b -= db*0.0001
                    g -= dg*0.0001
                    print(sum(exprs).subs(alpha,a).subs(beta,b).subs(gamma,g),a,b,g,da,db,dg)
                sys.exit()
                '''
                '''
                #beta*(alpha*triangle[0][i]+(1-alpha)*triangle[1][i])+(1-beta)*triangle[2][i]
                a, b, g, distance = Polyhedron.pinpoint(triangle, position, direction)
                if distance < 0.5:
                    min_distance = distance
                    point = tuple(a*triangle[0][i]+(1-a)*triangle[1][i] for i in range(3))
                    point = tuple(b*x+(1-b)*triangle[2][i] for i,x in enumerate(point))
                    output.append((g, point, face))
                    double_break = True
                '''
                func = lambda x: [(x[0]*triangle[0][i]+(1-x[0])*triangle[1][i])*x[1]+(1-x[1])*triangle[2][i]-(position[i]+direction[i]*x[2]) for i in range(3)]
                solution = fsolve(func, (0,0,0))
                alpha, beta, gamma = solution
                if not np.allclose(func(solution), [0,0,0]):
                    continue
                '''
                try:
                    solutions = nsolve(eqs, (alpha,beta,gamma), (1,1,1), dict=True)
                except ValueError:
                    continue
                alpha, beta, gamma = solutions[0][alpha], solutions[0][beta], solutions[0][gamma]
                '''
                #(beta*(alpha*triangle[0][i]+(1-alpha)*triangle[1][i])+(1-beta)*triangle[2][i]-position[i]+direction[i]*gamma)**2
                if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1:
                    point = tuple(alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3))
                    point = tuple(beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(point))
                    output.append((gamma, point, face_index))
                    break
        return sorted(output)
    def merge(self, other):
        edges1 = [frozenset([self.verts[index] for index in edge]) for edge in self.edges]
        edges2 = [frozenset([other.verts[index] for index in edge]) for edge in other.edges]
        faces1 = set([frozenset([edges1[index] for index in face]) for face in self.faces])
        faces2 = set([frozenset([edges2[index] for index in face]) for face in other.faces])
        edges1 = set(edges1)
        edges2 = set(edges2)
        faces_to_remove = faces1 & faces2
        edges_to_remove = set()
        for face in faces_to_remove:
            edges_to_remove.update(face)
        new_edges = set()
        new_faces = set()
        for e1 in edges1:
            for e2 in edges2:
                if e1 == e2:
                    edges_to_remove.add(e1)
                    temp = []
                elif len(e1&e2) and Polyhedron.colinear(list(e1|e2)):
                    edges_to_remove.add(e1)
                    edges_to_remove.add(e2)
                    new_edges.add(e1^e2)
                elif Polyhedron.colinear(list(e1|e2)):
                    e1 = sorted(tuple(e1))
                    e2 = sorted(tuple(e2))
                    interior_points = set()
                    for index,p in enumerate(e1):
                        a = np.array([[e2[1][i]-e2[0][i]] for i in range(3)])
                        b = np.array([p[i]-e2[0][i] for i in range(3)])
                        x, res, rank, s = np.linalg.lstsq(a, b)
                        if x[0] >= 0 and x[0] <= 1:
                            interior_points.add(p)
                    for index,p in enumerate(e2):
                        a = np.array([[e1[1][i]-e1[0][i]] for i in range(3)])
                        b = np.array([p[i]-e1[0][i] for i in range(3)])
                        x, res, rank, s = np.linalg.lstsq(a, b)
                        if x[0] >= 0 and x[0] <= 1:
                            interior_points.add(p)
                    if not len(set(e1)-interior_points):
                        edges_to_remove.add(frozenset(e1))
                        edges_to_remove.add(frozenset(e2))
                        '''
                        for face in faces1|faces2:
                            if frozenset(e1) in face or frozenset(e2) in face:
                                faces_to_remove.add(face)
                        '''
                        new_edges.add(frozenset([e2[0],e1[0]]))
                        new_edges.add(frozenset([e1[1],e2[1]]))
                    if not len(set(e2)-interior_points):
                        edges_to_remove.add(frozenset(e1))
                        edges_to_remove.add(frozenset(e2))
                        '''
                        for face in faces1|faces2:
                            if frozenset(e1) in face or frozenset(e2) in face:
                                faces_to_remove.add(face)
                        '''
                        new_edges.add(frozenset([e1[0],e2[0]]))
                        new_edges.add(frozenset([e2[1],e1[1]]))
                    e1 = frozenset(e1)


        for f1 in faces1:
            for f2 in faces2:
                if f1 == f2:
                    faces_to_remove.add(f1)
                    edges_to_remove.update(f1)
                elif Polyhedron.coplanar(list({point for edge in f1 for point in edge}|{point for edge in f2 for point in edge})):
                    cont = True
                    for e1 in f1:
                        for e2 in f2:
                            e1 = tuple(e1)
                            e2 = tuple(e2)
                            interior_points = set()
                            for index,p in enumerate(e1):
                                a = np.array([[e2[1][i]-e2[0][i]] for i in range(3)])
                                b = np.array([p[i]-e2[0][i] for i in range(3)])
                                x, res, rank, s = np.linalg.lstsq(a, b)
                                if x[0] >= 0 and x[0] <= 1 and np.allclose(np.dot(a, x), b, atol=0.1):
                                    interior_points.add(p)
                            for index,p in enumerate(e2):
                                a = np.array([[e1[1][i]-e1[0][i]] for i in range(3)])
                                b = np.array([p[i]-e1[0][i] for i in range(3)])
                                x, res, rank, s = np.linalg.lstsq(a, b)
                                if x[0] >= 0 and x[0] <= 1 and np.allclose(np.dot(a, x), b, atol=0.1):
                                    interior_points.add(p)
                            if len(interior_points):
                                #print('interior', interior_points)
                                cont = False
                    if cont:
                        continue
                    edges_to_remove.update(f1&f2)
                    faces_to_remove.add(f1)
                    faces_to_remove.add(f2)
                    f3 = f1 ^ f2
                    for e1 in f1:
                        for e2 in f2:
                            if len(e1&e2) == 1 and Polyhedron.colinear(list(e1|e2)):
                                f3 |= frozenset([e1^e2])
                                f3 -= frozenset([e1,e2])
                    for e1 in f3 & edges_to_remove:
                        for e2 in new_edges:
                            if len(e1&e2) == 1 and Polyhedron.colinear(list(e1|e2)):
                                f3 |= frozenset([e2])
                    f3 = f3 - edges_to_remove
                    visited = set()
                    try:
                        current = list(f3)[0]
                    except IndexError:
                        continue
                    while True:
                        visited.add(current)
                        for edge in f3:
                            if len(current&edge) == 1 and edge not in visited:
                                current = edge
                                break
                        else:
                            break
                    if len(visited) != len(f3):
                        #print(f1, f2)
                        #print(f1&f2)
                        #print(store)
                        new_faces.add(frozenset(visited))
                        new_faces.add(frozenset(f3-visited))
                    else:
                        new_faces.add(f3)
        temp = list(new_edges) 
        for e1 in temp:
            for e2 in temp:
                if e1 != e2 and Polyhedron.colinear(list(e1|e2)):
                    e1 = tuple(e1)
                    e2 = tuple(e2)
                    interior_points = set()
                    for index,p in enumerate(e1):
                        a = np.array([[e2[1][i]-e2[0][i]] for i in range(3)])
                        b = np.array([p[i]-e2[0][i] for i in range(3)])
                        x, res, rank, s = np.linalg.lstsq(a, b)
                        if x[0] >= 0 and x[0] <= 1:
                            interior_points.add(p)
                    for index,p in enumerate(e2):
                        a = np.array([[e1[1][i]-e1[0][i]] for i in range(3)])
                        b = np.array([p[i]-e1[0][i] for i in range(3)])
                        x, res, rank, s = np.linalg.lstsq(a, b)
                        if x[0] >= 0 and x[0] <= 1:
                            interior_points.add(p)
                    if len(interior_points):
                        e1 = frozenset(e1)
                        e2 = frozenset(e2)
                        try:
                            new_edges.remove(e1)
                        except KeyError:
                            pass
                        try:
                            new_edges.remove(e2)
                        except KeyError:
                            pass
                        e3 = (e1|e2)-interior_points
                        #print(e1,e2,e3,interior_points)
                        new_edges.add(e3)
                        new_faces = list(new_faces)
                        for i in range(len(new_faces)):
                            if e1 in new_faces[i]:
                                new_faces[i] = new_faces[i] - frozenset([e1])
                                new_faces[i] = new_faces[i]|frozenset([e3])
                            if e2 in new_faces[i]:
                                new_faces[i] = new_faces[i] - frozenset([e2])
                                new_faces[i] = new_faces[i]|frozenset([e3])
                        new_faces = set(new_faces)
                    e1 = frozenset(e1)

        

        #print('faces to remove/faces 2', faces2&faces_to_remove)
        edges = list((edges1|edges2|new_edges)-edges_to_remove)
        verts = list(set(self.verts)^set(other.verts))
        faces = (faces1|faces2|new_faces)-faces_to_remove
        faces = [frozenset([edges.index(edge) for edge in face]) for face in faces]
        edges = [frozenset([verts.index(point) for point in edge]) for edge in edges]
        #print(verts)
        #print(edges)
        #print(faces)
        #polyhedron = Polyhedron()
        self.verts = verts
        self.edges = edges
        self.faces = faces
        #return polyhedron
        #sys.exit()

def get_cube(displacement=(0,0,0), factors=(1,1,1), angles=(0,0,0)):
    points = [(0,0,0),(1,0,0),(1,1,0),(0,1,0)]
    points.extend([(p[0],p[1],1) for p in points])
    center = (0.5,0.5,0.5)
    points = shift(points, (-center[0], -center[1], -center[2]))
    points = scale(points, factors)
    points = rotate(points, angles)
    #points = shift(points, center)
    points = shift(points, displacement)
    edges = [(i,(i+1)%4) for i in range(4)]
    edges.extend([(4+i,4+(i+1)%4) for i in range(4)])
    edges.extend([(i,i+4) for i in range(4)])
    polyhedron = Polyhedron()
    polyhedron.verts = points
    polyhedron.edges = [frozenset(x) for x in edges]
    polyhedron.construct_faces()
    return polyhedron

class Crosshair():
    def __init__(self):
        self.pos = (0,0)
    def draw(self, pygame, screen):
        screen_width, screen_height = screen.get_size()
        p = [(self.pos[0]+screen_width/2+x, -self.pos[1]+screen_height/2) for x in (-10,10)]
        pygame.draw.line(screen, "white", p[0], p[1])
        p = [(self.pos[0]+screen_width/2, -self.pos[1]+screen_height/2+x) for x in (-10,10)]
        pygame.draw.line(screen, "white", p[0], p[1])

class Block:
    def __init__(self, width, height, depth, unit=1):
        self.size = (width,height,depth)
        self.block = get_cube((width/2,height/2,depth/2),factors=(width,height,depth))
        self.illumination_map = dict()
        self.contig = dict()
        self.polys = dict()
        self.segments = dict()
        self.unique = []
        self.select = [0,0,0,0]
        self.unit = unit
        self.select_size = [unit,unit,unit]
    def draw(self, pygame, screen):
        screen_width, screen_height = screen.get_size()
        for edge in self.block.edges:
            p1, p2 = tuple(self.block.verts[index] for index in edge)
            #print(p1,p2)
            p1, p2 = camera.project(p1), camera.project(p2)
            p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
            p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
            pygame.draw.line(screen, "white", p1, p2)
        for poly in self.polys.values():
            for face_index,face in enumerate(poly.faces):
                points = poly.circuit(face_index)
                centroid = tuple(sum(points[i][j] for i,x in enumerate(points))/len(points) for j in range(3))
                distance = sum((centroid[i]-camera.focal[i])**2 for i in range(3))
                points = [camera.project(x) for x in points]
                points = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in points]
                color = "white"
                pygame.draw.polygon(screen, color, points)
            for edge in poly.edges:
                p1, p2 = tuple(poly.verts[index] for index in edge)
                centroid = ((p1[0]+p2[0])/2,(p1[1]+p2[1])/2,(p1[2]+p2[2])/2)
                distance = sum((centroid[i]-camera.focal[i])**2 for i in range(3))
                #print(p1,p2)
                p1, p2 = camera.project(p1), camera.project(p2)
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                pygame.draw.line(screen, "black", p1, p2)
        width, height, depth = self.size
        if self.select[3] == 0:
            axes = get_cube((width/2,self.unit/2+self.select[1],depth/2), factors=(width,self.unit,depth))
        elif self.select[3] == 1:
            axes = get_cube((width/2,height/2,self.unit/2+self.select[2]), factors=(width,height,self.unit))
        elif self.select[3] == 2:
            axes = get_cube((self.unit/2+self.select[0],height/2,depth/2), factors=(self.unit,height,depth))
        for edge in axes.edges:
            p1, p2 = tuple(axes.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                pygame.draw.line(screen, "white", p1, p2)
        select_cube = get_cube((self.select[0]+self.select_size[0]/2,self.select[1]+self.select_size[1]/2,self.select[2]+self.select_size[2]/2), factors=self.select_size)
        for edge in select_cube.edges:
            p1, p2 = tuple(select_cube.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                pygame.draw.line(screen, (128,128,128), p1, p2)
        if self.select[1] > 0:
            shadow = get_cube((self.select[0]+self.select_size[0]/2,0+self.select[1]/2,self.select[2]+self.select_size[2]/2), factors=(self.select_size[0],self.select[1],self.select_size[2]))
            for edge in shadow.edges:
                p1, p2 = tuple(shadow.verts[index] for index in edge)
                p1, p2 = camera.project(p1), camera.project(p2)
                if p1 and p2:
                    p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                    p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                    pygame.draw.line(screen, (128,128,128), p1, p2)
        if self.select[0] > 0:
            shadow = get_cube((0+self.select[0]/2,self.select[1]+self.select_size[1]/2,self.select[2]+self.select_size[2]/2), factors=(self.select[0],self.select_size[1],self.select_size[2]))
            for edge in shadow.edges:
                p1, p2 = tuple(shadow.verts[index] for index in edge)
                p1, p2 = camera.project(p1), camera.project(p2)
                if p1 and p2:
                    p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                    p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                    pygame.draw.line(screen, (128,128,128), p1, p2)
    def add_poly(self,pos,size=(1,1,1),reconstruct=True):
        edited = set()
        i = pos[0]
        while i < pos[0]+size[0]:
            j = pos[1]
            while j < pos[1]+size[1]:
                k = pos[2]
                while k < pos[2]+size[2]:
                    self.contig[(i,j,k)] = {(i,j,k)}
                    for delta in [(self.unit,0,0),(0,self.unit,0),(0,0,self.unit)]:
                        for m in [-1,1]:
                            other = (i+m*delta[0],j+m*delta[1],k+m*delta[2])
                            if other in self.contig and self.contig[(i,j,k)] != self.contig[other]:
                                edited.add(frozenset(self.contig[other]))
                                self.contig[other].update(self.contig[(i,j,k)])
                                self.contig[(i,j,k)] = self.contig[other]
                                for other in self.contig[(i,j,k)]:
                                    self.contig[other] = self.contig[(i,j,k)]
                    k += self.unit
                j += self.unit
            i += self.unit
        i,j,k = pos
        contig = self.contig[(i,j,k)]
        segments = dict()
        for x in edited:
            if x in self.polys:
                del self.polys[x]
            if x in self.segments:
                for key in self.segments[x]:
                    if key not in segments:
                        segments[key] = 0
                    segments[key] += self.segments[x][key]
        for i1 in range(3):
            pos = [i,j,k]
            start = [x for x in pos]
            while pos[i1] < start[i1]+size[i1]:
                end = [x for x in pos]
                end[i1] += self.unit
                segment = frozenset([tuple(pos),tuple(end)])
                if segment not in segments:
                    segments[segment] = 0
                segments[segment] += 1
                pos[i1] += self.unit
            for i2 in range(3):
                if i1 < i2:
                    pos = [i,j,k]
                    pos[i1] += size[i1]
                    start = [x for x in pos]
                    while pos[i2] < start[i2]+size[i2]:
                        end = [x for x in pos]
                        end[i2] += self.unit
                        segment = frozenset([tuple(pos),tuple(end)])
                        if segment not in segments:
                            segments[segment] = 0
                        segments[segment] += 1
                        pos[i2] += self.unit
                    pos = [i,j,k]
                    pos[i2] += size[i2]
                    start = [x for x in pos]
                    while pos[i1] < start[i1]+size[i1]:
                        end = [x for x in pos]
                        end[i1] += self.unit
                        segment = frozenset([tuple(pos),tuple(end)])
                        if segment not in segments:
                            segments[segment] = 0
                        segments[segment] += 1
                        pos[i1] += self.unit
                    for i3 in range(3):
                        if i1 != i3 and i2 != i3:
                            pos = [i,j,k]
                            pos[i1] += size[i1]
                            pos[i2] += size[i2]
                            start = [x for x in pos]
                            while pos[i3] < start[i3]+size[i3]:
                                end = [x for x in pos]
                                end[i3] += self.unit
                                segment = frozenset([tuple(pos),tuple(end)])
                                if segment not in segments:
                                    segments[segment] = 0
                                segments[segment] += 1
                                pos[i3] += self.unit
        '''
        for d1 in [(1,0,0),(0,1,0),(0,0,1)]:
            segment = frozenset([(i,j,k),(i+d1[0],j+d1[1],k+d1[2])])
            if segment not in segments:
                segments[segment] = 0
            segments[segment] += 1
            for d2 in [(1,0,0),(0,1,0),(0,0,1)]:
                if d1 < d2:
                    segment = frozenset([(i+d1[0],j+d1[1],k+d1[2]),(i+d1[0]+d2[0],j+d1[1]+d2[1],k+d1[2]+d2[2])])
                    if segment not in segments:
                        segments[segment] = 0
                    segments[segment] += 1
                    segment = frozenset([(i+d2[0],j+d2[1],k+d2[2]),(i+d1[0]+d2[0],j+d1[1]+d2[1],k+d1[2]+d2[2])])
                    if segment not in segments:
                        segments[segment] = 0
                    segments[segment] += 1
                    for d3 in [(1,0,0),(0,1,0),(0,0,1)]:
                        if d1 != d3 and d2 != d3:
                            segment = frozenset([(i+d1[0]+d2[0],j+d1[1]+d2[1],k+d1[2]+d2[2]),(i+d1[0]+d2[0]+d3[0],j+d1[1]+d2[1]+d3[1],k+d1[2]+d2[2]+d3[2])])
                            if segment not in segments:
                                segments[segment] = 0
                            segments[segment] += 1
        '''
        self.segments[frozenset(contig)] = segments
        if reconstruct:
            self.construct_poly(frozenset(contig))
    def construct_poly(self, contig):
        segments = self.segments[contig]
        segment_map = dict()
        for segment in segments:
            if segments[segment] % 2 == 1:
                for point in segment:
                    if point not in segment_map:
                        segment_map[point] = set()
                    segment_map[point].add(segment)

        for key in segment_map:
            found = True
            while found:
                found = False
                for segment1 in segment_map[key]:
                    for segment2 in segment_map[key]:
                        if segment1 != segment2 and Polyhedron.colinear(list(segment1|segment2)):
                            segment3 = segment1 ^ segment2
                            found = True
                            break
                    if found:
                        break
                if found:
                    for point in segment1:
                        segment_map[point].remove(segment1)
                    for point in segment2:
                        segment_map[point].remove(segment2)
                    for point in segment3:
                        segment_map[point].add(segment3)
        edges = set()
        for key in list(segment_map.keys()):
            edges.update(segment_map[key])
            if len(segment_map[key]) == 0:
                del segment_map[key]
        verts = set()
        for edge in edges:
            verts.update(edge)
        face_map = dict()
        for key in segment_map:
            for edge1 in segment_map[key]:
                for edge2 in segment_map[key]:
                    if edge1 != edge2 and len(edge1&edge2):
                        if edge1 not in face_map:
                            face_map[edge1] = set()
                        face_map[edge1].add(frozenset([edge1,edge2]))
                        if edge2 not in face_map:
                            face_map[edge2] = set()
                        face_map[edge2].add(frozenset([edge1,edge2]))
        for key in face_map:
            found = True
            while found:
                found = False
                for face1 in face_map[key]:
                    for face2 in face_map[key]:
                        if face1 != face2 and Polyhedron.coplanar(list({point for edge in face1|face2 for point in edge})):
                            found = True
                            face3 = face1|face2
                            break
                    if found:
                        break
                if found:
                    for edge in face1:
                        face_map[edge].remove(face1)
                    for edge in face2:
                        face_map[edge].remove(face2) 
                    for edge in face3:
                        face_map[edge].add(face3)

        faces = set()
        for value in face_map.values():
            faces.update(value)
        for face in faces:
            if len(face) == 3:
                print(face)
                sys.exit()
        verts = list(verts)
        edges = list(edges)
        poly = Polyhedron()
        poly.faces = [frozenset([edges.index(edge) for edge in face]) for face in faces]
        poly.edges = [frozenset([verts.index(vert) for vert in edge]) for edge in edges]
        poly.verts = verts
        self.polys[contig] = poly
    def del_poly(self,pos,size=(1,1,1)):
        i,j,k = pos
        contig = self.contig[(i,j,k)]
        del self.polys[frozenset(contig)]
        removed = set()
        i = pos[0]
        while i < pos[0]+size[0]:
            j = pos[1]
            while j < pos[1]+size[1]:
                k = pos[2]
                while k < pos[2]+size[2]:
                    del self.contig[(i,j,k)]
                    contig.remove((i,j,k))
                    removed.add((i,j,k))
                    k += self.unit
                j += self.unit
            i += self.unit
        visited = set()
        queue = []
        triple_break = False
        for i,j,k in removed:
            for delta in [(self.unit,0,0),(0,self.unit,0),(0,0,self.unit)]:
                for m in [-1,1]:
                    other = (i+m*delta[0],j+m*delta[1],k+m*delta[2])
                    if other in contig:
                        queue.append(other)
                        triple_break = True
                        break
                if triple_break:
                    break
            if triple_break:
                break
        while len(queue):
            i,j,k = queue.pop()
            visited.add((i,j,k))
            for delta in [(self.unit,0,0),(0,self.unit,0),(0,0,self.unit)]:
                for m in [-1,1]:
                    other = (i+m*delta[0],j+m*delta[1],k+m*delta[2])
                    if other in contig and other not in visited:
                        queue.append(other)
        if len(visited) == len(contig):
            segments = self.segments[frozenset(contig|removed)]
            del self.segments[frozenset(contig|removed)]
            for i,j,k in removed:
                for d1 in [(self.unit,0,0),(0,self.unit,0),(0,0,self.unit)]:
                    segment = frozenset([(i,j,k),(i+d1[0],j+d1[1],k+d1[2])])
                    if segment not in segments:
                        segments[segment] = 0
                    segments[segment] -= 1
                    for d2 in [(self.unit,0,0),(0,self.unit,0),(0,0,self.unit)]:
                        if d1 < d2:
                            segment = frozenset([(i+d1[0],j+d1[1],k+d1[2]),(i+d1[0]+d2[0],j+d1[1]+d2[1],k+d1[2]+d2[2])])
                            if segment not in segments:
                                segments[segment] = 0
                            segments[segment] -= 1
                            segment = frozenset([(i+d2[0],j+d2[1],k+d2[2]),(i+d1[0]+d2[0],j+d1[1]+d2[1],k+d1[2]+d2[2])])
                            if segment not in segments:
                                segments[segment] = 0
                            segments[segment] += 1
                            for d3 in [(self.unit,0,0),(0,self.unit,0),(0,0,self.unit)]:
                                if d1 != d3 and d2 != d3:
                                    segment = frozenset([(i+d1[0]+d2[0],j+d1[1]+d2[1],k+d1[2]+d2[2]),(i+d1[0]+d2[0]+d3[0],j+d1[1]+d2[1]+d3[1],k+d1[2]+d2[2]+d3[2])])
                                    if segment not in segments:
                                        segments[segment] = 0
                                    segments[segment] -= 1
            self.segments[frozenset(contig)] = segments
            self.construct_poly(frozenset(contig))
        else:      
            del self.segments[frozenset(contig|removed)]
            for other in contig:
                del self.contig[other]
            edited = set()
            for i,j,k in contig:
                self.add_poly((i,j,k),reconstruct=False)
                edited.add((i,j,k))
            for contig in {frozenset(self.contig[x]) for x in edited}:
                self.construct_poly(contig)
    def subdivide(self, divisor):
        block = Block(self.size[0], self.size[1], self.size[2], self.unit/divisor)
        block.contig = self.contig
        block.select = self.select
        block.polys = self.polys 
        block.segments = dict()
        for key in self.segments:
            segments = dict()
            for segment in self.segments[key]:
                low, high = min(segment), max(segment)
                diff = ((high[0]-low[0])/divisor,(high[1]-low[1])/divisor,(high[2]-low[2])/divisor)
                point = low
                while point < high:
                    seg = frozenset([point, tuple(point[i]+diff[i] for i in range(3))])
                    segments[seg] = self.segments[key][segment]
                    point = tuple(point[i]+diff[i] for i in range(3))
            block.segments[key] = segments
        points = list(self.contig.keys())
        for point in points:
            contig = frozenset(self.contig[point])
            i = point[0]
            while i < point[0]+self.unit:
                j = point[1]
                while j < point[1]+self.unit:
                    k = point[2]
                    while k < point[2]+self.unit:
                        if (i,j,k) not in block.contig[point]:
                            block.contig[point].add((i,j,k))
                            block.contig[(i,j,k)] = block.contig[point]
                            if contig in block.segments:
                                new_contig = frozenset(block.contig[point])
                                block.segments[new_contig] = block.segments[contig]
                                block.polys[new_contig] = block.polys[contig]
                                del block.segments[contig]
                                del block.polys[contig]
                                contig = new_contig
                        k += block.unit
                    j += block.unit
                i += block.unit

        '''
        for key in self.contig:
            block.add_poly(key,tuple(self.unit for i in range(3)),reconstruct=False)
        for value in {frozenset(x) for x in block.contig.values()}:
            block.construct_poly(value)
        '''
        return block
    def select_by_void(self):
        i,j,k = self.select[:3]
        if not (i,j,k) in self.contig:
            return True
        for delta in [(self.unit,0,0),(0,self.unit,0),(0,0,self.unit)]:
            for m in [-1,1]:
                other = (i+m*delta[0],j+m*delta[1],k+m*delta[2])
                if not other in self.contig:
                    return True
        return False
        

    def contiguous(self):
        output = {(i,j,k):frozenset([(i,j,k)]) for (i,j,k) in self.polys}
        for i,j,k in output:
            for delta in [(1,0,0),(0,1,0),(0,0,1)]:
                other = (i+delta[0],j+delta[1],k+delta[2])
                if other in output:
                    output[(i,j,k)] |= output[other]
                    output[other] = output[(i,j,k)]
                other = (i-delta[0],j-delta[1],k-delta[2])
                if other in output:
                    output[(i,j,k)] |= output[other]
                    output[other] = output[(i,j,k)]
        return set(output.values())

    def ray_process(cube, position, direction):
        return cube.ray(position,direction)
    def ray(self, position, direction):
        hits = []
        inputs = []
        coords = []
        for i,j,k in self.polys:
            cube = self.polys[(i,j,k)]
            inputs.append((cube, position, direction))
            coords.append((i,j,k))
        '''
        with mp.Pool(8) as p:
            hits = p.starmap(Block.ray_process, inputs)
            hits = [x[0]+(coords[i],) for i,x in enumerate(hits) if len(x)]
        '''
        hits = [Block.ray_process(*args) for args in inputs]
        hits = [x[0]+(coords[i],) for i,x in enumerate(hits) if len(x)]
        hits = sorted(hits)
        #print(hits)
        return hits
    def illuminate(self, light, coord, face_index):
        i,j,k = coord
        cube = self.polys[(i,j,k)]
        for triangle in Polyhedron.triangulation(cube.circuit(face_index)):
            triangle = list(triangle)
            centroid = tuple((triangle[0][i]+triangle[1][i]+triangle[2][i])/3 for i in range(3))
            vec = tuple(centroid[i]-light[i] for i in range(3))
            mag = sqrt(dot(vec,vec))
            unit = (vec[0]/mag,vec[1]/mag,vec[2]/mag)
            hits = self.ray(light, unit)
            return hits[0]
                        
def lighting_helper(f,x):
    return f(*x)
def process_lighting(input_queue, output_queue):
    while True:
        to_map = []
        while True:
            try:
                f, args = input_queue.get(False)
                print(args)
                to_map.append((f,args))
            except queue.Empty:
                break
        with mp.Pool(8) as p:
            hits = p.starmap(lighting_helper, to_map)
            for hit in hits:
                output_queue.put(hit)
        
if __name__ == "__main__":
    import pygame
    camera = Camera()
    camera.zoom = 100
    camera.focal = (0,0,-100)
    cube = get_cube((0,0,-1))
    pygame.init()
    screen = pygame.display.set_mode((1280, 720))
    screen_width, screen_height = screen.get_size()
    clock = pygame.time.Clock()
    pygame.mouse.set_visible(False)
    pygame.mouse.set_pos((screen_width/2, screen_height/2))
    crosshair = Crosshair()
    mouse_down = False
    running = True
    framerate = 60
    light = (10,10,10)
    cube.triangles = {face_index:Polyhedron.triangulation(cube.circuit(face_index)) for face_index,face in enumerate(cube.faces)}
    light_input_queue = mp.Queue()
    light_output_queue = mp.Queue()
    #light_process = mp.Process(target=process_lighting, args=(light_input_queue,light_output_queue)) 
    #light_process.start()
    light_dict = dict()
    block = Block(3,3,3)
    #block = Block(5,5,5)
    #block.block[1][4][2] = 1
    dts = []
    count = 0
    #block.add_poly(0,0,0)
    i,j,k = 0,0,0
    space_down = False
    dir_mult1 = 1
    dir_mult2 = 1
    z_forward = True
    while running:
        '''
        i,j,k = (random.randrange(block.size[i]) for i in range(3))
        while (i,j,k) in block.polys:
            i,j,k = (random.randrange(block.size[i]) for i in range(3))
        i = count//height//depth%width
        j = count//depth%height
        k = count%depth
        if count < 100:
            width, height, depth = block.size
            delta = random.choice([(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0)])
            new_i = i + delta[0]
            new_j = j + delta[1]
            new_k = k + delta[2]
            while new_i < 0 or new_i >= width or new_j < 0 or new_j >= height or new_k < 0 or new_k >= depth or (new_i,new_j,new_k) in block.contig:
                delta = random.choice([(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0)])
                new_i = i + delta[0]
                new_j = j + delta[1]
                new_k = k + delta[2]
            i,j,k = new_i,new_j,new_k
            print(i,j,k)
            block.add_poly(i,j,k)
        count += 1
        '''
        #print(block.contig)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_a:
                    camera.origin, camera.focal, camera.x_vector, camera.y_vector = rotate([camera.origin, camera.focal, camera.x_vector, camera.y_vector], (0,pi/180*10,0))
                if event.key == pygame.K_d:
                    camera.origin, camera.focal, camera.x_vector, camera.y_vector = rotate([camera.origin, camera.focal, camera.x_vector, camera.y_vector], (0,-pi/180*10,0))
                if event.key == pygame.K_w:
                    camera.origin, camera.focal, camera.x_vector, camera.y_vector = rotate([camera.origin, camera.focal, camera.x_vector, camera.y_vector], (pi/180*10,0,0))
                if event.key == pygame.K_s:
                    camera.origin, camera.focal, camera.x_vector, camera.y_vector = rotate([camera.origin, camera.focal, camera.x_vector, camera.y_vector], (-pi/180*10,0,0))
                if event.key == pygame.K_q:
                    camera.origin, camera.focal, camera.x_vector, camera.y_vector = rotate([camera.origin, camera.focal, camera.x_vector, camera.y_vector], (0,0,-pi/180*10))
                if event.key == pygame.K_e:
                    camera.origin, camera.focal, camera.x_vector, camera.y_vector = rotate([camera.origin, camera.focal, camera.x_vector, camera.y_vector], (0,0,pi/180*10))
                if event.key == pygame.K_2:
                    block = block.subdivide(2)
                if event.key == pygame.K_3:
                    block = block.subdivide(3)
                dir_mult1 = 1
                if dot(camera.forward_vector(),(0,0,1)) < dot(camera.forward_vector(),(0,0,-1)):
                    dir_mult1 = -1
                dir_mult2 = 1
                if dot(camera.forward_vector(),(1,0,0)) < dot(camera.forward_vector(),(-1,0,0)):
                    dir_mult2 = -1
                dir_mult3 = 1
                if (dir_mult1 == -1 or dir_mult2 == -1) and dir_mult1 != dir_mult2:
                    dir_mult3 = -1
                direction = 1
                if event.key in {pygame.K_DOWN, pygame.K_LEFT}:
                    direction = -1
                z_forward = True
                if max(dot(camera.forward_vector(),(0,0,1)),dot(camera.forward_vector(),(0,0,-1))) < max(dot(camera.forward_vector(),(1,0,0)),dot(camera.forward_vector(),(-1,0,0))):
                    z_forward = False
                if event.key in {pygame.K_UP, pygame.K_DOWN}:
                    index = 2
                    if not z_forward:
                        index = 0
                    if not space_down:
                        if block.select[3] == 0:
                            block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                            while not block.select_by_void():
                                block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                        elif block.select[3] == 1:
                            block.select[1] += direction*block.unit
                            while not block.select_by_void():
                                block.select[1] += direction*block.unit
                        elif block.select[3] == 2:
                            block.select[1] += direction*block.unit
                            while not block.select_by_void():
                                block.select[1] += direction*block.unit
                        for i in range(3):
                            if block.select[i] < 0:
                                block.select[i] = 0
                            if block.select[i] > block.size[i]-block.unit:
                                block.select[i] = block.size[i]-block.unit
                    else:
                        if block.select[3] == 0:
                            block.select_size[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                        elif block.select[3] == 1:
                            block.select_size[1] += direction*block.unit
                        elif block.select[3] == 2:
                            block.select_size[1] += direction*block.unit
                        if block.select_size[index] == 0:
                            block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                        if block.select_size[1] == 0:
                            block.select_size[1] = direction*block.unit
                        for i in range(3):
                            if block.select[i]+block.select_size[i] > block.size[i]:
                                block.select_size[i] = block.size[i]-block.select[i]
                            if block.select[i]+block.select_size[i] < 0:
                                block.select_size[i] = -block.select[i]
                if event.key in {pygame.K_RIGHT, pygame.K_LEFT}:
                    index = 0
                    if not z_forward:
                        index = 2
                    if not space_down:
                        if block.select[3] == 0:
                            block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(2*int(z_forward)-1)*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                            while not block.select_by_void():
                                block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(2*int(z_forward)-1)*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                        elif block.select[3] == 1:
                            block.select[0] += direction*dir_mult1*block.unit
                            while not block.select_by_void():
                                block.select[0] += direction*dir_mult1*block.unit
                        elif block.select[3] == 2:
                            block.select[2] += direction*dir_mult2*block.unit
                            while not block.select_by_void():
                                block.select[2] += direction*dir_mult2*block.unit
                        for i in [0,2]:
                            if block.select[i] > block.size[i]-block.unit:
                                block.select[i] = block.size[i]-block.unit
                            if block.select[i] < 0:
                                block.select[i] = 0
                    else:
                        if block.select[3] == 0:
                            block.select_size[index] += direction*dir_mult1*dir_mult3*block.unit*(2*int(z_forward)-1)*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                            if block.select_size[index] == 0:
                                block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(2*int(z_forward)-1)*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                        elif block.select[3] == 1:
                            block.select_size[0] += direction*dir_mult1*block.unit
                            if block.select_size[0] == 0:
                                block.select_size[0] = direction*dir_mult1*block.unit
                        elif block.select[3] == 2:
                            block.select_size[2] -= direction*dir_mult2*block.unit
                            if block.select_size[2] == 0:
                                block.select_size[2] = -direction*dir_mult2*block.unit
                        for i in [0,2]:
                            if block.select[i]+block.select_size[i] > block.size[i]:
                                block.select_size[i] = block.size[i]-block.select[i]
                            if block.select[i]+block.select_size[i] < 0:
                                block.select_size[i] = -block.select[i]
                if event.key == pygame.K_LSHIFT:
                    block.select[3] = (block.select[3]-1)%3
                if event.key == pygame.K_RSHIFT:
                    block.select[3] = (block.select[3]+1)%3
                if event.key == pygame.K_z:
                    camera.zoom *= 1.1
                if event.key == pygame.K_x:
                    camera.zoom /= 1.1
                if event.key == pygame.K_SPACE:
                    space_down = True    
            if event.type == pygame.KEYUP:
                if event.key == pygame.K_SPACE:
                    space_down = False
                    to_add = False
                    i = block.select[0]
                    i,i_max = block.select[0],block.select[0]+block.select_size[0]
                    if block.select_size[0] < 0:
                        i,i_max = i_max,i
                    while i < i_max:
                        j,j_max = block.select[1],block.select[1]+block.select_size[1]
                        if block.select_size[1] < 0:
                            j,j_max = j_max,j
                        while j < j_max:
                            k,k_max = block.select[2],block.select[2]+block.select_size[2]
                            if block.select_size[2] < 0:
                                k,k_max = k_max,k
                            while k < k_max:
                                if not (i,j,k) in block.contig:
                                    to_add = True
                                k += block.unit
                            j += block.unit
                        i += block.unit
                    i = min(block.select[0],block.select[0]+block.select_size[0])
                    j = min(block.select[1],block.select[1]+block.select_size[1])
                    k = min(block.select[2],block.select[2]+block.select_size[2])
                    size = (abs(block.select_size[0]),abs(block.select_size[1]),abs(block.select_size[2]))
                    if to_add:
                        block.add_poly((i,j,k),size)
                    else:
                        block.del_poly((i,j,k),size)
                    block.select = [block.select[0]+block.select_size[0]-block.unit,block.select[1]+block.select_size[1]-block.unit,block.select[2]+block.select_size[2]-block.unit,block.select[3]]
                    block.select_size = [block.unit,block.unit,block.unit]
            if event.type == pygame.MOUSEMOTION:
                x,y = pygame.mouse.get_pos()
                x -= screen_width/2
                y = -(y-screen_height/2)
                crosshair.pos = (x,y)
                if mouse_down:
                    x,y = pygame.mouse.get_pos()
                    delta_x, delta_y = x-mousedown_x, y-mousedown_y
                    camera.move((camera.x_vector[0]*-delta_x/camera.zoom,camera.x_vector[1]*-delta_x/camera.zoom,camera.x_vector[2]*-delta_x/camera.zoom))
                    camera.move((camera.y_vector[0]*delta_y/camera.zoom,camera.y_vector[1]*delta_y/camera.zoom,camera.y_vector[2]*delta_y/camera.zoom))
                    mousedown_x, mousedown_y = x,y
            if event.type == pygame.MOUSEBUTTONDOWN:
                mousedown_x,mousedown_y = pygame.mouse.get_pos() 
                mouse_down = True
            if event.type == pygame.MOUSEBUTTONUP:
                mouse_down = False

        screen.fill("black")
        crosshair.draw(pygame, screen)
        for edge in cube.edges:
            p1, p2 = tuple(cube.verts[index] for index in edge)
            #print(p1,p2)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                pygame.draw.line(screen, "white", p1, p2)
        block.draw(pygame,screen)
        '''
        for face_index,face in enumerate(cube.faces):
            for triangle in cube.triangles[face_index]:
                triangle = list(triangle)
                centroid = tuple((triangle[0][i]+triangle[1][i]+triangle[2][i])/3 for i in range(3))
                vec = tuple(centroid[i]-light[i] for i in range(3))
                mag = sqrt(dot(vec,vec))
                unit = (vec[0]/mag,vec[1]/mag,vec[2]/mag)
                hits = cube.ray(light, unit)

                if hits[0][2] == face:
                    triangle = [camera.project(x) for x in triangle]
                    triangle = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in triangle]
                    print(hits[0][0])
                    pygame.draw.polygon(screen, (255*100/hits[0][0]**2,255*100/hits[0][0]**2,255*100/hits[0][0]**2), triangle)
        '''
        #block.illuminate(light)
        while True:
            try:
                hit = light_output_queue.get(False)
                print(hit)
                block.illumination_map[hit[3]][hit[2]] = 1/hit[0]**2
            except queue.Empty:
                break
        for key in light_dict:
            triangle = [camera.project(x) for x in key]
            triangle = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in triangle]
            pygame.draw.polygon(screen, light_dict[key], triangle)
        pygame.display.flip()
        dT = clock.tick(60) / 1000
        print(dT, len(block.polys), dir_mult1, dir_mult2, z_forward)
        dts.append(dT)
    #light_process.kill()
    plt.plot(dts)
    plt.ylabel('time deltas')
    plt.show()
    pygame.quit()
