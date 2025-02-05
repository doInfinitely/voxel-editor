from sympy import nsolve, solve, symbols, Eq, diff
import sympy
from math import sin, cos, pi, sqrt, copysign
import numpy as np
from decimal import Decimal
import random
from scipy.optimize import fsolve
import time
from multiprocessing import Pool
import math

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
def round_float(x):
    return float(round(x*1000)/1000)
def round_point(point):
    return tuple(round_float(x) for x in point)
def round_edge(edge):
    return frozenset(round_point(x) for x in edge)
def round_float_meter(x, meters=None):
    if meters is None:
        meters = {1}
    mindiff = (float('inf'), None)
    for y in meters:
        diff = abs(x/y - round(x/y))
        if diff < mindiff[0]:
            mindiff = (diff, round(x/y)*y)
    return float(mindiff[1])
def round_point_meter(x, meters=None):
    return tuple(round_float_meter(y, meters) for y in x)
def round_edge_meter(edge, meters=None):
    return frozenset(round_point_meter(x,meters) for x in edge)
def remeter(meters, new_meter, precision=0.001):
    for x in meters:
        if x != 1 and abs(new_meter/x - round(new_meter/x)) <= precision:
            return meters
    return {x for x in meters if abs(x/new_meter - round(x/new_meter)) > precision}|set([new_meter])
def distance(point1, point2):
    displacement = tuple(point2[i]-x for i,x in enumerate(point1))
    return sqrt(sum(x**2 for x in displacement))
def dot(vector1, vector2):
    return sum(x*vector2[i] for i,x in enumerate(vector1))
def cross3D(vector1, vector2):
    return tuple(vector1[(i+1)%3]*vector2[(i+2)%3]-vector1[(i+2)%3]*vector2[(i+1)%3] for i in range(3))

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

class Box:
    def __init__(self, x, y, z):
        self.x_min = x[0]
        self.x_max = x[1]
        self.y_min = y[0]
        self.y_max = y[1]
        self.z_min = z[0]
        self.z_max = z[1]
    def colinear(points):
        if len(points) < 3:
            return True
        points = list(points)
        for p in points[2:]:
            a = np.array([[float(points[1][i]-points[0][i])] for i in range(3)])
            b = np.array([float(p[i]-points[0][i]) for i in range(3)])
            x, res, rank, s = np.linalg.lstsq(a, b)
            if not np.allclose(np.dot(a, x), b, rtol=0.1):
                return False
        return True
    def coplanar(points, rtol=0.1):
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
        if Box.colinear(points):
            return True
        triple_break = False
        for i,x in enumerate(points):
            for j,y in enumerate(points):
                for k,z in enumerate(points):
                    if i != j and i != k and j != k and not Box.colinear([x,y,z]):
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
            x, res, rank, s = np.linalg.lstsq(a.astype('float'), b.astype('float'))
            #print(a, b, x, res)
            if not np.allclose(np.dot(a.astype('float'), x.astype('float')), b.astype('float'), rtol=rtol):
                return False
        return True
    def point_on_segment(edge, point):
        p1, p2 = edge
        alpha = symbols("alpha")
        eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], point[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            alpha = solutions[0][alpha]
            if alpha >= 0 and alpha <= 1:
                return True
        return False
    def intersect_segments(edge1, edge2):
        p1, p2 = edge1
        p3, p4 = edge2
        alpha, beta = symbols("alpha beta")
        eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], beta*p3[i]+(1-beta)*p4[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            try:
                alpha = solutions[0][alpha]
            except KeyError:
                alpha = None
            try:
                beta = solutions[0][beta]
            except KeyError:
                beta = None
            try:
                if alpha is not None and alpha >= 0 and alpha <= 1 and beta is not None and beta >= 0 and beta <= 1:
                    return tuple(alpha*p1[i]+(1-alpha)*p2[i] for i in range(3))
            except TypeError:
                    pass
                    #alpha = 0
                    #output.append(tuple(alpha*p1[i]+(1-alpha)*p2[i] for i in range(3)))
        return None
    def intersect(self, edge):
        p1, p2 = edge
        if p1[0] >= self.x_min and p1[0] <= self.x_max and p1[1] >= self.y_min and p1[1] <= self.y_max and p1[2] >= self.z_min and p1[2] <= self.z_max:
            return True
        if p2[0] >= self.x_min and p2[0] <= self.x_max and p2[1] >= self.y_min and p2[1] <= self.y_max and p2[2] >= self.z_min and p2[2] <= self.z_max:
            return True
        segment_shadow = ((p1[0],0,p1[2]),(p2[0],0,p2[2]))
        intersections = Polyhedron.circuit_intersect(segment_shadow, {((self.x_min,0,self.z_min),(self.x_max,0,self.z_min),(self.x_max,0,self.z_max),(self.x_min,0,self.z_max))})
        for intersection in intersections:
            y1 = p1[1]
            y2 = y1+(p2[1]-p1[1])*intersection[0]
            if y2 >= self.y_min and y2 <= self.y_max:
                return True
        return False

    def inside_polyhedron(poly, point):
        print('inside_polyhedron_start')
        intersections = 0
        for face_index,face in enumerate(poly.faces):
            circuits = poly.circuits(face_index)
            circuit = Polyhedron.circuit_cut(circuits)
            interior_circuits = circuits - set([Polyhedron.find_exterior_circuit(circuits)])
            if any(Polyhedron.inside_triangle(x,point) for x in Polyhedron.triangulate(circuit)):
                return True
            if circuit[0][2] == circuit[1][2] and circuit[0][2] == circuit[2][2] and circuit[0][2] > point[2]:
                projection = (point[0],point[1],circuit[0][2])
                if any(Polyhedron.inside_triangle(x,projection) for x in Polyhedron.triangulate(circuit)) and not any(Polyhedron.inside_triangle(y,projection) for x in interior_circuits for y in Polyhedron.triangulate(x)):
                    print('inside_polyhedron', circuit, interior_circuits)
                    intersections += 1
                if any(Polyhedron.inside_triangle(x,point) for x in Polyhedron.triangulate(circuit)):
                    return True
        print('inside_polyhedron_end')
        return intersections % 2 == 1
    
    def delete_circuit_helper(path_edges, start=None, previous=None, current=None, path=None, seen=None):
        #print(path)
        edge_lookup = dict()
        for edge in path_edges:
            for point in edge:
                if point not in edge_lookup:
                    edge_lookup[round_point(point)] = set()
                edge_lookup[round_point(point)].add(edge)
        if seen is None:
            seen = set()
        circuits = set()
        if start is None:
            #print(path_edges)
            #print('edge_lookup',edge_lookup)
            for edge in path_edges:
                path = []
                start = edge
                point = list(start)[0]
                path.append(point)
                temp = list(edge_lookup[round_point(point)] - set([start]))
                #print(path, temp)
                for i in range(len(temp)):
                    current = temp[i]
                    circuits.update(Box.delete_circuit_helper(path_edges, start, start, current, path, seen))
                del path[-1]
        else:
            #print(current, previous,start)
            if current == start:
                #print(path)
                return {tuple(path)}
            if current in seen:
                return set()
            seen.add(current)
            point = list(current - previous)[0]
            if point in path:
                return {}
            path.append(point)
            previous = current
            temp = list(edge_lookup[round_point(point)]-set([previous]))
            #print(path, temp)
            for i in range(len(temp)):
                if Box.coplanar(set(round_point(x) for x in path)|round_edge(temp[i])):
                    current = temp[i]
                    circuits.update(Box.delete_circuit_helper(path_edges, start, previous, current, path, seen))
            del path[-1]
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            for j,y in enumerate(circuits_list):
                if len(y)==len(x) and not len(set(x)-set(y)) and j > i and y in circuits:
                    circuits.remove(y)
        '''
                if len(y)>len(x) and not len(set(x)-set(y)) and y in circuits:
                    circuits.remove(y)
                if len(y)==len(x) and not len(set(x)-set(y)) and j > i and y in circuits:
                    circuits.remove(y)
                if x != y and not len(set(x)-set(y)):
                    y_r = tuple(reversed(y))
                    if (y[y.index(x[0]):]+y[:y.index(x[0])])[:len(x)] == x or (y_r[y_r.index(x[0]):]+y_r[:y_r.index(x[0])])[:len(x)] == x:
                        if y in circuits:
                            circuits.remove(y)
        '''
        return circuits

    def delete_preprocess(self, poly, face_index):
        start_time_face = time.time()
        face = poly.faces[face_index]
        #print(face)
        circuit = Polyhedron.circuit_cut(poly.circuits(face_index))
        #print('circuit', circuit)
        new_edges = set()
        path_edges = set()
        #print(face)
        #print('circuit', circuit)
        box = []
        box_i = None
        if circuit[0][0] == circuit[1][0] and circuit[0][0] == circuit[2][0] and circuit[0][0] >= self.x_min and circuit[0][0] <= self.x_max:
            box = [(circuit[0][0], self.y_min, self.z_min),(circuit[0][0], self.y_max, self.z_min),(circuit[0][0], self.y_max, self.z_max),(circuit[0][0], self.y_min, self.z_max)]
        elif circuit[0][1] == circuit[1][1] and circuit[0][1] == circuit[2][1] and circuit[0][1] >= self.y_min and circuit[0][1] <= self.y_max:
            box = [(self.x_min, circuit[0][1], self.z_min),(self.x_max, circuit[0][1], self.z_min),(self.x_max, circuit[0][1], self.z_max),(self.x_min, circuit[0][1], self.z_max)]
        elif circuit[0][2] == circuit[1][2] and circuit[0][2] == circuit[2][2] and circuit[0][2] >= self.z_min and circuit[0][2] <= self.z_max:
            box = [(self.x_min, self.y_min, circuit[0][2]),(self.x_max, self.y_min, circuit[0][2]),(self.x_max, self.y_max, circuit[0][2]),(self.x_min, self.y_max, circuit[0][2])]
        if all(x[0] >= self.x_min and x[0] <= self.x_max and x[1] >= self.y_min and x[1] <= self.y_max and x[2] >= self.z_min and x[2] <= self.z_max for i,x in enumerate(circuit)):
            double_break = False
            for p1 in circuit:
                if p1 in poly.verts:
                    for edge in {frozenset(poly.verts[index] for index in edge) for edge in poly.edges}:
                        intersection = 0
                        if p1 in edge:
                            p2 = list(edge - frozenset([p1]))[0]
                            if p2 not in circuit:
                                dim = [p1[i]!=p2[i] for i in range(3)].index(True)
                                vec = tuple(p2[i]-p1[i] for i in range(3))
                                for w in [(self.x_min,self.x_max),(self.y_min,self.y_max),(self.z_min,self.z_max)][dim]:
                                    if dim == 0:
                                        triangles = [[(w,self.y_min,self.z_min),(w,self.y_max,self.z_min),(w,self.y_max,self.z_max)],[(w,self.y_min,self.z_min),(w,self.y_min,self.z_max),(w,self.y_max,self.z_max)]]
                                    if dim == 1:
                                        triangles = [[(self.x_min,w,self.z_min),(self.x_max,w,self.z_min),(self.x_max,w,self.z_max)],[(self.x_min,w,self.z_min),(self.x_min,w,self.z_max),(self.x_max,w,self.z_max)]]
                                    if dim == 2:
                                            triangles = [[(self.x_min,self.y_min,w),(self.x_max,self.y_min,w),(self.x_max,self.y_max,w)],[(self.x_min,self.y_min,w),(self.x_min,self.y_max,w),(self.x_max,self.z_max,w)]]
                                    if vec[dim] > 0:
                                        if w > p1[dim]:
                                            projection = [p1[0],p1[1],p1[2]]
                                            projection[dim] = w
                                            projection = tuple(projection)
                                            if any(Polyhedron.inside_triangle(x,projection) for x in triangles):
                                                intersection += 1
                                    else:
                                        if w < p1[dim]:
                                            projection = [p1[0],p1[1],p1[2]]
                                            projection[dim] = w
                                            projection = tuple(projection)
                                            if any(Polyhedron.inside_triangle(x,projection) for x in triangles):
                                                intersection += 1
                        if intersection%2 == 0:
                            double_break = True
                            break
                    if double_break:
                        break
            else:
                print('face processing time', face_index, time.time()-start_time_face)
                return new_edges, path_edges, box_i
        intersections = []
        for i,x in enumerate(circuit):
            if not any(circuit[i-1] in y and x in y for y in poly.circuits(face_index)):
                continue
            p1 = circuit[i-1]
            p2 = x
            last_intersections = []
            for j, box_point in enumerate(box):
                if not Box.colinear(set([p1,p2,box[j-1],box_point])):
                    intersection = Box.intersect_segments(frozenset([p1,p2]),frozenset([box[j-1],box_point]))
                    #print("intersection", intersection, p1, p2)
                    if intersection is not None:
                        last_intersections.append(intersection)
                else:
                    if Box.point_on_segment(frozenset([box[j-1],box_point]),p1) and Box.point_on_segment(frozenset([box[j-1],box_point]),p2):
                        if p1 not in box:
                            last_intersections.append(p1)
                        if p2 not in box:
                            last_intersections.append(p2)
            #print('last intersections', last_intersections)
            if len(last_intersections) == 0:
                new_edges.add(frozenset([p1,p2]))
            elif len(last_intersections) == 1:
                if (self.x_min > p2[0] or p2[0] > self.x_max) or (self.y_min > p2[1] or p2[1] > self.y_max) or (self.z_min > p2[2] or p2[2] > self.z_max):
                    last_intersections[-1] = (last_intersections[-1],False)
                    #print("intersection p2", p2, last_intersections[0][0])
                    if len(frozenset([p2,last_intersections[0][0]])) == 2:
                        new_edges.add(frozenset([p2,last_intersections[0][0]]))
                    if (last_intersections[0][0], True) in intersections:
                        intersections[intersections.index((last_intersections[0][0], True))] = last_intersections[0]
                    elif (last_intersections[0][0], False) in intersections:
                        intersections[intersections.index((last_intersections[0][0], False))] = last_intersections[0]
                    elif (last_intersections[0][0], None) in intersections:
                        intersections[intersections.index((last_intersections[0][0], None))] = last_intersections[0]
                    else:
                        intersections.extend(last_intersections)
                elif (self.x_min > p1[0] or p1[0] > self.x_max) or (self.y_min > p1[1] or p1[1] > self.y_max) or (self.z_min > p1[2] or p1[2] > self.z_max):
                    last_intersections[-1] = (last_intersections[-1],True)
                    #print("intersection p1", p1, last_intersections[0][0])
                    if len(frozenset([p1,last_intersections[0][0]])) == 2:
                        new_edges.add(frozenset([p1,last_intersections[0][0]]))
                    if (last_intersections[0][0], True) in intersections:
                        intersections[intersections.index((last_intersections[0][0], True))] = last_intersections[0]
                    elif (last_intersections[0][0], False) in intersections:
                        #intersections[intersections.index((last_intersections[0][0], False))] = last_intersections[0]
                        pass
                    elif (last_intersections[0][0], None) in intersections:
                        intersections[intersections.index((last_intersections[0][0], None))] = last_intersections[0]
                    else:
                        intersections.extend(last_intersections)
            elif len(last_intersections) == 2:
                last_intersections = [y[1] for j,y in enumerate(sorted((distance(p1,x),x) for x in last_intersections))]
                if ((self.x_min > p1[0] or p1[0] > self.x_max) or (self.y_min > p1[1] or p1[1] > self.y_max) or (self.z_min > p1[2] or p1[2] > self.z_max)) and ((self.x_min > p2[0] or p2[0] > self.x_max) or (self.y_min > p2[1] or p2[1] > self.y_max) or (self.z_min > p2[2] or p2[2] > self.z_max)):
                    last_intersections[-1] = (last_intersections[-1], False)
                    last_intersections[-2] = (last_intersections[-2], True)
                    if len(frozenset([p1,last_intersections[0][0]])) == 2:
                        new_edges.add(frozenset([p1,last_intersections[0][0]]))
                    if len(frozenset([p2,last_intersections[-1][0]])) == 2:
                        new_edges.add(frozenset([p2,last_intersections[-1][0]]))
                elif (self.x_min > p1[0] or p1[0] > self.x_max) or (self.y_min > p1[1] or p1[1] > self.y_max) or (self.z_min > p1[2] or p1[2] > self.z_max):
                    last_intersections[-1] = (last_intersections[-1], None)
                    last_intersections[-2] = (last_intersections[-2], True)
                    if len(frozenset([p1,last_intersections[0][0]])) == 2:
                        new_edges.add(frozenset([p1,last_intersections[0][0]]))
                elif (self.x_min > p2[0] or p2[0] > self.x_max) or (self.y_min > p2[1] or p2[1] > self.y_max) or (self.z_min > p2[2] or p2[2] > self.z_max):
                    last_intersections[-1] = (last_intersections[-1], False)
                    last_intersections[-2] = (last_intersections[-2], None)
                    if len(frozenset([p2,last_intersections[-1][0]])) == 2:
                        new_edges.add(frozenset([p2,last_intersections[-1][0]]))
                else:
                    last_intersections[-1] = (last_intersections[-1], None)
                    last_intersections[-2] = (last_intersections[-2], None)
                    #if len(frozenset([last_intersections[-1][0],last_intersections[-2][0]])) == 2:
                    #    new_edges.add(frozenset([last_intersections[-1][0],last_intersections[-2][0]]))
                for intersection in last_intersections:
                    if (intersection[0], True) in intersections:
                        if intersection[1] is not None:
                            intersections[intersections.index((intersection[0], True))] = intersection
                    elif (intersection[0], False) in intersections:
                        if intersection[1] is not None:
                            intersections[intersections.index((intersection[0], False))] = intersection
                    elif (intersection[0], None) in intersections:
                        if intersection[1]:
                            intersections[intersections.index((intersection[0], None))] = intersection
                    else:
                        intersections.append(intersection)
        #print('WHO',new_edges)
        f = {frozenset(poly.verts[index] for index in poly.edges[edge_index]) for edge_index in face}
        for box_index, box_point in enumerate(box):
            edge1 = frozenset([box[box_index-1],box_point])
            for edge2 in f:
                for point in edge1:
                    if not Box.point_on_segment(round_edge(edge2), round_point(point)):
                        break
                else:
                    p1, p2 = edge2
                    edge3 = [tuple(float(y) for y in p1), sorted([(distance(p1,y),y) for y in edge1])[0][1]]
                    if distance(edge3[0],edge3[1]) > 0:
                        new_edges.add(frozenset(edge3))
                    edge3 = [tuple(float(y) for y in p2), sorted([(distance(p2,y),y) for y in edge1])[0][1]]
                    if distance(edge3[0],edge3[1]) > 0:
                        new_edges.add(frozenset(edge3))
                    break
        box_old = [x for x in box]
        print('intersections', intersections)
        if len(intersections) > 1:
            print(box)
            if circuit[0][0] == circuit[1][0] and circuit[0][0] == circuit[2][0]:
                if round_float(circuit[0][0]) == round_float(self.x_min):
                    box_i = 0
                if round_float(circuit[0][0]) == round_float(self.x_max):
                    box_i = 1
            elif circuit[0][1] == circuit[1][1] and circuit[0][1] == circuit[2][1]:
                if round_float(circuit[0][1]) == round_float(self.y_min):
                    box_i = 2
                if round_float(circuit[0][1]) == round_float(self.y_max):
                    box_i = 3
            elif circuit[0][2] == circuit[1][2] and circuit[0][2] == circuit[2][2]:
                if round_float(circuit[0][2]) == round_float(self.z_min):
                    box_i = 4
                if round_float(circuit[0][2]) == round_float(self.z_max):
                    box_i = 5
            for intersection in intersections:
                for j in range(len(box)):
                    #print(frozenset([box[j-1],box[j]]), intersection)
                    if Box.point_on_segment(round_edge(frozenset([box[j-1],box[j]])), round_point(intersection[0])) and round_point(box[j-1]) != round_point(intersection[0]) and round_point(box[j]) != round_point(intersection[0]):
                        box.insert(j,intersection[0])
                        break
            for i, intersection in enumerate(intersections):
                try:
                    index = box.index(intersection[0])
                    path1 = [intersection[0]]
                except ValueError:
                    index = box.index(round_point(intersection[0]))
                    path1 = [round_point(intersection[0])]
                inside1 = intersection[1]
                index += 1
                index %= len(box)
                while True:
                    path1.append(box[index])
                    if (round_point(box[index]),True) in [(round_point(x[0]),x[1])for x in intersections]:
                        inside2 = True
                        break
                    if (round_point(box[index]),False) in [(round_point(x[0]),x[1])for x in intersections]:
                        inside2 = False
                        break
                    if (round_point(box[index]),None) in [(round_point(x[0]),x[1])for x in intersections]:
                        inside2 = None
                        break
                    index += 1
                    index %= len(box)
                try:
                    index = box.index(intersection[0])
                    path2 = [intersection[0]]
                except ValueError:
                    index = box.index(round_point(intersection[0]))
                    path2 = [round_point(intersection[0])]
                inside3 = intersection[1]
                index -= 1
                index %= len(box)
                while True:
                    path2.append(box[index])
                    if (round_point(box[index]),True) in [(round_point(x[0]),x[1])for x in intersections]:
                        inside4 = True
                        break
                    if (round_point(box[index]),False) in [(round_point(x[0]),x[1])for x in intersections]:
                        inside4 = False
                        break
                    if (round_point(box[index]),None) in [(round_point(x[0]),x[1])for x in intersections]:
                        inside4 = None
                        break
                    index -= 1
                    index %= len(box)
                print("path1",path1,inside1,inside2)
                print("path2",path2,inside3,inside4)
                is_border_path1 = inside1 is not None and inside2 is not None and not inside1 and inside2
                is_border_path2 = inside3 is not None and inside4 is not None and not inside3 and inside4
                if path1[-1] == path2[-1] and is_border_path1 and is_border_path2:
                    distance1 = sum(distance(path1[i-1],x) for i,x in enumerate(path1[:-1]))
                    distance2 = sum(distance(path2[i-1],x) for i,x in enumerate(path2[:-1]))
                    if distance2 > distance1:
                        is_border_path1 = False
                    else:
                        is_border_path2 = False
                for path, is_border_path in [(path1, is_border_path1), (path2, is_border_path2)]:
                    if is_border_path:
                        for k,z in enumerate(path[:-1]):
                            point = tuple((z[j]+path[k+1][j])/2 for j in range(3))
                            if not Box.inside_polyhedron(poly,point):
                                break
                        else:
                            for j,x in enumerate(path[:-1]):
                                edge1 = frozenset([x,path[j+1]])
                                for edge2 in new_edges:
                                    for point in edge1:
                                        if not Box.point_on_segment(round_edge(edge2), round_point(point)):
                                            break
                                    else:
                                        p1, p2 = edge2
                                        edge3 = [p1, sorted([(distance(p1,y),y) for y in edge1])[0][1]]
                                        if distance(edge3[0],edge3[1]) > 0:
                                            new_edges.add(frozenset(edge3))
                                        edge3 = [p2, sorted([(distance(p2,y),y) for y in edge1])[0][1]]
                                        if distance(edge3[0],edge3[1]) > 0:
                                            new_edges.add(frozenset(edge3))
                                        break
                                else:
                                    if len(edge1) == 2:
                                        new_edges.add(edge1)
                                        path_edges.add(edge1)
            '''
            for i, intersection1 in enumerate(intersections):
                if intersection1[1] is None:
                    continue
                for j, intersection2 in enumerate(intersections[i+1:]+intersections):
                    if intersection2[1] is None:
                        continue
                    break
                path = (intersections[i:]+intersections)[:min(len(intersections),j+2)]
                #print('temp',path)
                for k, intersection3 in enumerate(path):
                    edge1 = frozenset([path[k-1][0], intersection3[0]])
                    edge1 = frozenset([tuple(float(round(y*1000)/1000) for y in x) for x in edge1])
                    p1, p2 = edge1
                    if not any(Box.point_on_segment(edge2, tuple((p1[l]+p2[l])/2 for l in range(3))) for edge2 in new_edges|{frozenset([poly.verts[index] for index in poly.edges[edge_index]]) for edge_index in face}):
                        continue
                    for edge2 in new_edges:
                        for point in edge1:
                            if not Box.point_on_segment(edge2, point):
                                break
                        else:
                            if edge2 in new_edges:
                                new_edges.remove(edge2)
                                p1, p2 = edge2
                                edge3 = [tuple(float(x) for x in p1), sorted([(distance(p1,x),x) for x in edge1])[0][1]]
                                if distance(edge3[0],edge3[1]) > 0:
                                    new_edges.add(frozenset(edge3))
                                edge3 = [tuple(float(x) for x in p2), sorted([(distance(p2,x),x) for x in edge1])[0][1]]
                                if distance(edge3[0],edge3[1]) > 0:
                                    new_edges.add(frozenset(edge3))
                            else:
                                new_edges.add(edge1)
                            break
                    else:
                        new_edges.add(edge1)
                    for edge2 in {frozenset([poly.verts[index] for index in poly.edges[edge_index]]) for edge_index in face}:
                        for point in edge1:
                            if not Box.point_on_segment(edge2, point):
                                break
                        else:
                            new_edges.add(edge1)
                            break
            '''
            #if all(any(Box.point_on_segment(frozenset([x[j-1],y]), box_point) for x in poly.circuits(face_index) for j,y in enumerate(x)) for box_point in box_old):
            if all((box_point, None) in intersections for box_point in box_old) and not (len(circuit) == 4 and not len(set(circuit)-set(box_old))):
                for k, intersection3 in enumerate(intersections):
                    edge1 = frozenset([intersections[k-1][0], intersection3[0]])
                    p1, p2 = edge1
                    if not any(Box.point_on_segment(round_edge(edge2), round_point(tuple((p1[l]+p2[l])/2 for l in range(3)))) for edge2 in new_edges|{frozenset([poly.verts[index] for index in poly.edges[edge_index]]) for edge_index in face}):
                        continue
                    for edge2 in new_edges|{frozenset([poly.verts[index] for index in poly.edges[edge_index]]) for edge_index in face}:
                        for point in edge1:
                            if not Box.point_on_segment(round_edge(edge2), round_point(point)):
                                break
                        else:
                            if edge2 in new_edges:
                                new_edges.remove(edge2)
                                p1, p2 = edge2
                                edge3 = [p1, sorted([(distance(p1,x),x) for x in edge1])[0][1]]
                                if distance(edge3[0],edge3[1]) > 0:
                                    new_edges.add(frozenset(edge3))
                                edge3 = [p2, sorted([(distance(p2,x),x) for x in edge1])[0][1]]
                                if distance(edge3[0],edge3[1]) > 0:
                                    new_edges.add(frozenset(edge3))
                            else:
                                new_edges.add(edge1)
                                path_edges.add(edge1)
                            break
                    else:
                        new_edges.add(edge1)
                        path_edges.add(edge1)
            '''
            else:
                for i in range(len(box)):
                    if not any(Polyhedron.inside_triangle(y,tuple((box[i-1][j]+box[i][j])/2 for j in range(3))) for y in Polyhedron.triangulate(circuit)):
                        break
                else:
                    for i, intersection in enumerate(intersections):
                        index = box.index(intersection[0])
                        edge = frozenset([tuple(float(round(x*1000)/1000) for x in box[index]),tuple(float(round(x*1000)/1000) for x in box[(index+1)%len(box)])])
                        if edge in new_edges:
                            new_edges.remove(edge)
                        else:
                            new_edges.add(edge)
                        edge = frozenset([tuple(float(round(x*1000)/1000) for x in box[index]),tuple(float(round(x*1000)/1000) for x in box[(index-1)%len(box)])])
                        if edge in new_edges:
                            new_edges.remove(edge)
                        else:
                            new_edges.add(edge)
            while True:
                combine = False
                for edge1 in new_edges:
                    for edge2 in new_edges:
                        edge1_int = frozenset(tuple(int(round(x*1000)) for x in point) for point in edge1)
                        edge2_int = frozenset(tuple(int(round(x*1000)) for x in point) for point in edge2)
                        #print(edge1, edge2)
                        if len(edge1_int&edge2_int)==1 and Box.colinear(edge1|edge2):
                            new_edges.add(frozenset(tuple(x/1000 for x in point) for point in edge1_int^edge2_int))
                            new_edges.remove(edge1)
                            new_edges.remove(edge2)
                            combine = True
                            break
                    if combine:
                        break
                else:
                    break
                '''
        elif len(intersections) == 0:
            #print('box',box)
            if len(box) and any([Polyhedron.inside_triangle(x,box[0]) for x in Polyhedron.triangulate(circuit)]):
                for j, box_point in enumerate(box):
                    new_edges.add(frozenset([box[j-1],box_point]))
                if circuit[0][0] == circuit[1][0] and circuit[0][0] == circuit[2][0]:
                    if round_float(circuit[0][0]) == round_float(self.x_min):
                        box_i = 0
                    if round_float(circuit[0][0]) == round_float(self.x_max):
                        box_i = 1
                elif circuit[0][1] == circuit[1][1] and circuit[0][1] == circuit[2][1]:
                    if round_float(circuit[0][1]) == round_float(self.y_min):
                        box_i = 2
                    if round_float(circuit[0][1]) == round_float(self.y_max):
                        box_i = 3
                elif circuit[0][2] == circuit[1][2] and circuit[0][2] == circuit[2][2]:
                    if round_float(circuit[0][2]) == round_float(self.z_min):
                        box_i = 4
                    if round_float(circuit[0][2]) == round_float(self.z_max):
                        box_i = 5
            for i,x in enumerate(circuit):
                if circuit[i-1] in poly.verts and x in poly.verts and frozenset([poly.verts.index(circuit[i-1]),poly.verts.index(x)]) in poly.edges:
                    new_edges.add(frozenset([circuit[i-1],x]))
        print('face processing time', poly.circuits(face_index), time.time()-start_time_face)
        return new_edges, path_edges, box_i
    def delete(self, poly):
        start_time = time.time()
        box_faces = []
        box_faces.extend([{(x,y,z) for y in (self.y_min,self.y_max) for z in (self.z_min,self.z_max)} for x in (self.x_min,self.x_max)])
        box_faces.extend([{(x,y,z) for x in (self.x_min,self.x_max) for z in (self.z_min,self.z_max)} for y in (self.y_min,self.y_max)])
        box_faces.extend([{(x,y,z) for x in (self.x_min,self.x_max) for y in (self.y_min,self.y_max)} for z in (self.z_min,self.z_max)])
        if not Box.inside_polyhedron(poly, ((self.x_min+self.x_max)/2,(self.y_min+self.y_max)/2,(self.z_min+self.z_max)/2)):
            for edge in poly.edges:
                edge = frozenset(poly.verts[index] for index in edge)
                if self.intersect(edge) and not any(all(Box.colinear(frozenset([(x1,y1,z1),(x2,y2,z2)])|frozenset([point])) for point in edge) for x1 in (self.x_min,self.x_max) for y1 in (self.y_min,self.y_max) for z1 in (self.z_min,self.z_max) for x2 in (self.x_min,self.x_max) for y2 in (self.y_min,self.y_max) for z2 in (self.z_min,self.z_max) if sum([x1 != x2,y1 != y2, z1 != z2]) == 1) and not any(Box.coplanar(edge|x) for x in box_faces):
                    print('intersect', edge)
                    triple_break = False
                    for p1 in [x for x in edge if self.x_min <= x[0] and x[0] <= self.x_max and self.y_min <= x[1] and x[1] <= self.y_max and self.z_min <= x[2] and x[2] <= self.z_max]:
                        p2 = list(edge-frozenset([p1]))[0]
                        vec = tuple(p2[i]-p1[i] for i in range(3))
                        for i in range(20):
                            if vec[0]/2**i == 0 and vec[1]/2**i == 0 and vec[2]/2**i == 0:
                                break
                            p3 = tuple(p1[i]+vec[i]/2**i for i in range(3))
                            if self.x_min <= p3[0] and p3[0] <= self.x_max and self.y_min <= p3[1] and p3[1] <= self.y_max and self.z_min <= p3[2] and p3[2] <= self.z_max:
                                triple_break = True
                                break
                        if triple_break:
                            break
                    if triple_break:
                        break
            else:
                return poly
        with Pool(8) as p:
            output = p.starmap(self.delete_preprocess, [(poly,i) for i in range(len(poly.faces))])
            faces = [x[0] for x in output]
            print('POOL')
            for face in faces:
                print(face)
            temp = [x[1] for x in output]
            box_map = [None for i in range(6)]
            print([x[2] for x in output])
            for j,y in enumerate([x[2] for x in output]):
                if y is not None:
                    box_map[y] = poly.faces[j]
            print('box_map', box_map)
        path_edges = set()
        for x in temp:
            path_edges |= x
        print('preprocessing time', time.time()-start_time)
        old_faces = [face for face_index,face in enumerate(poly.faces) if len(faces[face_index])]
        old_face_indices = [poly.faces.index(x) for x in old_faces]
        faces = [face for face in faces if len(face)]
        if all(x is None for x in box_map):
            box_map = [False for i in range(6)]
        for face_index,face in enumerate(box_map):
            if face_index == 0:
                box = [(self.x_min, self.y_min, self.z_min),(self.x_min, self.y_max, self.z_min),(self.x_min, self.y_max, self.z_max),(self.x_min, self.y_min, self.z_max)]
            if face_index == 1:
                box = [(self.x_max, self.y_min, self.z_min),(self.x_max, self.y_max, self.z_min),(self.x_max, self.y_max, self.z_max),(self.x_max, self.y_min, self.z_max)]
            if face_index == 2:
                box = [(self.x_min, self.y_min, self.z_min),(self.x_max, self.y_min, self.z_min),(self.x_max, self.y_min, self.z_max),(self.x_min, self.y_min, self.z_max)]
            if face_index == 3:
                box = [(self.x_min, self.y_max, self.z_min),(self.x_max, self.y_max, self.z_min),(self.x_max, self.y_max, self.z_max),(self.x_min, self.y_max, self.z_max)]
            if face_index == 4:
                box = [(self.x_min, self.y_min, self.z_min),(self.x_max, self.y_min, self.z_min),(self.x_max, self.y_max, self.z_min),(self.x_min, self.y_max, self.z_min)] 
            if face_index == 5:
                box = [(self.x_min, self.y_min, self.z_max),(self.x_max, self.y_min, self.z_max),(self.x_max, self.y_max, self.z_max),(self.x_min, self.y_max, self.z_max)]
            if face is None:
                if all(Box.inside_polyhedron(poly, x) for x in box):
                    print('inside', box)
                    faces.append(set())
                    for i,x in enumerate(box):
                        faces[-1].add(round_edge(frozenset([box[i-1],x])))
            elif face != False and any(not any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(Polyhedron.find_exterior_circuit(poly.circuits(poly.faces.index(face))))) for y in box) and all(Box.inside_polyhedron(poly,y) and not any(Polyhedron.inside_triangle(x,y) for face_index,face in enumerate(poly.faces) for x in Polyhedron.triangulate(Polyhedron.circuit_cut(poly.circuits(face_index)))) for y in box):# or any(Box.point_on_segment(frozenset(poly.verts[w] for w in poly.edges[z]), y) for z in face) for y in box):
                for i,x in enumerate(box):
                    edge1 = frozenset([box[i-1],x])
                    if face in old_faces:
                        for edge2 in faces[old_faces.index(face)]:
                            if Box.point_on_segment(edge2,box[i-1]) and Box.point_on_segment(edge2,x):
                                faces[old_faces.index(face)].remove(edge2)
                                p1, p2 = edge2
                                edge3 = [tuple(float(x) for x in p1), sorted([(distance(p1,x),x) for x in edge1])[0][1]]
                                if distance(edge3[0],edge3[1]) > 0:
                                    faces[old_faces.index(face)].add(frozenset(edge3))
                                edge3 = [tuple(float(x) for x in p2), sorted([(distance(p2,x),x) for x in edge1])[0][1]]
                                if distance(edge3[0],edge3[1]) > 0:
                                    faces[old_faces.index(face)].add(frozenset(edge3))
                                break
                        else:
                            faces[old_faces.index(face)].add(edge1)
            '''
            elif face != False:
                if face_index == 0 or face_index == 1:
                    edges = {frozenset(x) for x in zip([(self.x_min, self.y_min, self.z_min),(self.x_min, self.y_max, self.z_min),(self.x_min, self.y_max, self.z_max),(self.x_min, self.y_min, self.z_max)],[(self.x_max, self.y_min, self.z_min),(self.x_max, self.y_max, self.z_min),(self.x_max, self.y_max, self.z_max),(self.x_max, self.y_min, self.z_max)])}
                if face_index == 2 or face_index == 3:
                    edges = {frozenset(x) for x in zip([(self.x_min, self.y_min, self.z_min),(self.x_max, self.y_min, self.z_min),(self.x_max, self.y_min, self.z_max),(self.x_min, self.y_min, self.z_max)],[(self.x_min, self.y_max, self.z_min),(self.x_max, self.y_max, self.z_min),(self.x_max, self.y_max, self.z_max),(self.x_min, self.y_max, self.z_max)])}
                if face_index == 4 or face_index == 5:
                    edges = {frozenset(x) for x in zip([(self.x_min, self.y_min, self.z_min),(self.x_max, self.y_min, self.z_min),(self.x_max, self.y_max, self.z_min),(self.x_min, self.y_max, self.z_min)],[(self.x_min, self.y_min, self.z_max),(self.x_max, self.y_min, self.z_max),(self.x_max, self.y_max, self.z_max),(self.x_min, self.y_max, self.z_max)])}
                triple_break = False
                for p1 in box:
                    for edge in edges:
                        if p1 in edge:
                            p2 = list(edge-frozenset([p1]))[0]
                            vec = tuple(p2[i]-p1[i] for i in range(3))
                            for i in range(20):
                                if vec[0]/2**i == 0 and vec[1]/2**i == 0 and vec[2]/2**i == 0:
                                    break
                                if Box.inside_polyhedron(poly, tuple(p1[i]+vec[i]/2**i for i in range(3))):
                                    triple_break = True
                                    break
                            if triple_break:
                                break
                    else:
                        for edge in list(path_edges):
                            if edge in faces[old_faces.index(face)]:
                                faces[old_faces.index(face)].remove(edge)
                    if triple_break:
                        break
            '''
        circuits = Box.delete_circuit_helper(path_edges)
        #print('path_edges', path_edges)
        print('circuits', circuits)
        for circuit in circuits:
            for face_index, face in enumerate(faces):
                if face_index >= len(old_face_indices):
                    break
                if Box.coplanar(set(round_point(x) for x in circuit)|{round_point(point) for edge in face for point in edge}):
                    #print('combine', circuit, poly.circuits(old_face_indices[face_index]))
                    for i,x in enumerate(circuit):
                        exterior_circuit = Polyhedron.find_exterior_circuit(poly.circuits(old_face_indices[face_index]))
                        if frozenset([circuit[i-1],x]) in face:
                            if not any(Box.point_on_segment(frozenset([exterior_circuit[j-1],y]),circuit[i-1]) and Box.point_on_segment(frozenset([exterior_circuit[j-1],y]),x) for j,y in enumerate(exterior_circuit)):
                                #print('howdy', frozenset([circuit[i-1],x]))
                                face.remove(frozenset([circuit[i-1],x]))
                        else:
                            #if not any(Box.point_on_segment(frozenset([exterior_circuit[j-1],y]),circuit[i-1]) for j,y in enumerate(exterior_circuit)) or not any(Box.point_on_segment(frozenset([exterior_circuit[j-1],y]),x) for j,y in enumerate(exterior_circuit)):
                            face.add(frozenset([circuit[i-1],x]))
                    break
            else:
                faces.append({frozenset([circuit[i-1],x]) for i,x in enumerate(circuit)})
        for face in faces:
            print(face)
        while True:
            for face in faces:
                for edge1 in face:
                    combine = False
                    for edge2 in face:
                        if len(edge1&edge2)==1 and Box.colinear(edge1|edge2):
                            combine = True
                            break
                    if combine:
                        print('combine', edge1, edge2)
                        face.remove(edge1)
                        face.remove(edge2)
                        point_map = dict()
                        for point in edge1:
                            point_map[round_point(point)] = point
                        for point in edge2:
                            point_map[round_point(point)] = point
                        if len(round_edge(edge1)^round_edge(edge2)) == 2:
                            face.add(frozenset(point_map[point] for point in round_edge(edge1)^round_edge(edge2)))
                        else:
                            maxi = (0, None)
                            for x in edge1:
                                for y in edge2:
                                    if distance(x,y) > maxi[0]:
                                        maxi = (distance(x,y),frozenset([x,y]))
                            face.add(maxi[1])
                        #print("combine",edge1^edge2)
                        break
                if combine:
                    break
            else:
                break
        faces = [face for face in faces if len(face)]
        verts = list(set(point for face in faces for edge in face for point in edge))
        edges = list(set(frozenset(point for point in edge) for face in faces for edge in face))
        #print('BLAH',edges)
        faces = list(set(frozenset(edges.index(frozenset(point for point in edge)) for edge in face) for face in faces))
        edges = [frozenset(verts.index(point) for point in edge) for edge in edges]
        new_poly = Polyhedron()
        new_poly.verts = verts
        new_poly.edges = edges
        new_poly.faces = faces
        for face_index, face in enumerate(new_poly.faces):
            pass
            #print('face', [frozenset(verts[y] for y in edges[x]) for x in face])
            #new_poly.circuits(face_index)
        print("delete time", time.time()-start_time)
        return new_poly
        '''
        faces = [negative[face] for face in negative]
        verts = list(set(point for face in faces for edge in face for point in edge))
        edges = list(set(edge for face in faces for edge in face))
        faces = [frozenset(edges.index(edge) for edge in face) for face in faces]
        edges = [frozenset(verts.index(point) for point in edge) for edge in edges]
        anti_poly = Polyhedron()
        anti_poly.verts = verts
        anti_poly.edges = edges
        anti_poly.faces = faces
        return new_poly, anti_poly
        '''
    
    def add_circuit_helper(face, start=None, previous=None, current=None, path=None):
        #print(path)
        edge_lookup = dict()
        for edge in face:
            for point in edge:
                if point not in edge_lookup:
                    edge_lookup[round_point(point)] = set()
                edge_lookup[round_point(point)].add(edge)
        circuits = set()
        if start is None:
            #print(path_edges)
            #print('edge_lookup',edge_lookup)
            for edge in face:
                path = []
                start = edge
                point = list(start)[0]
                path.append(point)
                temp = list(edge_lookup[round_point(point)] - set([start]))
                #print(path, temp)
                for i in range(len(temp)):
                    current = temp[i]
                    circuits.update(Box.add_circuit_helper(face, start, start, current, path))
                del path[-1]
        else:
            #print(current, previous,start)
            if current == start:
                #print(path)
                return {tuple(path)}
            point = list(current - previous)[0]
            if point in path:
                return {}
            path.append(point)
            previous = current
            temp = list(edge_lookup[round_point(point)]-set([previous]))
            #print(path, temp)
            for i in range(len(temp)):
                if Box.coplanar(set(round_point(x) for x in path)|round_edge(temp[i])):
                    current = temp[i]
                    circuits.update(Box.add_circuit_helper(face, start, previous, current, path))
            del path[-1]
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            for j,y in enumerate(circuits_list):
                if len(y)>len(x) and not len(set(x)-set(y)) and y in circuits:
                    circuits.remove(y)
                if len(y)==len(x) and not len(set(x)-set(y)) and j > i and y in circuits:
                    circuits.remove(y)
                if x != y and not len(set(x)-set(y)) and j > i and y in circuits:
                    y_r = tuple(reversed(y))
                    if (y[y.index(x[0]):]+y[:y.index(x[0])])[:len(x)] == x or (y_r[y_r.index(x[0]):]+y_r[:y_r.index(x[0])])[:len(x)] == x:
                        circuits.remove(y)
        return circuits
    def add(self, poly):
        edges = [frozenset(poly.verts[x] for x in edge) for edge in poly.edges]
        faces = [frozenset(edges[x] for x in face) for face in poly.faces]
        new_edges = []
        #0 x_min
        #1 x_max
        #2 y_min
        #3 y_max
        #4 z_min
        #5 z_max
        new_faces = [set() for i in range(6)]
        new_edges.append(frozenset([(self.x_min,self.y_min,self.z_min),(self.x_max,self.y_min,self.z_min)]))
        new_faces[2].add(new_edges[-1])
        new_faces[4].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_min,self.y_max,self.z_min),(self.x_max,self.y_max,self.z_min)]))
        new_faces[3].add(new_edges[-1])
        new_faces[4].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_min,self.y_min,self.z_max),(self.x_max,self.y_min,self.z_max)]))
        new_faces[2].add(new_edges[-1])
        new_faces[5].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_min,self.y_max,self.z_max),(self.x_max,self.y_max,self.z_max)]))
        new_faces[3].add(new_edges[-1])
        new_faces[5].add(new_edges[-1])

        new_edges.append(frozenset([(self.x_min,self.y_min,self.z_min),(self.x_min,self.y_max,self.z_min)]))
        new_faces[0].add(new_edges[-1])
        new_faces[4].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_max,self.y_min,self.z_min),(self.x_max,self.y_max,self.z_min)]))
        new_faces[1].add(new_edges[-1])
        new_faces[4].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_min,self.y_min,self.z_max),(self.x_min,self.y_max,self.z_max)]))
        new_faces[0].add(new_edges[-1])
        new_faces[5].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_max,self.y_min,self.z_max),(self.x_max,self.y_max,self.z_max)]))
        new_faces[1].add(new_edges[-1])
        new_faces[5].add(new_edges[-1])

        new_edges.append(frozenset([(self.x_min,self.y_min,self.z_min),(self.x_min,self.y_min,self.z_max)]))
        new_faces[0].add(new_edges[-1])
        new_faces[2].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_max,self.y_min,self.z_min),(self.x_max,self.y_min,self.z_max)]))
        new_faces[1].add(new_edges[-1])
        new_faces[2].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_min,self.y_max,self.z_min),(self.x_min,self.y_max,self.z_max)]))
        new_faces[0].add(new_edges[-1])
        new_faces[3].add(new_edges[-1])
        new_edges.append(frozenset([(self.x_max,self.y_max,self.z_min),(self.x_max,self.y_max,self.z_max)]))
        new_faces[1].add(new_edges[-1])
        new_faces[3].add(new_edges[-1])

        #print('faces')
        #for face in faces:
        #    print(face)
        #print('new_edges', new_edges)
        mapping = dict()
        for face_index1, face1 in enumerate(faces):
            circuit = Polyhedron.circuit_cut(poly.circuits(face_index1))
            if any(self.intersect(edge) for edge in face1):
                for face_index2, face2 in enumerate(new_faces):
                    print(circuit, poly.circuits(face_index1), face1)
                    if Box.coplanar({round_point(point) for edge in face1 for point in edge}|{round_point(point) for edge in face2 for point in edge}) and any(Box.intersect_segments(frozenset([circuit[i-1],x]),y) is not None for i,x in enumerate(circuit) for y in face2 if any(circuit[i-1] in poly.verts and x in poly.verts and poly.verts.index(circuit[i-1]) in z and poly.verts.index(x) in z for z in poly.edges)):
                        temp = set(new_edges)-new_faces[face_index2//2]-new_faces[face_index2//2+1]
                        double_break = False
                        triple_break = False
                        print(face1)
                        for p1 in circuit:
                            if p1 in poly.verts:
                                for edge1 in edges:
                                    if p1 in edge1:
                                        p2 = list(edge1-frozenset([p1]))[0]
                                        if p2 not in {point for edge in face1 for point in edge}:
                                            for edge2 in temp:
                                                p3, p4 = edge2
                                                if (Box.point_on_segment(edge2, p1) and Box.point_on_segment(edge2, p2)) or (Box.point_on_segment(edge1, p3) and Box.point_on_segment(edge1, p4)):
                                                    triple_break = True
                                                    break
                                    if triple_break:
                                        break
                                for edge1 in edges:
                                    p2 = list(edge1-frozenset([p1]))[0]
                                    if p2 not in {point for edge in face1 for point in edge} and self.intersect(edge1):# and not any(Box.colinear(edge1|edge2) for edge2 in new_edges):
                                        double_break = True
                                        break
                            if double_break or triple_break:
                                break
                        if double_break or triple_break:
                            if face_index2 not in mapping:
                                mapping[face_index2] = set()
                            mapping[face_index2].add(face1)
            for face_index2, face2 in enumerate(new_faces):
                if Box.coplanar({round_point(point) for edge in face1 for point in edge}|{round_point(point) for edge in face2 for point in edge}) and (any(Box.intersect_segments(frozenset([circuit[i-1],x]),y) is not None for i,x in enumerate(circuit) for y in face2 if any(circuit[i-1] in poly.verts and x in poly.verts and poly.verts.index(circuit[i-1]) in z and poly.verts.index(x) in z for z in poly.edges)) or all(any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit)) for y in {point for edge in face2 for point in edge}) or all(x[0] >= self.x_min and x[0] <= self.x_max and x[1] >= self.y_min and x[1] <= self.y_max and x[2] >= self.z_min and x[2] <= self.z_max for x in circuit)):
                    if all(any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit)) for y in {point for edge in face2 for point in edge}) or all(x[0] >= self.x_min and x[0] <= self.x_max and x[1] >= self.y_min and x[1] <= self.y_max and x[2] >= self.z_min and x[2] <= self.z_max for x in circuit):
                        if face_index2 not in mapping:
                            mapping[face_index2] = set()
                        mapping[face_index2].add(face1)
        #if not len(mapping) and len(poly.verts):
        #    return poly
        print('mapping')
        for key in mapping:
            print(new_faces[key])
            print(key, mapping[key])
        for face_index in mapping:
            for face2 in mapping[face_index]:
                for edge1 in face2:
                    for edge2 in new_faces[face_index]:
                        for point in edge2:
                            if not Box.point_on_segment(round_edge(edge1),round_point(point)):
                                break
                        else:
                            p1, p2 = edge1
                            new_faces[face_index].remove(edge2)
                            edge3 = [p1, sorted([(distance(p1,x),x) for x in edge2])[0][1]]
                            if distance(edge3[0],edge3[1]) > 0:
                                new_faces[face_index].add(frozenset(edge3))
                            edge3 = [p2, sorted([(distance(p2,x),x) for x in edge2])[0][1]]
                            if distance(edge3[0],edge3[1]) > 0:
                                new_faces[face_index].add(frozenset(edge3))
                            break
                        for point in edge1:
                            if not Box.point_on_segment(round_edge(edge2),round_point(point)):
                                break
                        else:
                            p1, p2 = edge2
                            new_faces[face_index].remove(edge2)
                            edge3 = [p1, sorted([(distance(p1,x),x) for x in edge1])[0][1]]
                            if distance(edge3[0],edge3[1]) > 0:
                                new_faces[face_index].add(frozenset(edge3))
                            edge3 = [p2, sorted([(distance(p2,x),x) for x in edge1])[0][1]]
                            if distance(edge3[0],edge3[1]) > 0:
                                new_faces[face_index].add(frozenset(edge3))
                            break
                        if Box.colinear(edge1|edge2):
                            p1, p2 = edge1
                            p3, p4 = edge2
                            if Box.point_on_segment(edge2, p1) and Box.point_on_segment(edge1, p3):
                                new_faces[face_index].remove(edge2)
                                p5 = sorted([(distance(p2,x),x) for x in edge2])[0][1]
                                p6 = sorted([(distance(p4,x),x) for x in edge1])[0][1]
                                new_faces[face_index].add(frozenset([p2,p5]))
                                new_faces[face_index].add(frozenset([p4,p6]))
                                break
                            if Box.point_on_segment(edge2, p2) and  Box.point_on_segment(edge1, p3):
                                new_faces[face_index].remove(edge2)
                                p5 = sorted([(distance(p1,x),x) for x in edge2])[0][1]
                                p6 = sorted([(distance(p4,x),x) for x in edge1])[0][1]
                                new_faces[face_index].add(frozenset([p1,p5]))
                                new_faces[face_index].add(frozenset([p4,p6]))
                                break
                            if Box.point_on_segment(edge2, p1) and  Box.point_on_segment(edge1, p4):
                                new_faces[face_index].remove(edge2)
                                p5 = sorted([(distance(p2,x),x) for x in edge2])[0][1]
                                p6 = sorted([(distance(p3,x),x) for x in edge1])[0][1]
                                new_faces[face_index].add(frozenset([p2,p5]))
                                new_faces[face_index].add(frozenset([p3,p6]))
                                break
                            if Box.point_on_segment(edge2, p2) and  Box.point_on_segment(edge1, p4):
                                new_faces[face_index].remove(edge2)
                                p5 = sorted([(distance(p1,x),x) for x in edge2])[0][1]
                                p6 = sorted([(distance(p3,x),x) for x in edge1])[0][1]
                                new_faces[face_index].add(frozenset([p1,p5]))
                                new_faces[face_index].add(frozenset([p3,p6]))
                                break
                    else:
                        new_faces[face_index].add(edge1)
            faces = [face for face in faces if face not in mapping[face_index]]
        #print('len(new_faces)', [len(x) for x in new_faces]) 
        new_faces = [face for face in new_faces if len(face)]
        faces = [set(face) for face in faces]
        while True:
            for face in faces:
                for edge1 in face:
                    combine = False
                    for edge2 in face:
                        if len(edge1&edge2)==1 and Box.colinear(edge1|edge2):
                            combine = True
                            break
                    if combine:
                        #print('combine', edge1, edge2)
                        face.remove(edge1)
                        face.remove(edge2)
                        point_map = dict()
                        for point in edge1:
                            point_map[round_point(point)] = point
                        for point in edge2:
                            point_map[round_point(point)] = point
                        if len(round_edge(edge1)^round_edge(edge2)) == 2:
                            face.add(frozenset(point_map[point] for point in round_edge(edge1)^round_edge(edge2)))
                        else:
                            maxi = (0, None)
                            for x in edge1:
                                for y in edge2:
                                    if distance(x,y) > maxi[0]:
                                        maxi = (distance(x,y),frozenset([x,y]))
                            face.add(maxi[1])
                        #print("combine",edge1^edge2)
                        break
                if combine:
                    break
            else:
                break
        for face in list(new_faces):
            circuits = Box.add_circuit_helper(face) 
            if len(circuits) > 1 and Polyhedron.find_exterior_circuit(circuits) is None:
                circuits = list(circuits)
                for circuit in circuits[1:]:
                    new_faces.append(set())
                    for i,x in enumerate(circuit):
                        new_faces[-1].add(frozenset([circuit[i-1],x]))
                        face.remove(frozenset([circuit[i-1],x]))
        faces = faces + new_faces
        faces = [frozenset(face) for face in faces]
        verts = list(set(tuple(float(x) for x in point) for face in faces for edge in face for point in edge))
        edges = list(set(edge for face in faces for edge in face))
        faces = [frozenset([edges.index(edge) for edge in face]) for face in faces]
        edges = [frozenset([verts.index(tuple(float(x) for x in point)) for point in edge]) for edge in edges]
        new_poly = Polyhedron()
        new_poly.verts = verts
        new_poly.edges = edges
        new_poly.faces = faces
        while True:
            double_break = False
            for i,x in enumerate(new_poly.verts):
                for j,y in enumerate(new_poly.verts[i+1:]):
                    if distance(x,y) < 0.001:
                        del new_poly.verts[i+1+j]
                        for k,z in enumerate(new_poly.edges):
                            p1, p2 = z
                            if p1 == i+1+j:
                                p1 = i
                            elif p1 > i+1+j:
                                p1 -= 1
                            if p2 == i+1+j:
                                p2 = i
                            elif p2 > i+1+j:
                                p2 -= 1
                            new_poly.edges[k] = frozenset([p1,p2])
                        double_break = True
                        break
                if double_break:
                    break
            else:
                break
        while True:
            double_break = False
            for i,x in enumerate(new_poly.edges):
                for j,y in enumerate(new_poly.edges[i+1:]):
                    if {new_poly.verts[z] for z in x} == {new_poly.verts[z] for z in y}:
                        del new_poly.edges[i+1+j]
                        for k,z in enumerate(new_poly.faces):
                            temp = list(z)
                            for l,w in enumerate(temp):
                                if w > i+1+j:
                                    temp[l] -= 1
                                elif w == i+1+j:
                                    temp[l] = i
                            new_poly.faces[k] = frozenset(temp)
                        double_break = True
                        break
                if double_break:
                    break
            else:
                break
        #new_poly.verts = [round_point(x) for x in new_poly.verts]
        #for face_index, face in enumerate(new_poly.faces):
        #    print(face_index, face, set(frozenset(verts[y] for y in edges[x]) for x in face))
        #    print(new_poly.circuits(face_index))
        return new_poly
    '''
    def add(self, poly):
        new_edges = set()
        new_edges.add(frozenset([(self.x_min,self.y_min,self.z_min),(self.x_max,self.y_min,self.z_min)]))
        new_edges.add(frozenset([(self.x_min,self.y_max,self.z_min),(self.x_max,self.y_max,self.z_min)]))
        new_edges.add(frozenset([(self.x_min,self.y_min,self.z_max),(self.x_max,self.y_min,self.z_max)]))
        new_edges.add(frozenset([(self.x_min,self.y_max,self.z_max),(self.x_max,self.y_max,self.z_max)]))

        new_edges.add(frozenset([(self.x_min,self.y_min,self.z_min),(self.x_min,self.y_max,self.z_min)]))
        new_edges.add(frozenset([(self.x_max,self.y_min,self.z_min),(self.x_max,self.y_max,self.z_min)]))
        new_edges.add(frozenset([(self.x_min,self.y_min,self.z_max),(self.x_min,self.y_max,self.z_max)]))
        new_edges.add(frozenset([(self.x_max,self.y_min,self.z_max),(self.x_max,self.y_max,self.z_max)]))
    
        new_edges.add(frozenset([(self.x_min,self.y_min,self.z_min),(self.x_min,self.y_min,self.z_max)]))
        new_edges.add(frozenset([(self.x_max,self.y_min,self.z_min),(self.x_max,self.y_min,self.z_max)]))
        new_edges.add(frozenset([(self.x_min,self.y_max,self.z_min),(self.x_min,self.y_max,self.z_max)]))
        new_edges.add(frozenset([(self.x_max,self.y_max,self.z_min),(self.x_max,self.y_max,self.z_max)]))

        edges = [frozenset(poly.verts[x] for x in edge) for edge in poly.edges]
        print('OLD EDGES', edges)
        faces = [set(edges[x] for x in face) for face in poly.faces]
        print('OLD FACES')
        for face in faces:
            print(face)
        edges = set(edges)
        del_edges = {edge for edge in edges if edge in new_edges}
        print("DEL_EDGES", del_edges)
        new_new_edges = set()
        for edge1 in edges:
            for edge2 in new_edges:
                for point in edge1:
                    if not Box.point_on_segment(edge2,point):
                        break
                else:
                    print("del",edge2)
                    del_edges.add(edge2)
                    del_edges.add(edge1)
                    p1, p2 = edge2
                    edge3 = [tuple(round(x*1000)/1000 for x in p1), sorted([(distance(p1,x),tuple(round(y*1000)/1000 for y in x)) for x in edge1])[0][1]]
                    print('edge3', edge3)
                    if distance(edge3[0],edge3[1]) > 0:
                        new_new_edges.add(frozenset(edge3))
                        for face in faces:
                            if edge1 in face or edge2 in face:
                                face.add(frozenset(edge3))
                    edge3 = [tuple(round(x*1000)/1000 for x in p2), sorted([(distance(p2,x),tuple(round(y*1000)/1000 for y in x)) for x in edge1])[0][1]]
                    print('edge3', edge3)
                    if distance(edge3[0],edge3[1]) > 0:
                        new_new_edges.add(frozenset(edge3))
                        for face in faces:
                            if edge1 in face or edge2 in face:
                                face.add(frozenset(edge3))
                for point in edge2:
                    if not Box.point_on_segment(edge1,point):
                        break
                else:
                    print("del",edge1)
                    del_edges.add(edge2)
                    del_edges.add(edge1)
                    p1, p2 = edge1
                    edge3 = [tuple(round(x*1000)/1000 for x in p1), sorted([(distance(p1,x),tuple(round(y*1000)/1000 for y in x)) for x in edge1])[0][1]]
                    print('edge3', edge3)
                    if distance(edge3[0],edge3[1]) > 0:
                        new_new_edges.add(frozenset(edge3))
                        for face in faces:
                            if edge1 in face or edge2 in face:
                                face.add(frozenset(edge3))
                    edge3 = [tuple(round(x*1000)/1000 for x in p2), sorted([(distance(p2,x),tuple(round(y*1000)/1000 for y in x)) for x in edge1])[0][1]]
                    print('edge3', edge3)
                    if distance(edge3[0],edge3[1]) > 0:
                        new_new_edges.add(frozenset(edge3))
                        for face in faces:
                            if edge1 in face or edge2 in face:
                                face.add(frozenset(edge3))
        print("del edges",del_edges)
        print("new new",new_new_edges)
        for edge in del_edges:
            for face in faces:
                if edge in face:
                    face.remove(edge)
        edges -= del_edges
        new_edges -= del_edges
        new_edges |= new_new_edges
        combine = False
        while True:
            for edge1 in new_edges:
                combine = False
                for edge2 in edges:
                    edge1_int = frozenset(tuple(int(round(x*1000)) for x in point) for point in edge1)
                    edge2_int = frozenset(tuple(int(round(x*1000)) for x in point) for point in edge2)
                    print(edge1_int, edge2_int, len(edge1_int&edge2_int)==1, Box.colinear(edge1|edge2))
                    if len(edge1_int&edge2_int)==1 and Box.colinear(edge1|edge2):
                        combine = True
                        break
                if combine:
                    edges.remove(edge2)
                    for face in faces:
                        if edge2 in face:
                            face.remove(edge2)
                            face.add(frozenset(tuple(x/1000 for x in point) for point in edge1_int^edge2_int))
                    new_edges.remove(edge1)
                    print("combine",edge1^edge2)
                    new_edges.add(frozenset(tuple(x/1000 for x in point) for point in edge1_int^edge2_int))
                    break
            else:
                break
        for edge1 in new_edges:
            for face in faces:
                if len(face) > 2:
                    if Box.coplanar(set(point for edge in face for point in edge)|set(edge1)):
                        face.add(frozenset(tuple(round(x*1000)/1000 for x in point) for point in edge1))
        new_faces = [set() for i in range(6)]
        del_faces = set([frozenset()])
        faces = [frozenset(face) for face in faces]
        for face in faces:
            new_face_bool = [True for i in range(6)]
            for edge in face:
                edge = list(edge)
                if edge[0][0] != self.x_min or edge[1][0] != self.x_min:
                    new_face_bool[0] = False
                if edge[0][0] != self.x_max or edge[1][0] != self.x_max:
                    new_face_bool[1] = False
                if edge[0][1] != self.y_min or edge[1][1] != self.y_min:
                    new_face_bool[2] = False
                if edge[0][1] != self.y_max or edge[1][1] != self.y_max:
                    new_face_bool[3] = False
                if edge[0][2] != self.z_min or edge[1][2] != self.z_min:
                    new_face_bool[4] = False
                if edge[0][2] != self.z_max or edge[1][2] != self.z_max:
                    new_face_bool[5] = False

            for edge in face:
                if not len(self.intersect(edge)):
                    for i in range(6):
                        new_face_bool[i] = False
            
            if any(new_face_bool):
                new_faces[new_face_bool.index(True)].update({frozenset(tuple(round(z*1000)/1000 for z in y) for y in x) for x in face})
                del_faces.add(face)

        print('new_edges', new_edges)
        for edge in new_edges|edges:
            edge = list(edge)
            edge = [tuple(round(y*1000)/1000 for y in x) for x in edge]
            if edge[0][0] == self.x_min and edge[1][0] == self.x_min:
                new_faces[0].add(frozenset(edge))
            if edge[0][0] == self.x_max and edge[1][0] == self.x_max:
                new_faces[1].add(frozenset(edge))
            if edge[0][1] == self.y_min and edge[1][1] == self.y_min:
                new_faces[2].add(frozenset(edge))
            if edge[0][1] == self.y_max and edge[1][1] == self.y_max:
                new_faces[3].add(frozenset(edge))
            if edge[0][2] == self.z_min and edge[1][2] == self.z_min:
                new_faces[4].add(frozenset(edge))
            if edge[0][2] == self.z_max and edge[1][2] == self.z_max:
                new_faces[5].add(frozenset(edge))
        print('faces ', faces)

        for face1 in new_faces:
            if not len(face1):
                continue
            print('face1', face1)
            edge_lookup = dict()
            for edge in face1:
                for point in edge:
                    if tuple(int(round(x*1000)) for x in point) not in edge_lookup:
                        edge_lookup[tuple(int(round(x*1000)) for x in point)] = set()
                    edge_lookup[tuple(int(round(x*1000)) for x in point)].add(edge)
            for point in edge_lookup:
                print(point, edge_lookup[point])
            try:
                start = list(face1)[0]
                previous = start
                point = list(start)[0]
                circuit = []
                circuit.append(point)
                current = list(edge_lookup[tuple(int(round(x*1000)) for x in point)] - set([start]))[0]
                while current != start:
                    print("current",current)
                    point = list(current - previous)[0]
                    circuit.append(point)
                    previous = current
                    for p in edge_lookup:
                        if current in edge_lookup[p]:
                            edge_lookup[p].remove(current)
                    current = list(edge_lookup[tuple(int(round(x*1000)) for x in point)])[0]
            except IndexError:
                print('LOL', point, edge_lookup[tuple(int(round(x*1000)) for x in point)])
                continue
                circuit = []
            print('circuit', circuit)
            for face_index, face2 in enumerate(faces):
                if face2 in del_faces:
                    continue
                print('face2',face2)
                circuits = poly.circuits(face_index)
                exterior_circuit = Polyhedron.find_exterior_circuit(circuits)
                print('exterior circuit', exterior_circuit, circuit)
                if exterior_circuit is not None:
                    for c in circuits:
                        if c != exterior_circuit:
                            for i,x in enumerate(c):
                                face1.add(frozenset([c[i-1],x]))
                    circuit_in_circuit = True
                    if Box.coplanar(set(circuit)|set(exterior_circuit)):
                        for y in circuit:
                            if not any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(exterior_circuit)):
                                circuit_in_circuit = False
                                break
                        else:
                            temp = {x for x in face1}
                            for j,y in enumerate(circuit):
                                face1.remove(frozenset([circuit[j-1],y]))
                            exterior_circuit = list(exterior_circuit)
                            for j,y in enumerate(exterior_circuit):
                                for edge in temp:
                                    if Box.colinear(edge|{exterior_circuit[j-1],y}) and Box.point_on_segment(edge,exterior_circuit[j-1]):
                                        exterior_circuit[j-1] = sorted([(distance(y,x),x) for x in edge])[0][1]
                                    if Box.colinear(edge|{exterior_circuit[j-1],y}) and Box.point_on_segment(edge,y):
                                        exterior_circuit[j] = sorted([(distance(exterior_circuit[j-1],x),x) for x in edge])[0][1]
                            print('EXT CIRCUIT', exterior_circuit)
                            for j,y in enumerate(exterior_circuit):
                                face1.add(frozenset([exterior_circuit[j-1],y]))
                            pass
                        if not circuit_in_circuit:
                            for y in exterior_circuit:
                                if not any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit)):
                                    intersections = []
                                    for k,z in enumerate(exterior_circuit):
                                        p1 = exterior_circuit[k-1]
                                        p2 = z
                                        last_intersections = []
                                        for i,x in enumerate(circuit):
                                            if not Box.colinear(set([circuit[i-1],x,p1,p2])):
                                                for intersection in Box.intersect_segments(frozenset([circuit[i-1],x]),frozenset([p1,p2])):
                                                    print("intersection", intersection, p1, p2)
                                                    last_intersections.append(intersection)
                                        if len(last_intersections) == 0:
                                            new_edges.add(frozenset([p1,p2]))
                                        elif len(last_intersections) == 1:
                                            if (self.x_min > p2[0] or p2[0] > self.x_max) or (self.y_min > p2[1] or p2[1] > self.y_max) or (self.z_min > p2[2] or p2[2] > self.z_max):
                                                last_intersections[-1] = (last_intersections[-1],True)
                                                print("intersection p2", p2, last_intersections[0][0])
                                                if len(frozenset([p2,last_intersections[0][0]])) == 2:
                                                    new_edges.add(frozenset([p2,last_intersections[0][0]]))
                                                intersections.extend(last_intersections)
                                            elif (self.x_min > p1[0] or p1[0] > self.x_max) or (self.y_min > p1[1] or p1[1] > self.y_max) or (self.z_min > p1[2] or p1[2] > self.z_max):
                                                last_intersections[-1] = (last_intersections[-1],False)
                                                print("intersection p1", p1, last_intersections[0][0])
                                                if len(frozenset([p1,last_intersections[0][0]])) == 2:
                                                    new_edges.add(frozenset([p1,last_intersections[0][0]]))
                                                intersections.extend(last_intersections)
                                        elif len(last_intersections) == 2:
                                            last_intersections = [y[1] for j,y in enumerate(sorted((distance(p1,x),x) for x in last_intersections))]
                                            if (self.x_min > p1[0] or p1[0] > self.x_max) or (self.y_min > p1[1] or p1[1] > self.y_max) or (self.z_min > p1[2] or p1[2] > self.z_max):
                                                if len(frozenset([p1,last_intersections[0]])) == 2:
                                                    new_edges.add(frozenset([p1,last_intersections[0]]))
                                                intersections.append((last_intersections[0],False))
                                            if (self.x_min > p2[0] or p2[0] > self.x_max) or (self.y_min > p2[1] or p2[1] > self.y_max) or (self.z_min > p2[2] or p2[2] > self.z_max):
                                                if len(frozenset([p2,last_intersections[1][0]])) == 2:
                                                    new_edges.add(frozenset([p2,last_intersections[1]]))
                                                intersections.append((last_intersections[1],True))
                                    for intersection in intersections:
                                        for j in range(len(circuit)):
                                            if Box.point_on_segment(frozenset([circuit[j-1],circuit[j]]), intersection[0]):
                                                circuit.insert(j,intersection[0])
                                                break
                                    for i, intersection in enumerate(intersections):
                                        index = circuit.index(intersection[0])
                                        index += 1
                                        index %= len(circuit)
                                        path1 = [intersection[0]]
                                        while True:
                                            path1.append(circuit[index])
                                            if (circuit[index],True) in intersections or (circuit[index],False) in intersections:
                                                break
                                            index += 1
                                            index %= len(circuit)
                                        index = circuit.index(intersection[0])
                                        index -= 1
                                        index %= len(circuit)
                                        path2 = [intersection[0]]
                                        while True:
                                            path2.append(circuit[index])
                                            if (circuit[index],True) in intersections or (circuit[index],False) in intersections:
                                                break
                                            index -= 1
                                            index %= len(circuit)
                                        if not intersection[1]:
                                            if path1[-1] == intersections[(i+1)%len(intersections)][0] and path2[-1] == intersections[(i+1)%len(intersections)][0]:
                                                distance1 = sum(distance(y,path1[j+1]) for j,y in enumerate(path1[:-1]))
                                                distance2 = sum(distance(y,path2[j+1]) for j,y in enumerate(path2[:-1]))
                                                if distance2 < distance1:
                                                    path1, path2 = path2, path1
                                            if path1[-1] == intersections[(i+1)%len(intersections)][0]:
                                                for j,x in enumerate(path1[:-1]):
                                                    new_edges.add(frozenset([x,path1[j+1]]))
                                            elif path2[-1] == intersections[(i+1)%len(intersections)][0]:
                                                for j,x in enumerate(path2[:-1]):
                                                    new_edges.add(frozenset([x,path2[j+1]]))
                                    break
                        del_faces.add(face2)
        for face in faces:
            print('HEY', face)
        faces = [set(face) for face in faces]
        temp = set()
        for face1 in new_faces:
            for face2 in faces:
                if len(face2) and Box.coplanar(set([point for point in edge for edge in face1])|set([point for point in edge for edge in face2])):
                    print('face1', face1)
                    face2.update(face1)
                    break
            else:
                temp.add(frozenset(face1))
        for face in faces:
            while True:
                remove = False
                for edge1 in face:
                    for edge2 in face:
                        if edge1 != edge2:
                            for point in edge2:
                                if not Box.point_on_segment(edge1, point):
                                    break
                            else:
                                remove = True
                                break
                    if remove:
                        break
                if remove:            
                    face.remove(edge2)
                else:
                    break
        faces.extend([face for face in temp])
        faces = faces + new_faces
        faces = [frozenset(face) for face in faces]
        faces = [face for face in faces if face not in del_faces and len(face)]
        verts = list(set(point for face in faces for edge in face for point in edge))
        edges = list(set(edge for face in faces for edge in face))
        for face in faces:
            print(face)
        for face in del_faces:
            print('del face', face)
        faces = [frozenset([edges.index(edge) for edge in face]) for face in faces]
        edges = [frozenset(verts.index(point) for point in edge) for edge in edges]
        new_poly = Polyhedron()
        new_poly.verts = verts
        new_poly.edges = edges
        new_poly.faces = faces
        for face_index, face in enumerate(new_poly.faces):
            print(face_index, face, set(frozenset(verts[y] for y in edges[x]) for x in face))
            print(new_poly.circuits(face_index))
        return new_poly
        '''
             
class Polyhedron:
    def __init__(self):
        self.verts = []
        self.edges = []
        self.faces = []
    def circuits(self, face_index, start=None, previous=None, current=None, path=None, old_circuits=None):
        #print(path, previous, current)
        face = self.faces[face_index]
        #print([self.edges[x] for x in face])
        if path is not None:
            for x in old_circuits:
                if not len(set(path)-set(x)):
                    return set()
        edge_lookup = dict()
        for edge_index in face:
            edge = self.edges[edge_index]
            for index in edge:
                if index not in edge_lookup:
                    edge_lookup[index] = set()
                edge_lookup[index].add(edge)
        circuits = set()
        if start is None:
            seen = set()
            for edge in face:
                if frozenset(self.verts[index] for index in self.edges[edge]) in seen:
                    continue
                path = []
                start = edge
                point = list(self.edges[start])[0]
                path.append(self.verts[point])
                temp = edge_lookup[point] - set([self.edges[start]])
                for y in temp:
                    current = self.edges.index(y)
                    output = self.circuits(face_index, start, start, current, path, circuits)
                    circuits.update(output)
                    for circuit in output:
                        for i,x in enumerate(circuit):
                            seen.add(frozenset([circuit[i-1],x]))
        else:
            if current == start:
                return {tuple(path)}
            point = list(self.edges[current] - self.edges[previous])[0]
            if self.verts[point] in path:
                return set()
            path.append(self.verts[point])
            previous = current
            temp = list(edge_lookup[point]-set([self.edges[previous]]))
            for y in temp:
                if Box.coplanar(set([round_point(x) for x in path])|{round_point(self.verts[x]) for x in y}):
                    current = self.edges.index(y)
                    circuits.update(self.circuits(face_index, start, previous, current, path, old_circuits|circuits))
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            for j,y in enumerate(circuits_list[i+1:]):
                if x != y and not len(set(x)-set(y)):
                    y_r = tuple(reversed(y))
                    if (y[y.index(x[0]):]+y[:y.index(x[0])])[:len(x)] == x or (y_r[y_r.index(x[0]):]+y_r[:y_r.index(x[0])])[:len(x)] == x:
                        if y in circuits:
                            circuits.remove(y)
        return circuits
    def inside_triangle(triangle, point):
        if point in triangle:
            return True
        triangle = list(triangle)
        alpha, beta = symbols("alpha beta")
        #print(triangle)
        exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
        exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
        eqs = [Eq(x,point[i]) for i,x in enumerate(exprs)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            try:
                alpha = solutions[0][alpha]
            except KeyError:
                alpha = 0
            try:
                beta = solutions[0][beta]
            except KeyError:
                beta = 0
            try:
                if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1:
                    return True
            except TypeError:
                return True
        return False
    def cross_product_triplet(a,b,c):
        v1 = tuple(b[i]-a[i] for i in range(3))
        v2 = tuple(c[i]-b[i] for i in range(3))
        return cross3D(v1,v2)
    # Find the convex angles of a circuit
    def convex_angles(circuit):
        v1 = tuple(circuit[1][i]-circuit[0][i] for i in range(3))
        v2 = tuple(circuit[2][i]-circuit[0][i] for i in range(3))
        normal = cross3D(v1,v2)
        cross_product_dot_normal = []
        for i in range(len(circuit)):
            cross_product_dot_normal.append(dot(Polyhedron.cross_product_triplet(circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]),normal))
        signs = [copysign(1,x) for x in cross_product_dot_normal]
        return [x==signs[i-1] or x == signs[(i+1)%len(signs)] for i,x in enumerate(signs)]


    # A circuit is a representation of a Jordan Polygon
    def clip_ear(circuit):
        is_convex = Polyhedron.convex_angles(circuit)
        #print(circuit, is_convex)
        for is_convex in [is_convex, [not x for x in is_convex]]:
            for i in range(len(circuit)):
                if is_convex[i]:
                    for j,y in enumerate(circuit):
                        if j != i and j != (i+1)%len(circuit) and j != (i-1)%len(circuit):
                            if Polyhedron.inside_triangle([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]],y) and y not in [circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]]:
                                #print([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]],y)
                                break
                    else:
                        return tuple([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]]), tuple([x for k,x in enumerate(circuit) if k != i])
    def triangulate(circuit):
        #print(circuit)
        output = []
        remainder = circuit
        while len(remainder) > 3:
            ear, remainder = Polyhedron.clip_ear(remainder)
            output.append(ear)
        output.append(remainder)
        return output
    def find_exterior_circuit(circuits):
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            triangulation = Polyhedron.triangulate(x)
            for j,y in enumerate(circuits_list):
                if j != i:
                    double_break = False
                    for point in y:
                        if not any(Polyhedron.inside_triangle(z,point) for z in triangulation):
                            double_break = True
                            break
                    if double_break:
                        break
            else:
                return x
        return None
    def circuit_intersect(segment, circuits):
        p1, p2 = segment
        output = []
        for circuit in circuits:
            for i,x in enumerate(circuit):
                p3 = circuit[i-1]
                p4 = x
                alpha, beta = symbols("alpha beta")
                eqs = [Eq(alpha*p2[i]+(1-alpha)*p1[i], beta*p3[i]+(1-beta)*p4[i]) for i in range(3)]
                solutions = solve(eqs, dict=True)
                if len(solutions):
                    try:
                        alpha = solutions[0][alpha]
                    except KeyError:
                        alpha = None
                    try:
                        beta = solutions[0][beta]
                    except KeyError:
                        beta = None
                    try:
                        if (alpha is None or (alpha >= 0 and alpha <= 1)) and (beta is None or (beta >= 0 and beta <= 1)):
                            output.append((alpha,tuple(alpha*p2[i]+(1-alpha)*p1[i] for i in range(3)),(p3,p4),circuit))
                    except TypeError:
                        pass
        return output
    # Converts a polygon with one or more holes into a single circuit by forming cuts from the holes to the exterior
    def circuit_cut(circuits):
        #print(circuits)
        exterior = Polyhedron.find_exterior_circuit(circuits)
        if exterior is None:
            return None
        output = [x for x in circuits if x != exterior]
        output.append(exterior)
        while len(output) > 1:
            interior = output[0]
            for x in interior:
                double_break = False
                for y in output[-1]:
                    segment = (x,y)
                    if not len(Polyhedron.circuit_intersect(segment, output[:-1])):
                        double_break = True
                        break
                if double_break:
                    break
            intersections = Polyhedron.circuit_intersect(segment, output)
            intersections.sort()
            for intersection in intersections:
                if intersection[3] == output[-1]:
                    first_exterior_intersection = intersection
                    break
            for intersection in intersections:
                if intersection[3] == output[-1]:
                    continue
                last_interior_intersection = intersection
            for i,x in enumerate(output[-1]):
                if x in first_exterior_intersection[2] and output[-1][i-1] in first_exterior_intersection[2]:
                    break
            for j,y in enumerate(last_interior_intersection[3]):
                if y in last_interior_intersection[2] and last_interior_intersection[3][j-1] in last_interior_intersection[2]:
                    break
            output[-1] = output[-1][:i] + tuple([first_exterior_intersection[1], last_interior_intersection[1]]) + tuple(reversed(last_interior_intersection[3][:j])) + tuple(reversed(last_interior_intersection[3][j:])) + tuple([last_interior_intersection[1], first_exterior_intersection[1]]) + output[-1][i:]
            output.remove(last_interior_intersection[3])
            output[-1] = list(output[-1])
            i = 1
            while i < len(output[-1])+1:
                if output[-1][i%len(output[-1])] == output[-1][i-1]:
                    del output[-1][i%len(output[-1])]
                else:
                    i += 1
            output[-1] = tuple(output[-1])
        return output[0]
    def is_inside(self, point):
        #print('is_inside', point)
        if any(Polyhedron.inside_triangle(x,point) for face_index,face in enumerate(self.faces) for x in Polyhedron.triangulate(Polyhedron.circuit_cut(self.circuits(face_index)))):
            return True
        vec = (random.random()*2-1,random.random()*2-1,random.random()*2-1)
        output = []
        for face_index,face in enumerate(self.faces):
            circuit = Polyhedron.circuit_cut(self.circuits(face_index))
            min_distance, mini = float("inf"), None
            for triangle in Polyhedron.triangulate(circuit):
                triangle = list(triangle)
                alpha, beta, gamma = symbols("alpha beta gamma")
                exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
                exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
                eqs = [Eq(expr,point[i]+vec[i]*gamma) for i,expr in enumerate(exprs)]
                #func = lambda x: [(x[0]*triangle[0][i]+(1-x[0])*triangle[1][i])*x[1]+(1-x[1])*triangle[2][i]-(point[i]+vec[i]*x[2]) for i in range(3)]
                #solution = fsolve(func, (0,0,0))
                #alpha, beta, gamma = solution
                #if not np.allclose(func(solution), [0.0,0.0,0.0]):
                #    continue
                solutions = solve(eqs, dict=True)
                if len(solutions):
                    alpha, beta, gamma = solutions[0][alpha], solutions[0][beta], solutions[0][gamma]
                #(beta*(alpha*triangle[0][i]+(1-alpha)*triangle[1][i])+(1-beta)*triangle[2][i]-position[i]+direction[i]*gamma)**2
                    if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1 and gamma > 0:
                        point = tuple(alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3))
                        point = tuple(beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(point))
                        if Polyhedron.inside_triangle(triangle,point):
                            output.append((gamma, point, face_index))
                        break
        #print('is_inside', output)
        return len(output)%2==1
    def round_verts(self, round_point):
        self.verts = [round_point(vert) for vert in self.verts]
        while True:
            double_break = False
            for i,x in enumerate(self.verts):
                for j,y in enumerate(self.verts[i+1:]):
                    if distance(x,y) < 0.001:
                        del self.verts[i+1+j]
                        for k,z in enumerate(self.edges):
                            p1, p2 = z
                            if p1 == i+1+j:
                                p1 = i
                            elif p1 > i+1+j:
                                p1 -= 1
                            if p2 == i+1+j:
                                p2 = i
                            elif p2 > i+1+j:
                                p2 -= 1
                            self.edges[k] = frozenset([p1,p2])
                        double_break = True
                        break
                if double_break:
                    break
            else:
                break
        while True:
            double_break = False
            for i,x in enumerate(self.edges):
                for j,y in enumerate(self.edges[i+1:]):
                    if {self.verts[z] for z in x} == {self.verts[z] for z in y}:
                        del self.edges[i+1+j]
                        for k,z in enumerate(self.faces):
                            temp = list(z)
                            for l,w in enumerate(temp):
                                if w > i+1+j:
                                    temp[l] -= 1
                                elif w == i+1+j:
                                    temp[l] = i
                            self.faces[k] = frozenset(temp)
                        double_break = True
                        break
                if double_break:
                    break
            else:
                break
        while True:
            double_break = False
            for i,x in enumerate(self.faces):
                for j,y in enumerate(self.faces[i+1:]):
                    if {frozenset(self.verts[w] for w in self.edges[z]) for z in x} == {frozenset(self.verts[w] for w in self.edges[z]) for z in y}:
                        del self.faces[i+1+j]
                        double_break = True
                        break
                if double_break:
                    break
            else:
                break
    def circuit_overlap(circuits1, circuits2):
        circuit1 = Polyhedron.circuit_cut(circuits1)
        circuit2 = Polyhedron.circuit_cut(circuits2)
        if any(Box.intersect_segments(frozenset([circuit1[i-1],x]),frozenset([circuit2[j-1],y])) for i,x in enumerate(circuit1) for j,y in enumerate(circuit2)):
            return True
        if any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit1) for y in circuit2):
            return True
        if any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit2) for y in circuit1):
            return True
        return False

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
    faces = [{0,1,2,3},{4,5,6,7},{8,4,9,0},{10,6,11,2},{1,10,5,9},{3,11,7,8}]
    polyhedron.faces = [frozenset(face) for face in faces]
    return polyhedron

class Block:
    def __init__(self, width, height, depth, unit=1):
        self.size = (width, height, depth)
        self.block = get_cube((width/2,height/2,depth/2),factors=(width,height,depth))
        self.poly = Polyhedron()
        #self.poly = get_cube((width/2,height/2,depth/2),factors=(width,height,depth))
        self.select = [0,0,0,0]
        self.unit = unit
        self.select_size = [unit,unit,unit]
        self.meters = {1}
        self.polygons = set()
    def flip(self):
        self.polygons = set()
        print()
        print("verts",self.poly.verts)
        for face_index,face in enumerate(self.poly.faces):
            print("face",{frozenset(self.poly.verts[index] for index in self.poly.edges[edge_index]) for edge_index in face})
            points = Polyhedron.circuit_cut(self.poly.circuits(face_index))
            print("face",points)
            print('face', face)
            print()
            self.polygons.add(tuple(points))
    def draw(self, pygame, screen):
        screen_width, screen_height = screen.get_size()
        for edge in self.block.edges:
            p1, p2 = tuple(self.block.verts[index] for index in edge)
            #print(p1,p2)
            p1, p2 = camera.project(p1), camera.project(p2)
            p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
            p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
            pygame.draw.line(screen, "white", p1, p2)
        #print('faces', self.poly.faces)
        #start_time = time.time()
        '''
        with Pool(8) as p:
            circuits = p.map(self.poly.circuits, [face_index for face_index,face in enumerate(self.poly.faces)])
            points = p.map(Polyhedron.circuit_cut, circuits)
            for i,x in enumerate(points):
                points[i] = p.map(camera.project, points[i])
            for face_index,face in enumerate(self.poly.faces):
                points[face_index] = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in points[face_index]]
                color = "white"
                pygame.draw.polygon(screen, color, points[face_index])
        '''
        for points in self.polygons:
            points = [camera.project(x) for x in points]
            points = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in points]
            color = "white"
            pygame.draw.polygon(screen, color, points)
        '''
        for face_index,face in enumerate(self.poly.faces):
            #print("face",face, self.poly.circuits(face_index), {frozenset(self.poly.verts[index] for index in self.poly.edges[edge_index]) for edge_index in face})
            points = Polyhedron.circuit_cut(self.poly.circuits(face_index))
            #points = Polyhedron.circuit_cut(c[face_index])
            #print(points)
            #print({frozenset(self.poly.verts[index] for index in self.poly.edges[edge_index]) for edge_index in face})
            #centroid = tuple(sum(points[i][j] for i,x in enumerate(points))/len(points) for j in range(3))
            #distance = sum((centroid[i]-camera.focal[i])**2 for i in range(3))
            points = [camera.project(x) for x in points]
            points = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in points]
            color = "white"
            pygame.draw.polygon(screen, color, points)
        '''
        #print('face drawing time', time.time()-start_time)
        min_select = tuple(min(block.select[i],block.select[i]+block.select_size[i]) for i in range(3))
        max_select = tuple(max(block.select[i],block.select[i]+block.select_size[i]) for i in range(3))
        meter = 1
        for x in self.meters:
            meter *= x
        for d1 in range(3):
            for mult in [1,-1]:
                d2, d3 = [i for i in range(3) if i != d1]
                points = [[[0 for i in range(3)] for j in range(abs(int(self.select_size[d3]/meter)))] for k in range(abs(int(self.select_size[d2]/meter)))]
                for i in range(len(points)):
                    for j in range(len(points[i])):
                        points[i][j][d2] = min_select[d2] + meter*(i+.5)
                        points[i][j][d3] = min_select[d3] + meter*(j+.5)
                        points[i][j][d1] = mult*float('inf')
                for face_index,face in enumerate(self.poly.faces):
                    circuit = Polyhedron.circuit_cut(self.poly.circuits(face_index))
                    if circuit[0][d1] == circuit[1][d1] and circuit[1][d1] == circuit[2][d1]:
                        for point_row in points:
                            for point in point_row:
                                projection = list(point)
                                projection[d1] = circuit[0][d1]
                                projection = tuple(projection)
                                if mult == 1:
                                    select = max_select
                                else:
                                    select = min_select
                                if mult*circuit[0][d1] < mult*point[d1] and mult*circuit[0][d1] >= mult*select[d1] and any(Polyhedron.inside_triangle(x,projection) for x in Polyhedron.triangulate(circuit)):
                                   point[d1] = circuit[0][d1]
                seen = set()
                components = []
                for i in range(len(points)):
                    for j in range(len(points[i])):
                        if (i,j) in seen:
                            continue
                        if points[i][j][d1] == float('inf') or points[i][j][d1] == -float('inf'):
                            continue
                        component = {(i,j)}
                        seen.add((i,j))
                        queue = [(i-1,j),(i,j-1),(i+1,j),(i,j+1)]
                        while len(queue):
                            k,l = queue.pop()
                            if (k,l) in seen:
                                continue
                            try:
                                if k >= 0 and l >= 0 and points[i][j][d1] == points[k][l][d1]:
                                    component.add((k,l))
                                    seen.add((k,l))
                                    queue.extend([(k-1,l),(k,l-1),(k+1,l),(k,l+1)])
                            except IndexError:
                                pass
                        components.append(component)
                for component in components:
                    point_set = set()
                    for i in range(len(points)):
                        for j in range(len(points[i])):
                            if (i,j) in component:
                                if (i-1,j) not in component or (i,j-1) not in component or (i-1,j-1) not in component:
                                    point = list(points[i][j])
                                    point[d2] = min_select[d2] + i*meter
                                    point[d3] = min_select[d3] + j*meter
                                    point = tuple(point)
                                    point_set.add(point)
                                if (i-1,j) not in component or (i,j+1) not in component or (i-1,j+1) not in component:
                                    point = list(points[i][j])
                                    point[d2] = min_select[d2] + i*meter
                                    point[d3] = min_select[d3] + (j+1)*meter
                                    point = tuple(point)
                                    point_set.add(point)
                                if (i+1,j) not in component or (i,j-1) not in component or (i+1,j-1) not in component:
                                    point = list(points[i][j])
                                    point[d2] = min_select[d2] + (i+1)*meter
                                    point[d3] = min_select[d3] + j*meter
                                    point = tuple(point)
                                    point_set.add(point)
                                if (i+1,j) not in component or (i,j+1) not in component or (i+1,j+1) not in component:
                                    point = list(points[i][j])
                                    point[d2] = min_select[d2] + (i+1)*meter
                                    point[d3] = min_select[d3] + (j+1)*meter
                                    point = tuple(point)
                                    point_set.add(point)
                    point = tuple(float('inf') for i in range(3))
                    for x in point_set:
                        if (x[d2],x[d3]) < (point[d2],point[d3]):
                            point = x
                    path = []
                    '''
                    while True:
                        triple_break = False
                        for x in point_set:
                            for y in point_set:
                                if x != y:
                                    for z in point_set:
                                        if x != z and y != z:
                                            if Box.colinear({x,y,z}):
                                                point_set.remove(sorted([(distance(x,y),z), (distance(x,z),y), (distance(y,z),x)])[-1][1])
                                                triple_break = True
                                                break
                                if triple_break:
                                    break
                            if triple_break:
                                break
                        else:
                            break
                    '''
                    point_list = list(point_set)
                    rotated_point_mapping = {x:x for x in point_list}
                    reverse_point_mapping = {rotated_point_mapping[key]:key for key in rotated_point_mapping}
                    angle = [0,0,0]
                    while not len(path) or point != path[0]:
                        path.append(point)
                        keys = [key for key in rotated_point_mapping]
                        if len(path) > 1:
                            angle[d1] += math.atan2(reverse_point_mapping[path[-1]][d3],reverse_point_mapping[path[-1]][d2])
                            rotated_point_mapping = {x:rotated_point_mapping[keys[i]] for i,x in enumerate(rotate([(y[0]-reverse_point_mapping[path[-1]][0],y[1]-reverse_point_mapping[path[-1]][1],y[2]-reverse_point_mapping[path[-1]][2]) for y in keys], angle))}
                        else:
                            rotated_point_mapping = {(x[0]-reverse_point_mapping[path[-1]][0],x[1]-reverse_point_mapping[path[-1]][1],x[2]-reverse_point_mapping[path[-1]][2]):rotated_point_mapping[x] for x in keys}
                        reverse_point_mapping = {rotated_point_mapping[key]:key for key in rotated_point_mapping}
                        gift_wrap = []
                        for key in rotated_point_mapping:
                            origin = (0,0,0)
                            if key == origin:
                                continue
                            gift_wrap.append((math.atan2(key[d3],key[d2]),distance(key,origin),rotated_point_mapping[key]))
                        maxi = (-float('inf'),-float('inf'),None)
                        for x in gift_wrap:
                            if x > maxi:
                                maxi = x
                        point = maxi[2]
                    i = 0
                    while i < len(path):
                        if Box.colinear({path[i-1],path[i],path[(i+1)%len(path)]}):
                            del path[i]
                        else:
                            i += 1
                    path = [camera.project(x) for x in path]
                    path = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in path]
                    color = "gray"
                    pygame.draw.polygon(screen, color, path)
        for edge in self.poly.edges:
            p1, p2 = tuple(self.poly.verts[index] for index in edge)
            #print(p1,p2)
            p1, p2 = camera.project(p1), camera.project(p2)
            p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
            p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
            pygame.draw.line(screen, "black", p1, p2)
        '''
        for d1 in range(3):
            d2, d3 = [i for i in range(3) if i != d1]
            points = [[0 for i in range(3)] for j in range(4)]
            points[0][d2] = min_select[d2]
            points[0][d3] = min_select[d3]
            points[0][d1] = float('inf')
            points[1][d2] = max_select[d2]
            points[1][d3] = min_select[d3]
            points[1][d1] = float('inf')
            points[2][d2] = max_select[d2]
            points[2][d3] = max_select[d3]
            points[2][d1] = float('inf')
            points[3][d2] = min_select[d2]
            points[3][d3] = max_select[d3]
            points[3][d1] = float('inf')
            circuits = [None for i in range(4)]
            for face_index,face in enumerate(self.poly.faces):
                circuit = Polyhedron.circuit_cut(self.poly.circuits(face_index))
                if circuit[0][d1] == circuit[1][d1] and circuit[1][d1] == circuit[2][d1]:
                    for index, point in enumerate(points):
                        projection = list(point)
                        projection[d1] = circuit[0][d1]
                        projection = tuple(projection)
                        if circuit[0][d1] < point[d1] and circuit[0][d1] >= max_select[d1] and any(Polyhedron.inside_triangle(x,projection) for x in Polyhedron.triangulate(circuit)):
                            point[d1] = circuit[0][d1]
                            circuits[index] = circuit
            for point_index,point in enumerate(points):
                if point_index > 0 and points[point_index-1][d1] == point[d1]:
                    continue
                if point[d1] == float('inf'):
                    continue
                shadow = []
                last_j = None
                last_circuit = None
                last_intersection = None
                shadow.append(tuple(point))
                for i in range(1,5):
                    if points[(point_index+i)%4][d1] == point[d1]:
                        if points[(point_index+i-1)%4][d1] < point[d1]:
                            projection = list(points[(point_index+i)%4])
                            projection[d1] = last_circuit[d1]
                            projection = tuple(projection)
                            temp = list(points[(point_index+i-1)%4])
                            temp[d1] = projection[d1]
                            temp = tuple(temp)
                            for k,z in enumerate(last_circuit):
                                intersection = Polyhedron.intersect_segments(frozenset([projection,temp]),frozenset([last_circuit[k-1],z]))
                                if intersection is not None:
                                    break
                            j = last_j
                            if last_circuit[j-1] == last_intersection or last_circuit[j][d2] >= min_select[d2] and last_circuit[j][d2] <= max_select[d2] and last_circuit[j][d3] >= min_select[d3] and last_circuit[j][d3] <= max_select[d3]:
                                while j != k:
                                    temp = list(last_circuit[j])
                                    temp[d1] = point[d1]
                                    temp = tuple(temp)
                                    shadow.append(temp)
                                    j += 1
                                    j %= len(last_circuit)
                            elif last_circuit[j] == last_intersection or last_circuit[j-1][d2] >= min_select[d2] and last_circuit[j-1][d2] <= max_select[d2] and last_circuit[j-1][d3] >= min_select[d3] and last_circuit[j-1][d3] <= max_select[d3]:
                                j -= 1
                                while j != k:
                                    temp = list(last_circuit[j])
                                    temp[d1] = point[d1]
                                    temp = tuple(temp)
                                    shadow.append(temp)
                                    j -= 1
                                    j %= len(last_circuit)
                            shadow.append(intersection)
                        if points[(point_index+i-1)%4][d1] > point[d1]:
                            temp = list(points[(point_index+i-1)%4])
                            temp[d1] = point[d1]
                            temp = tuple(temp)
                            for k,z in enumerate(circuits[point_index]):
                                intersection = Box.intersect_segments(frozenset([tuple(points[(point_index+i)%4]),temp]),frozenset([circuits[point_index][k-1],z]))
                                if intersection is not None:
                                    break
                            j = last_j
                            angle_sum1 = 0
                            for l,w in enumerate(circuits[point_index]):
                                vec1 = tuple(w[i]-circuits[point_index][l-1][i] for i in range(3))
                                vec2 = tuple(circuits[point_index][(l+1)%len(circuits[point_index])][i]-w[i] for i in range(3))
                                angle_sum1 += math.atan((vec2[d2]-vec1[d2])/(vec2[d3]-vec1[d3]))
                            angle_sum2 = 0
                            for l,w in enumerate(points):
                                vec1 = tuple(w[i]-points[l-1][i] for i in range(3))
                                vec2 = tuple(points[(l+1)%len(points)][i]-w[i] for i in range(3))
                                angle_sum2 += math.atan((vec2[d2]-vec1[d2])/(vec2[d3]-vec1[d3]))
                            vec3 = tuple(points[(point_index+i)%4][i]-points[(point_index+i-1)%4][i] for i in range(3))
                            if copysign(1,angle_sum1) == copysign(1,angle_sum2):
                                if circuits[point_index][j]!=last_intersection and (circuits[point_index][j][d2] >= min_select[d2] and circuits[point_index][j][d2] <= max_select[d2] and circuits[point_index][j][d3] >= min_select[d3] and circuits[point_index][j][d3] <= max_select[d3]):
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j += 1
                                        j %= len(circuits[point_index])
                                elif circuits[point_index][j-1]!=last_intersection and (circuits[point_index][j-1][d2] >= min_select[d2] and circuits[point_index][j-1][d2] <= max_select[d2] and circuits[point_index][j-1][d3] >= min_select[d3] and circuits[point_index][j-1][d3] <= max_select[d3]):
                                    j -= 1
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j -= 1
                                        j %= len(circuits[point_index])
                            else:
                                if circuits[point_index][j]!=last_intersection and (circuits[point_index][j][d2] >= min_select[d2] and circuits[point_index][j][d2] <= max_select[d2] and circuits[point_index][j][d3] >= min_select[d3] and circuits[point_index][j][d3] <= max_select[d3]):
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j -= 1
                                        j %= len(circuits[point_index])
                                elif circuits[point_index][j-1]!=last_intersection and (circuits[point_index][j-1][d2] >= min_select[d2] and circuits[point_index][j-1][d2] <= max_select[d2] and circuits[point_index][j-1][d3] >= min_select[d3] and circuits[point_index][j-1][d3] <= max_select[d3]):
                                    j += 1
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j += 1
                                        j %= len(circuits[point_index])
                            if intersection != shadow[-1]:
                                shadow.append(intersection)
                        shadow.append(tuple(points[(point_index+i)%4]))
                    if points[(point_index+i)%4][d1] > point[d1] and points[(point_index+i-1)%4][d1] == point[d1]:
                        projection = list(points[(point_index+i)%4])
                        projection[d1] = point[d1]
                        projection = tuple(projection)
                        for j,y in enumerate(circuits[point_index]):
                            intersection = Box.intersect_segments(frozenset([projection,shadow[-1]]),frozenset([circuits[point_index][j-1],y]))
                            if intersection is not None:
                                break
                        if intersection != shadow[-1]:
                            shadow.append(intersection)
                        last_j = j
                        last_intersection = intersection
                    if points[(point_index+i)%4][d1] < point[d1] and points[(point_index+i-1)%4][d1] == point[d1]:
                        end_point = list(points[(point_index+i)%4])
                        end_point[d1] = point[d1]
                        end_point = tuple(end_point)
                        for face_index,face in enumerate(self.poly.faces):
                            circuit = Polyhedron.circuit_cut(self.poly.circuits(face_index))
                            if circuit[0][d1] == circuit[1][d1] and circuit[1][d1] == circuit[2][d1] and circuit[0][d1] >= max_select[d1] and circuit[0][d1] < point[d1]:
                                projection = list(points[(point_index+i)%4])
                                projection[d1] = circuits[point_index][d1]
                                projection = tuple(projection)
                                temp = list(shadow[-1])
                                temp[d1] = projection[d1]
                                temp = tuple(temp)
                                if any(Polyhedron.inside_triangle(x,projection) for x in Polyhedron.triangulate(circuit)):
                                    for j,y in enumerate(circuits[point_index]):
                                        intersection = Box.intersect_segments(frozenset([projection,temp]),frozenset([circuits[point_index][j-1],y]))
                                        if intersection is not None:
                                            break
                                    temp = list(intersection)
                                    temp[d1] = point[d1]
                                    temp = tuple(temp)
                                    if distance(temp,shadow[-1]) < distance(end_point,shadow[-1]):
                                        end_point = temp
                                        last_circuit = circuit
                                        last_j = j
                                        last_intersection = intersection
                        shadow.append(end_point)
                        for j,y in enumerate(circuits[point_index]):
                            intersection = Box.intersect_segments(frozenset([tuple(point),shadow[-1]]),frozenset([circuits[point_index][j-1],y]))
                            if intersection is not None:
                                break
                        shadow.append(intersection)
                del shadow[-1]
                if shadow[-1] == shadow[0]:
                    del shadow[-1]
                if len(shadow) > 1:
                    print(points, shadow, shadow[-1]==shadow[-2])
                if len({x for x in circuits[point_index] if x[d2] >= min_select[d2] and x[d2] <= max_select[d2] and x[d3] >= min_select[d3] and x[d3] <= max_select[d3]}) > 2:
                    if len(shadow) > 2:
                        shadow = [camera.project(x) for x in shadow]
                        shadow = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in shadow]
                        color = "gray"
                        pygame.draw.polygon(screen, color, shadow)
            points = [[0 for i in range(3)] for j in range(4)]
            points[0][d2] = min_select[d2]
            points[0][d3] = min_select[d3]
            points[0][d1] = -float('inf')
            points[1][d2] = max_select[d2]
            points[1][d3] = min_select[d3]
            points[1][d1] = -float('inf')
            points[2][d2] = max_select[d2]
            points[2][d3] = max_select[d3]
            points[2][d1] = -float('inf')
            points[3][d2] = min_select[d2]
            points[3][d3] = max_select[d3]
            points[3][d1] = -float('inf')
            circuits = [None for i in range(4)]
            for face_index,face in enumerate(self.poly.faces):
                circuit = Polyhedron.circuit_cut(self.poly.circuits(face_index))
                if circuit[0][d1] == circuit[1][d1] and circuit[1][d1] == circuit[2][d1]:
                    for index, point in enumerate(points):
                        projection = list(point)
                        projection[d1] = circuit[0][d1]
                        projection = tuple(projection)
                        if circuit[0][d1] > point[d1] and circuit[0][d1] <= min_select[d1] and any(Polyhedron.inside_triangle(x,projection) for x in Polyhedron.triangulate(circuit)):
                               point[d1] = circuit[0][d1]
                               circuits[index] = circuit
            for point_index,point in enumerate(points):
                if point_index > 0 and points[point_index-1][d1] == point[d1]:
                    continue
                if point[d1] == -float('inf'):
                    continue
                shadow = []
                last_j = None
                last_circuit = None
                shadow.append(tuple(point))
                for i in range(1,5):
                    if points[(point_index+i)%4][d1] == point[d1]:
                        if points[(point_index+i-1)%4][d1] > point[d1]:
                            projection = list(points[(point_index+i)%4])
                            projection[d1] = last_circuit[d1]
                            projection = tuple(projection)
                            temp = list(points[(point_index+i-1)%4])
                            temp[d1] = projection[d1]
                            temp = tuple(temp)
                            for k,z in enumerate(last_circuit):
                                intersection = Polyhedron.intersect_segments(frozenset([projection,temp]),frozenset([last_circuit[k-1],z]))
                                if intersection is not None:
                                    break
                            j = last_j
                            if last_circuit[j][d2] >= min_select[d2] and last_circuit[j][d2] <= max_select[d2] and last_circuit[j][d3] >= min_select[d3] and last_circuit[j][d3] <= max_select[d3]:
                                while j != k:
                                    temp = list(last_circuit[j])
                                    temp[d1] = point[d1]
                                    shadow.append(temp)
                                    temp = tuple(temp)
                                    j += 1
                                    j %= len(last_circuit)
                            elif last_circuit[j-1][d2] >= min_select[d2] and last_circuit[j-1][d2] <= max_select[d2] and last_circuit[j-1][d3] >= min_select[d3] and last_circuit[j-1][d3] <= max_select[d3]:
                                j -= 1
                                while j != k:
                                    temp = list(last_circuit[j])
                                    temp[d1] = point[d1]
                                    shadow.append(temp)
                                    temp = tuple(temp)
                                    j -= 1
                                    j %= len(last_circuit)
                            shadow.append(intersection)
                        if points[(point_index+i-1)%4][d1] < point[d1]:
                            temp = list(points[(point_index+i-1)%4])
                            temp[d1] = point[d1]
                            temp = tuple(temp)
                            for k,z in enumerate(circuits[point_index]):
                                intersection = Box.intersect_segments(frozenset([tuple(points[(point_index+i)%4]),temp]),frozenset([circuits[point_index][k-1],z]))
                                if intersection is not None:
                                    break
                            j = last_j
                            angle_sum1 = 0
                            for l,w in enumerate(circuits[point_index]):
                                vec1 = tuple(w[i]-circuits[point_index][l-1][i] for i in range(3))
                                vec2 = tuple(circuits[point_index][(l+1)%len(circuits[point_index])][i]-w[i] for i in range(3))
                                angle_sum1 += math.atan((vec2[d2]-vec1[d2])/(vec2[d3]-vec1[d3]))
                            angle_sum2 = 0
                            for l,w in enumerate(points):
                                vec1 = tuple(w[i]-points[l-1][i] for i in range(3))
                                vec2 = tuple(points[(l+1)%len(points)][i]-w[i] for i in range(3))
                                angle_sum2 += math.atan((vec2[d2]-vec1[d2])/(vec2[d3]-vec1[d3]))
                            if copysign(1,angle_sum1) == copysign(1,angle_sum2):
                                if circuits[point_index][j]!=last_intersection and (circuits[point_index][j][d2] >= min_select[d2] and circuits[point_index][j][d2] <= max_select[d2] and circuits[point_index][j][d3] >= min_select[d3] and circuits[point_index][j][d3] <= max_select[d3]):
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j += 1
                                        j %= len(circuits[point_index])
                                elif circuits[point_index][j-1]!=last_intersection and (circuits[point_index][j-1][d2] >= min_select[d2] and circuits[point_index][j-1][d2] <= max_select[d2] and circuits[point_index][j-1][d3] >= min_select[d3] and circuits[point_index][j-1][d3] <= max_select[d3]):
                                    j -= 1
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j -= 1
                                        j %= len(circuits[point_index])
                            else:
                                if circuits[point_index][j]!=last_intersection and (circuits[point_index][j][d2] >= min_select[d2] and circuits[point_index][j][d2] <= max_select[d2] and circuits[point_index][j][d3] >= min_select[d3] and circuits[point_index][j][d3] <= max_select[d3]):
                                    j -= 1
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j -= 1
                                        j %= len(circuits[point_index])
                                elif circuits[point_index][j-1]!=last_intersection and (circuits[point_index][j-1][d2] >= min_select[d2] and circuits[point_index][j-1][d2] <= max_select[d2] and circuits[point_index][j-1][d3] >= min_select[d3] and circuits[point_index][j-1][d3] <= max_select[d3]):
                                    while j != k:
                                        temp = list(circuits[point_index][j])
                                        temp[d1] = point[d1]
                                        temp = tuple(temp)
                                        shadow.append(temp)
                                        j += 1
                                        j %= len(circuits[point_index])
                            if intersection != shadow[-1]:
                                shadow.append(intersection)
                        shadow.append(tuple(points[(point_index+i)%4]))
                    if points[(point_index+i)%4][d1] < point[d1] and points[(point_index+i-1)%4][d1] == point[d1]:
                        projection = list(points[(point_index+i)%4])
                        projection[d1] = point[d1]
                        projection = tuple(projection)
                        for j,y in enumerate(circuits[point_index]):
                            intersection = Box.intersect_segments(frozenset([projection,shadow[-1]]),frozenset([circuits[point_index][j-1],y]))
                            if intersection is not None:
                                break
                        shadow.append(intersection)
                        last_j = j
                        last_intersection = intersection
                    if points[(point_index+i)%4][d1] > point[d1] and points[(point_index+i-1)%4][d1] == point[d1]:
                        end_point = list(points[(point_index+i)%4])
                        end_point[d1] = point[d1]
                        end_point = tuple(end_point)
                        for face_index,face in enumerate(self.poly.faces):
                            circuit = Polyhedron.circuit_cut(self.poly.circuits(face_index))
                            if circuit[0][d1] == circuit[1][d1] and circuit[1][d1] == circuit[2][d1] and circuit[0][d1] <= min_select[d1] and circuit[0][d1] > point[d1]:
                                projection = list(points[(point_index+i)%4])
                                projection[d1] = circuits[point_index][d1]
                                projection = tuple(projection)
                                temp = list(shadow[-1])
                                temp[d1] = projection[d1]
                                temp = tuple(temp)
                                if any(Polyhedron.inside_triangle(x,projection) for x in Polyhedron.triangulate(circuit)):
                                    for j,y in enumerate(circuits[point_index]):
                                        intersection = Box.intersect_segments(frozenset([projection,temp]),frozenset([circuits[point_index][j-1],y]))
                                        if intersection is not None:
                                            break
                                    temp = list(intersection)
                                    temp[d1] = point[d1]
                                    temp = tuple(temp)
                                    if distance(temp,shadow[-1]) < distance(end_point,shadow[-1]):
                                        end_point = temp
                                        last_circuit = circuit
                                        last_j = j
                        shadow.append(end_point)
                        for j,y in enumerate(circuits[point_index]):
                            intersection = Box.intersect_segments(frozenset([tuple(point),shadow[-1]]),frozenset([circuits[point_index][j-1],y]))
                            if intersection is not None:
                                break
                        shadow.append(intersection)
                del shadow[-1]
                if shadow[-1] == shadow[0]:
                    del shadow[-1]
                if len(shadow) > 1:
                    print(points, shadow, shadow[-1]==shadow[-2])
                if len({x for x in circuits[point_index] if x[d2] >= min_select[d2] and x[d2] <= max_select[d2] and x[d3] >= min_select[d3] and x[d3] <= max_select[d3] and x[d1] == point[d1]}) > 2:
                    if len(shadow) > 2:
                        shadow = [camera.project(x) for x in shadow]
                        shadow = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in shadow]
                        color = "gray"
                        pygame.draw.polygon(screen, color, shadow)
        '''
        width, height, depth = self.size
        if self.select[3] == 0:
            axes = get_cube((width/2,abs(self.select_size[1])/2+min(self.select[1],self.select[1]+self.select_size[1]),depth/2), factors=(width,abs(self.select_size[1]),depth))
        elif self.select[3] == 1:
            axes = get_cube((width/2,height/2,abs(self.select_size[2])/2+min(self.select[2],self.select[2]+self.select_size[2])), factors=(width,height,abs(self.select_size[2])))
        elif self.select[3] == 2:
            axes = get_cube((abs(self.select_size[0])/2+min(self.select[0],self.select[0]+self.select_size[0]),height/2,depth/2), factors=(abs(self.select_size[0]),height,depth))
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
        if self.select[2] > 0:
            shadow =  get_cube((self.select[0]+self.select_size[0]/2,self.select[1]+self.select_size[1]/2,0+self.select[2]/2), factors=(self.select_size[0],self.select_size[1],self.select[2]))
            for edge in shadow.edges:
                p1, p2 = tuple(shadow.verts[index] for index in edge)
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

    def reunit(self, unit):
        block = Block(self.size[0], self.size[1], self.size[2], unit)
        block.select = self.select
        block.poly = self.poly
        block.meters = self.meters
        block.polygons = self.polygons
        return block
    def select_by_void(self):
        if not len(self.poly.verts):
            return True
        points = [(x,y,z) for x in (self.select[0],self.select[0]+self.unit) for y in (self.select[1],self.select[1]+self.unit) for z in (self.select[2],self.select[2]+self.unit)]
        for point in points:
            if not Box.inside_polyhedron(self.poly, point):
                return True
        for vert in self.poly.verts:
            if self.select[0] <= vert[0] and vert[0] <= self.select[0]+self.unit and self.select[1] <= vert[1] and vert[1] <= self.select[1]+self.unit and self.select[2] <= vert[2] and vert[2] <= self.select[2]+self.unit:
                return True
        edges = {frozenset([x,y]) for x in points for y in points if sum(z!=y[k] for k,z in enumerate(x))==1}
        for edge1 in edges:
            for edge2 in [frozenset(self.poly.verts[index] for index in edge) for edge in self.poly.edges]:
                if all(Box.point_on_segment(edge1, point) for point in edge2) or all(Box.point_on_segment(edge2, point) for point in edge1):
                    return True
        faces = []
        faces.extend([[(x,y,z) for y in (self.select[1],self.select[1]+self.unit) for z in (self.select[2],self.select[2]+self.unit)] for x in (self.select[0],self.select[0]+self.unit)])
        faces[-1][-1], faces[-1][-2] = faces[-1][-2], faces[-1][-1]
        faces[-2][-1], faces[-2][-2] = faces[-2][-2], faces[-2][-1]
        faces.extend([[(x,y,z) for x in (self.select[0],self.select[0]+self.unit) for z in (self.select[2],self.select[2]+self.unit)] for y in (self.select[1],self.select[1]+self.unit)])
        faces[-1][-1], faces[-1][-2] = faces[-1][-2], faces[-1][-1]
        faces[-2][-1], faces[-2][-2] = faces[-2][-2], faces[-2][-1]
        faces.extend([[(x,y,z) for x in (self.select[0],self.select[0]+self.unit) for y in (self.select[1],self.select[1]+self.unit)] for z in (self.select[2],self.select[2]+self.unit)])
        faces[-1][-1], faces[-1][-2] = faces[-1][-2], faces[-1][-1]
        faces[-2][-1], faces[-2][-2] = faces[-2][-2], faces[-2][-1]
        for face1 in faces:
            for face_index, face2 in enumerate(self.poly.faces):
                if Box.coplanar(set(face1)|{self.poly.verts[index] for edge_index in face2 for index in self.poly.edges[edge_index]}) and Polyhedron.circuit_overlap((face1,),self.poly.circuits(face_index)):
                    return True
        return False


points = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
faces = [{0,1,2,3},{4,5,6,7},{8,4,9,0},{10,6,11,2},{1,10,5,9},{3,11,7,8}]
'''
cube = Polyhedron()
cube.verts = points
cube.edges = [frozenset(x) for x in edges]
cube.faces = [frozenset(x) for x in faces]
box = Box((0.25,0.75),(-1,2),(-1,2))

cube = box.delete(cube)
print([face for face in cube.faces])
#box = Box((0.25,0.75),(-1,2),(-1,2))
box = Box((0.25,0.75),(0,1),(0,1))
#box = Box((0.25,0.75),(-0.5,1),(-0.5,1))
#box = Box((0.25,0.75),(0,0.25),(0,0.25))
cube = box.add(cube)
'''
if __name__ == "__main__":
    import pygame
    camera = Camera()
    camera.zoom = 100
    camera.focal = (0,0,-100)
    pygame.init()
    screen = pygame.display.set_mode((1280, 720))
    screen_width, screen_height = screen.get_size()
    clock = pygame.time.Clock()
    pygame.mouse.set_visible(False)
    running = True
    mouse_down = False
    space_down = False
    meta_down = False
    delete = False
    dir_mult1 = 1
    dir_mult2 = 1
    dir_mult3 = 1
    z_forward = True
    block = Block(3,3,3)
    meters = {1}
    while running:
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
                divisor_keys = [pygame.K_2,pygame.K_3,pygame.K_4,pygame.K_5,pygame.K_6,pygame.K_7,pygame.K_8,pygame.K_9,pygame.K_0]
                if event.key in divisor_keys:
                    if not meta_down:
                        block = block.reunit(block.unit/(divisor_keys.index(event.key)+2))
                        meters = remeter(meters, block.unit)
                        block.meters = meters
                    else:
                        block = block.reunit(block.unit*(divisor_keys.index(event.key)+2))
                        block.select[0] = block.select[0]//block.unit*block.unit
                        block.select[1] = block.select[1]//block.unit*block.unit
                        block.select[2] = block.select[2]//block.unit*block.unit
                if event.key == pygame.K_1:
                    block = block.reunit(1)
                if event.key == pygame.K_BACKSPACE:
                    delete = not delete
                    #filename = input("Input filename: ")
                    #pickle.dump(block, open(filename,'wb'))
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
                            if block.select_size[index] == 0:
                                block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                                if block.select_size[index] == block.unit:
                                    block.select[index] -= block.unit
                                    block.select_size[index] += block.unit
                                else:
                                    block.select[index] += block.unit
                                    block.select_size[index] -= block.unit
                        elif block.select[3] == 1:
                            block.select_size[1] += direction*block.unit
                            if block.select_size[1] == 0:
                                block.select_size[1] = direction*block.unit
                                if block.select_size[1] == block.unit:
                                    block.select[1] -= block.unit
                                    block.select_size[1] += block.unit
                                else:
                                    block.select[1] += block.unit
                                    block.select_size[1] -= block.unit
                        elif block.select[3] == 2:
                            block.select_size[1] += direction*block.unit
                            if block.select_size[1] == 0:
                                block.select_size[1] = direction*block.unit
                                if block.select_size[1] == block.unit:
                                    block.select[1] -= block.unit
                                    block.select_size[1] += block.unit
                                else:
                                    block.select[1] += block.unit
                                    block.select_size[1] -= block.unit
                        for i in range(3):
                            if block.select_size[i] == 0:
                                block.select_size[i] = direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                                if block.select_size[i] == block.unit:
                                    block.select[i] -= block.unit
                                else:
                                    block.select[i] += block.unit
                                    block.select_size[i] -= block.unit
                            if block.select[i]+block.select_size[i] > block.size[i]:
                                block.select_size[i] = block.size[i]-block.select[i]
                                if block.select_size[i] == 0:
                                    block.select_size[i] -= block.unit
                            if block.select[i]+block.select_size[i] < 0:
                                block.select_size[i] = -block.select[i]
                                if block.select_size[i] == 0:
                                    block.select_size[i] += block.unit
                            if block.select[i] < 0:
                                block.select[i] += block.unit
                            if block.select[i] > block.size[i]:
                                block.select[i] -= block.unit
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
                            block.select[2] -= direction*dir_mult2*block.unit
                            while not block.select_by_void():
                                block.select[2] -= direction*dir_mult2*block.unit
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
                                if block.select_size[index] == block.unit:
                                    block.select[index] -= block.unit
                                    block.select_size[index] += block.unit
                                else:
                                    block.select[index] += block.unit
                                    block.select_size[index] -= block.unit
                        elif block.select[3] == 1:
                            block.select_size[0] += direction*dir_mult1*block.unit
                            if block.select_size[0] == 0:
                                block.select_size[0] = -direction*dir_mult1*block.unit
                                if block.select_size[0] == block.unit:
                                    block.select[0] -= block.unit
                                    block.select_size[0] += block.unit
                                else:
                                    block.select[0] += block.unit
                                    block.select_size[0] -= block.unit
                        elif block.select[3] == 2:
                            block.select_size[2] -= direction*dir_mult2*block.unit
                            if block.select_size[2] == 0:
                                block.select_size[2] = -direction*dir_mult2*block.unit
                                if block.select_size[2] == block.unit:
                                    block.select[2] -= block.unit
                                    block.select_size[2] += block.unit
                                else:
                                    block.select[2] += block.unit
                                    block.select_size[2] -= block.unit
                        for i in [0,2]:
                            if block.select[i]+block.select_size[i] > block.size[i]:
                                block.select_size[i] = block.size[i]-block.select[i]
                            if block.select[i]+block.select_size[i] < 0:
                                block.select_size[i] = -block.select[i]
                                if block.select_size[i] == 0:
                                    block.select_size[i] += block.unit
                            if block.select[i] < 0:
                                block.select[i] += block.unit
                            if block.select[i] > block.size[i]:
                                block.select[i] -= block.unit
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
                if event.key in {pygame.K_LMETA, pygame.K_RMETA}:
                    meta_down = True
            if event.type == pygame.KEYUP:
                if event.key == pygame.K_SPACE:
                    space_down = False
                    i,i_max = block.select[0],block.select[0]+block.select_size[0]
                    if block.select_size[0] < 0:
                        i,i_max = i_max,i
                    j,j_max = block.select[1],block.select[1]+block.select_size[1]
                    if block.select_size[1] < 0:
                        j,j_max = j_max,j
                    k,k_max = block.select[2],block.select[2]+block.select_size[2]
                    if block.select_size[2] < 0:
                        k,k_max = k_max,k
                    i = min(block.select[0],block.select[0]+block.select_size[0])
                    j = min(block.select[1],block.select[1]+block.select_size[1])
                    k = min(block.select[2],block.select[2]+block.select_size[2])
                    size = (abs(block.select_size[0]),abs(block.select_size[1]),abs(block.select_size[2]))
                    box = Box((float(i),float(i+size[0])),(float(j),float(j+size[1])),(float(k),float(k+size[2])))
                    #print('box',(i,i+size[0]),(j,j+size[1]),(k,k+size[2]))
                    block.poly = box.delete(block.poly)
                    if not delete:
                        block.poly = box.add(block.poly)
                    block.poly.round_verts(lambda x: round_point_meter(x, meters))
                    block.flip()
                    print(meters)
                    block.select = [block.select[0]+block.select_size[0]-(block.unit if block.select_size[0] > 0 else 0),block.select[1]+block.select_size[1]-(block.unit if block.select_size[1] > 0 else 0),block.select[2]+block.select_size[2]-(block.unit if block.select_size[2] > 0 else 0),block.select[3]]
                    block.select_size = [block.unit,block.unit,block.unit]
                if event.key in {pygame.K_LMETA, pygame.K_RMETA}:
                    meta_down = False
            if event.type == pygame.MOUSEMOTION:
                x,y = pygame.mouse.get_pos()
                x -= screen_width/2
                y = -(y-screen_height/2)
                #crosshair.pos = (x,y)
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
        '''
        for edge in cube.edges:
            p1, p2 = tuple(cube.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "white", p1, p2)
        '''
        block.draw(pygame,screen)
        pygame.display.flip()
