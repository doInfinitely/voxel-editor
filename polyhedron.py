from utility import *
from camera import Camera
from math import sin, cos, pi, sqrt, copysign, atan2
import numpy as np
import math
from sympy import nsolve, solve, symbols, Eq, diff
import sympy
import errno
import os
import signal
import functools
from scipy.spatial import ConvexHull
from scipy.spatial._qhull import QhullError
import subprocess
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
import sys
from itertools import permutations
import networkx as nx
from multiprocessing import Pool

def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator

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
            for start in face:
                if frozenset(self.verts[index] for index in self.edges[start]) in seen:
                    continue
                path = []
                point = list(self.edges[start])[0]
                temp = edge_lookup[point] - set([self.edges[start]])
                for y in temp:
                    current = self.edges.index(y)
                    output = self.circuits(face_index, start, start, current, path + [self.verts[point]], circuits)
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
            previous = current
            temp = list(edge_lookup[point]-set([self.edges[previous]]))
            for y in temp:
                current = self.edges.index(y)
                circuits.update(self.circuits(face_index, start, previous, current, path + [self.verts[point]], old_circuits|circuits))
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            for j,y in enumerate(circuits_list[i+1:]):
                if not len(set(x)-set(y)):
                    y_r = tuple(reversed(y))
                    if (y[y.index(x[0]):]+y[:y.index(x[0])])[:len(x)] == x or (y_r[y_r.index(x[0]):]+y_r[:y_r.index(x[0])])[:len(x)] == x:
                        if y in circuits:
                            circuits.remove(y)
        return circuits
    def colinear(points, rtol=0.0001):
        points = list(points)
        index = 1
        while index < len(points) and distance(points[0],points[index]) < rtol:
            index += 1
        if index >= len(points):
            return True
        if len(points) < 3:
            return True
        for p in points[1:index]+points[index+1:]:
            a = np.array([[float(points[index][i]-points[0][i])] for i in range(3)])
            b = np.array([float(p[i]-points[0][i]) for i in range(3)])
            x, res, rank, s = np.linalg.lstsq(a, b)
            if not distance(np.dot(a, x), b) <= rtol:
                return False
        return True
    def coplanar(points, rtol=0.0001):
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
        if Polyhedron.colinear(points, rtol):
            return True
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
            x, res, rank, s = np.linalg.lstsq(a.astype('float'), b.astype('float'))
            #print(a, b, x, res)
            if not distance(np.dot(a.astype('float'), x.astype('float')), b.astype('float')) <= rtol:
                return False
        return True
    def point_on_segment(edge, point, rtol=0.0001):
        p1, p2 = edge
        '''
        alpha = symbols("alpha")
        eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], point[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            alpha = solutions[0][alpha]
            if alpha >= 0 and alpha <= 1:
                return True
        return False
        '''
        a = np.array([[p1[i], p2[i]] for i in range(3)]+[[1,1]])
        b = np.array(list(point)+[1])
        x, res, rank, s = np.linalg.lstsq(a.astype('float'), b.astype('float'))
        #print(distance(point, tuple(x[0]*p1[i]+x[1]*p2[i] for i in range(3))), x[0]+x[1])
        if round_float(x[0]) < 0 or round_float(x[1]) < 0 or not np.allclose(x[0]+x[1], 1, rtol=rtol) or distance(point, tuple(x[0]*p1[i]+x[1]*p2[i] for i in range(3))) > 0.0001:
            return False
        return True
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
                    return tuple(float(alpha*p1[i]+(1-alpha)*p2[i]) for i in range(3))
            except TypeError:
                    pass
                    #alpha = 0
                    #output.append(tuple(alpha*p1[i]+(1-alpha)*p2[i] for i in range(3)))
        return None

    @timeout(1)
    def timed_solve(*args, **kwargs):
        return solve(*args, **kwargs)
    def inside_triangle(triangle, point):
        if point in triangle:
            return True
        triangle = list(triangle)
        #print(triangle)
        '''
        alpha, beta = symbols("alpha beta")
        #print(triangle)
        exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
        exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
        eqs = [Eq(x,point[i]) for i,x in enumerate(exprs)]
        try:
            solutions = Polyhedron.timed_solve(eqs, dict=True)
        except TimeoutError:
            return False
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
        '''
        a = np.array([[triangle[j][i] for j in range(3)] for i in range(3)] + [[1 for j in range(3)]]).astype('float')
        b = np.array([point[i] for i in range(3)] + [1]).astype('float')
        x = np.linalg.lstsq(a,b)[0]
        dist = distance(point, tuple(sum(x[j]*triangle[j][i] for j in range(3)) for i in range(3)))
        #print(dist, x)
        if dist <= 0.00001 and all(round_float(x[i]) >= 0 for i in range(3)):
            return True
        return False
    def intersect_triangle(segment, triangle):
        p1, p2 = segment
        triangle = list(triangle)
        '''
        alpha, beta, gamma = symbols("alpha beta gamma")
        exprs = [beta*triangle[0][i]+(1-beta)*triangle[1][i] for i in range(3)]
        exprs = [gamma*x+(1-gamma)*triangle[2][i] for i,x in enumerate(exprs)]
        eqs = [Eq(x,(1-alpha)*p1[i]+alpha*p2[i]) for i,x in enumerate(exprs)]
        try:
            solutions = Polyhedron.timed_solve(eqs, dict=True)
        except TimeoutError:
            return None
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
                gamma = solutions[0][gamma]
            except KeyError:
                gamma = None
            try:
                if (alpha is None or (alpha >= 0 and alpha <= 1)) and (beta is None or (beta >= 0 and beta <= 1)) and (beta is None or (beta >= 0 and beta <= 1)) and (gamma is None or (gamma >= 0 and gamma <= 1)):
                    #print('triangle intersect', tuple(alpha*p2[i]+(1-alpha)*p1[i] for i in range(3)))
                    return (alpha,tuple(alpha*p2[i]+(1-alpha)*p1[i] for i in range(3)))
            except TypeError:
                pass
        '''
        a = np.array([[triangle[j][i] for j in range(3)]+[p1[i]-p2[i]] for i in range(3)] + [[1 for j in range(3)]+[0]]).astype('float')
        b = np.array([p1[i] for i in range(3)] + [1]).astype('float')
        x = np.linalg.lstsq(a,b)[0]
        point = tuple(float(x[3]*p2[i]+(1-x[3])*p1[i]) for i in range(3))
        point_prime = tuple(float(sum(x[j]*triangle[j][i] for j in range(3))) for i in range(3))
        if round_float(x[0]) >= 0 and round_float(x[1]) >= 0 and round_float(x[2]) and round_float(x[3]) >= 0 and round_float(x[3]) <= 1 and distance(np.matmul(a,x)[:3], p1) <= 0.00001:
            if Polyhedron.inside_triangle(triangle, point):
                return (float(x[3]), point)
        return None
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
    def clockwise_angles(planar):
        output = []
        #temp = []
        for i in range(len(planar)):
            planar = [tuple(y[j]-planar[i][j] for j in range(3)) for y in planar]
            angle = (0,0,-math.atan2(planar[i][1]-planar[i-1][1],planar[i][0]-planar[i-1][0]))
            planar = rotate(planar, angle)
            theta = math.atan2(planar[(i+1)%len(planar)][1],planar[(i+1)%len(planar)][0])
            output.append(theta < 0)
            #print(planar, theta/math.pi*180)
            #temp.append(math.acos(dot(planar[i-1],planar[(i+1)%len(planar)])/distance((0,0,0), planar[i-1])/distance((0,0,0),planar[(i+1)%len(planar)]))/math.pi*180)
        #print(temp)
        return output
    def make_planar(circuit):
        vec1 = tuple(circuit[1][i]-circuit[0][i] for i in range(3))
        vec2 = tuple(circuit[2][i]-circuit[0][i] for i in range(3))
        vec0 = cross3D(vec1,vec2)
        vec1 = (0,vec0[1],vec0[2])
        vec2 = (vec0[0],0,vec0[2])
        angles = [0,0,0]
        try:
            angles[0] = -math.acos(max(min(1,vec1[2]/distance(vec1)),-1))
        except ZeroDivisionError:
            pass
        try:
            angles[1] = -math.acos(max(min(1,vec2[2]/distance(vec2)),-1))
        except ZeroDivisionError:
            pass
        #print(angles)
        return tuple((x[0],x[1],0) for x in rotate(circuit, angles))
        '''
        p1 = (distance(circuit[0],circuit[1]),0,0)
        angle = [0,0, math.asin((circuit[1][1]-circuit[0][1])/p1[0])]
        p1 = rotate([p1], angle)[0]
        #print(p1,tuple(circuit[1][i]-circuit[0][i] for i in range(3)))
        angle[1] = math.asin(min(max((circuit[1][2]-circuit[0][2])/p1[0],1),-1))
        p1 = rotate([p1], (0,angle[1],0))[0]
        #print(p1,tuple(circuit[1][i]-circuit[0][i] for i in range(3)))
        angle = [-x for x in angle]
        p2 = rotate([tuple(circuit[2][i]-circuit[0][i] for i in range(3))], angle)[0]
        #p2 = rotate([tuple(circuit[2][i]-circuit[0][i] for i in range(3))], (0,0,angle[2]))[0]
        #p2 = rotate([p2], (0,angle[1],0))[0]
        #print(p2, rotate([tuple(circuit[2][i]-circuit[0][i] for i in range(3))], angle)[0])
        angle[0] = -math.atan2(p2[2]-p1[2],p2[1]-p1[1])
        #print(angle)
        #p2 = rotate([p2], (angle[0],0,0))[0]
        #print(p2, tuple(circuit[2][i]-circuit[0][i] for i in range(3)))
        #print(rotate([tuple(x[i]-circuit[0][i] for i in range(3)) for x in circuit], angle))
        output = rotate([tuple(x[i]-circuit[0][i] for i in range(3)) for x in circuit], (0,angle[1],angle[2]))
        output = rotate(output, (angle[0],0,0))
        #print(temp)
        return tuple((y[0],y[1],y[2]) for y in output)
        '''
    def is_clockwise(planar):
        summa = 0
        #print(planar) 
        for i in range(len(planar)):
            planar = [tuple(y[j]-planar[i][j] for j in range(3)) for y in planar]
            angle = (0,0,-math.atan2(planar[i][1]-planar[i-1][1],planar[i][0]-planar[i-1][0]))
            planar = rotate(planar, angle)
            theta = math.atan2(planar[(i+1)%len(planar)][1],planar[(i+1)%len(planar)][0])
            summa += theta
        #print(summa)
        return summa < 0
    def make_clockwise(circuits):
        return tuple(x if Polyhedron.is_clockwise(Polyhedron.make_planar(x)) else tuple(reversed(x)) for x in circuits)
        
    # A circuit is a representation of a Jordan Polygon
    def clip_ear(circuit):
        planar = Polyhedron.make_planar(circuit)
        #print(planar)
        if not Polyhedron.is_clockwise(planar):
            planar = tuple(reversed(planar))
            circuit = tuple(reversed(circuit))
        #print("clockwise_angles")
        is_convex = Polyhedron.clockwise_angles(planar)
        #print(circuit, is_convex)
        #print(tuple(reversed(circuit)), Polyhedron.clockwise_angles(tuple(reversed(circuit))))
        for i in range(len(circuit)):
            if is_convex[i]:
                for j,y in enumerate(circuit):
                    if j != i and j != (i+1)%len(circuit) and j != (i-1)%len(circuit):
                        if Polyhedron.inside_triangle([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]],y) and round_point(y) not in [round_point(z) for z in [circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]]]:
                            #print([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]],y)
                            break
                else:
                    remainder = [x for k,x in enumerate(circuit) if k != i]
                    remainder = [x for k,x in enumerate(remainder) if x != remainder[k-1]]
                    return tuple([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]]), tuple(remainder)
        print('hey')
        1/0
        sys.exit(1)
    def triangulate(circuit):
        #print('triangulate', circuit)
        output = []
        remainder = circuit
        while len(remainder) > 3:
            ear, remainder = Polyhedron.clip_ear(remainder)
            #print("ear",ear)
            output.append(ear)
        output.append(remainder)
        #print("ears", output)
        return output
    def find_exterior_circuit(circuits):
        #print(circuits)
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            triangulation = Polyhedron.triangulate(x)
            for j,y in enumerate(circuits_list):
                if j != i:
                    double_break = False
                    for point in y:
                        #print(triangulation)
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
                            output.append((alpha,tuple(float(alpha*p2[i]+(1-alpha)*p1[i]) for i in range(3)),(p3,p4),circuit))
                    except TypeError:
                        pass
        return output
    def face_intersect(self, segment):
        output = []
        for face_index, face in enumerate(self.faces):
            for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(Polyhedron.make_clockwise(self.circuits(face_index)))):
                intersect = Polyhedron.intersect_triangle(segment, triangle)
                if intersect:
                    #print('inter', segment, triangle, intersect)
                    output.append(intersect+tuple([face_index]))
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
                    if not len(Polyhedron.circuit_intersect(segment, output[1:-1])):
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
            #print("circuit_cut")
            #print(output[-1][:i])
            #print(tuple([first_exterior_intersection[1], last_interior_intersection[1]]))
            #print(tuple(reversed(last_interior_intersection[3][:j])))
            #print(tuple(reversed(last_interior_intersection[3][j:])))
            #print(tuple([last_interior_intersection[1], first_exterior_intersection[1]]))
            #print(output[-1][i:])
            output[-1] = output[-1][:i] + tuple([first_exterior_intersection[1], last_interior_intersection[1]]) + tuple(reversed(last_interior_intersection[3][j:] + last_interior_intersection[3][:j])) + tuple([last_interior_intersection[1], first_exterior_intersection[1]]) + output[-1][i:]
            output.remove(last_interior_intersection[3])
            output[-1] = list(output[-1])
            i = 1
            while i < len(output[-1])+1:
                if round_point(output[-1][i%len(output[-1])]) == round_point(output[-1][i-1]):
                    del output[-1][i%len(output[-1])]
                else:
                    i += 1
            output[-1] = tuple(output[-1])
            #print(output[0])
            #print("end circuit_cut")
        point_map = dict()
        for point in output[0]:
            if round_point(point) in point_map:
                if distance(point, round_point(point)) < distance(point_map[round_point(point)], round_point(point)):
                    point_map[round_point(point)] = point
            else:
                point_map[round_point(point)] = point
        return tuple(point_map[round_point(x)] for x in output[0])
    def in_faces(self, point):
        output = []
        for face_index,face in enumerate(self.faces):
            for x in Polyhedron.triangulate(Polyhedron.circuit_cut(Polyhedron.make_clockwise(self.circuits(face_index)))):
                if Polyhedron.inside_triangle(x,point):
                    output.append(face_index)
                    break
        return output
    def is_inside(self, point, in_faces_check=True):
        #print('is_inside', point)
        if in_faces_check and len(self.in_faces(point)):
            return True
        vec = (random.random()*2-1,random.random()*2-1,random.random()*2-1)
        output = []
        for face_index,face in enumerate(self.faces):
            circuit = Polyhedron.circuit_cut(Polyhedron.make_clockwise(self.circuits(face_index)))
            for triangle in Polyhedron.triangulate(circuit):
                triangle = list(triangle)
                '''
                alpha, beta, gamma = symbols("alpha beta gamma")
                exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
                exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
                eqs = [Eq(expr,point[i]+vec[i]*gamma) for i,expr in enumerate(exprs)]
                '''
                #func = lambda x: [(x[0]*triangle[0][i]+(1-x[0])*triangle[1][i])*x[1]+(1-x[1])*triangle[2][i]-(point[i]+vec[i]*x[2]) for i in range(3)]
                #solution = fsolve(func, (0,0,0))
                #alpha, beta, gamma = solution
                #if not np.allclose(func(solution), [0.0,0.0,0.0]):
                #    continue
                a = np.array([[triangle[j][i] for j in range(3)]+[-vec[i]] for i in range(3)] + [[1 for j in range(3)]+[0]]).astype('float')
                b = np.array([point[i] for i in range(3)] + [1]).astype('float')
                x = np.linalg.lstsq(a,b)[0]
                #print(x)
                if round_float(x[0]) >= 0 and round_float(x[0]) <= 1 and round_float(x[1]) >= 0 and round_float(x[1]) <= 1 and round_float(x[2]) >= 0 and round_float(x[2]) <= 1 and x[3] > 0:
                    p = tuple(sum(x[j]*triangle[j][i] for j in range(3)) for i in range(3))
                    if Polyhedron.inside_triangle(triangle,p):
                        output.append((x[3], p, face_index))
                        break
                '''
                try:
                    solutions = Polyhedron.timed_solve(eqs, dict=True)
                except TimeoutError:
                    continue
                if len(solutions):
                    alpha, beta, gamma = solutions[0][alpha], solutions[0][beta], solutions[0][gamma]
                #(beta*(alpha*triangle[0][i]+(1-alpha)*triangle[1][i])+(1-beta)*triangle[2][i]-position[i]+direction[i]*gamma)**2
                    if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1 and gamma > 0:
                        p = tuple(alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3))
                        p = tuple(beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(p))
                        #if Polyhedron.inside_triangle(triangle,p):
                        output.append((gamma, p, face_index))
                        break
                '''
        #print('is_inside', output, vec, point)
        return len(output)%2==1
    def round_verts(self, rounder):
        """Round vertices and clean up duplicate vertices, edges, and faces"""
        
        # Round all vertices
        for i in range(len(self.verts)):
            self.verts[i] = rounder(self.verts[i])
        
        # Remove duplicate vertices and update edge indices
        while True:
            double_break = False
            for i in range(len(self.verts)):
                for j in range(i + 1, len(self.verts)):
                    if distance(self.verts[i], self.verts[j]) < 1e-10:  # Essentially same vertex
                        # Remove vertex j
                        self.verts.pop(j)
                        
                        # Update all edge indices
                        k = 0
                        while k < len(self.edges):
                            edge = self.edges[k]
                            new_edge = set()
                            for vertex_idx in edge:
                                if vertex_idx == j:
                                    new_edge.add(i)  # Replace j with i
                                elif vertex_idx > j:
                                    new_edge.add(vertex_idx - 1)  # Shift down indices > j
                                else:
                                    new_edge.add(vertex_idx)  # Keep unchanged
                            
                            # Check if edge becomes degenerate (same vertex twice)
                            if len(new_edge) == 1:
                                # Remove degenerate edge
                                self.edges.pop(k)
                                
                                # Update face indices that reference this edge
                                for l in range(len(self.faces)):
                                    face = list(self.faces[l])
                                    new_face = []
                                    for edge_idx in face:
                                        if edge_idx > k:
                                            new_face.append(edge_idx - 1)  # Shift down
                                        elif edge_idx != k:  # Skip the removed edge
                                            new_face.append(edge_idx)
                                    self.faces[l] = set(new_face)
                            else:
                                self.edges[k] = new_edge
                                k += 1
                        
                        double_break = True
                        break
                if double_break:
                    break
            if not double_break:
                break
        
        # Remove duplicate edges
        while True:
            double_break = False
            for i in range(len(self.edges)):
                # Convert edge to frozenset of actual vertex positions for comparison
                try:
                    edge1_verts = frozenset(tuple(self.verts[w]) for w in self.edges[i])
                except IndexError:
                    # Edge references non-existent vertex, remove it
                    self.edges.pop(i)
                    # Update face indices
                    for k in range(len(self.faces)):
                        face = list(self.faces[k])
                        new_face = []
                        for edge_idx in face:
                            if edge_idx > i:
                                new_face.append(edge_idx - 1)
                            elif edge_idx != i:
                                new_face.append(edge_idx)
                        self.faces[k] = set(new_face)
                    double_break = True
                    break
                    
                for j in range(i + 1, len(self.edges)):
                    try:
                        edge2_verts = frozenset(tuple(self.verts[w]) for w in self.edges[j])
                    except IndexError:
                        # Edge references non-existent vertex, remove it
                        self.edges.pop(j)
                        # Update face indices
                        for k in range(len(self.faces)):
                            face = list(self.faces[k])
                            new_face = []
                            for edge_idx in face:
                                if edge_idx > j:
                                    new_face.append(edge_idx - 1)
                                elif edge_idx != j:
                                    new_face.append(edge_idx)
                            self.faces[k] = set(new_face)
                        double_break = True
                        break
                        
                    if edge1_verts == edge2_verts:
                        # Remove duplicate edge j
                        self.edges.pop(j)
                        
                        # Update face indices that reference edge j to use edge i instead
                        for k in range(len(self.faces)):
                            face = list(self.faces[k])
                            new_face = []
                            for edge_idx in face:
                                if edge_idx == j:
                                    new_face.append(i)  # Replace with i
                                elif edge_idx > j:
                                    new_face.append(edge_idx - 1)  # Shift down
                                else:
                                    new_face.append(edge_idx)  # Keep unchanged
                            self.faces[k] = set(new_face)
                        
                        double_break = True
                        break
                if double_break:
                    break
            if not double_break:
                break
        
        # Remove duplicate faces
        while True:
            double_break = False
            for i in range(len(self.faces)):
                try:
                    # Convert face to frozenset of frozensets of vertex positions
                    face1_edges = frozenset(
                        frozenset(tuple(self.verts[w]) for w in self.edges[z]) 
                        for z in self.faces[i]
                    )
                except IndexError:
                    # Face references non-existent edge, remove it
                    self.faces.pop(i)
                    double_break = True
                    break
                    
                for j in range(i + 1, len(self.faces)):
                    try:
                        face2_edges = frozenset(
                            frozenset(tuple(self.verts[w]) for w in self.edges[z]) 
                            for z in self.faces[j]
                        )
                    except IndexError:
                        # Face references non-existent edge, remove it
                        self.faces.pop(j)
                        double_break = True
                        break
                        
                    if face1_edges == face2_edges:
                        # Remove duplicate face
                        self.faces.pop(j)
                        double_break = True
                        break
                if double_break:
                    break
            if not double_break:
                break
        
        # Remove empty faces
        self.faces = [face for face in self.faces if len(face) > 0] 
    def dump(self, filename=None):
        output = ""
        output += "ply\n"
        output += "format ascii 1.0\n"
        output += "element vertex " + str(len(self.verts)) + "\n";\
        output += "property float x\n"
        output += "property float y\n"
        output += "property float z\n"
        output += "element edge " + str(len(self.edges)) + "\n"
        output += "property list uchar uint vertex_indices\n";
        output += "element face " + str(len(self.faces)) + "\n";
        output += "property list uchar uint edge_indices\n";
        output += "end_header\n";
        for vert in self.verts:
            output += ' '.join([str(x) for x in vert]) + "\n"
        for edge in self.edges:
            output += "2 "
            output += ' '.join([str(index) for index in edge]) + "\n"
        for face in self.faces:
            output += str(len(face)) + " "
            output += ' '.join([str(index) for index in face]) + "\n"
        if filename is not None:
            with open(filename, 'w') as f:
                f.write(output)
        return output
    # Updated version that can be used as both instance method and class method
    def loads(self, text=None):
        """Load polyhedron data from PLY format string
        
        Can be called as:
        - poly.loads(text)  # instance method
        - Polyhedron.loads(text)  # class method (returns new instance)
        """
        # If called as class method (self is actually the text string)
        if isinstance(self, str) and text is None:
            text = self
            poly = Polyhedron()  # Create new instance
        else:
            poly = self  # Use existing instance
        
        element_index = -1
        vert_count = 0
        edge_count = 0  
        face_count = 0
        elements_vector = []
        elements = {}
        
        # Clear existing data
        poly.verts = []
        poly.edges = []
        poly.faces = []
        
        lines = text.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
                
            split_line = line.split()
            
            # Check for end of header
            if len(split_line) >= 1 and split_line[0] == "end_header":
                element_index = 0
                continue
            
            # Parse header
            if element_index == -1:
                if len(split_line) >= 3 and split_line[0] == "element":
                    element_name = split_line[1]
                    elements_vector.append(element_name)
                    element_count = int(split_line[2])
                    
                    if element_name == "vertex":
                        vert_count = element_count
                    elif element_name == "edge":
                        edge_count = element_count
                    elif element_name == "face":
                        face_count = element_count
                        
                    # Initialize properties list for this element
                    elements[element_name] = []
                    
                elif len(split_line) >= 3 and split_line[0] == "property":
                    current_element = elements_vector[-1] if elements_vector else None
                    if current_element:
                        if split_line[1] == "float":
                            # Simple property like "property float x"
                            elements[current_element].append(split_line[2])
                        elif len(split_line) >= 5 and split_line[1:4] == ["list", "uchar", "uint"]:
                            # List property like "property list uchar uint vertex_indices"
                            elements[current_element].append(split_line[4])
            
            # Parse data
            elif element_index < len(elements_vector):
                element = elements_vector[element_index]
                
                if element == "vertex" and len(split_line) >= 3:
                    # Parse vertex coordinates
                    vert = [float(split_line[0]), float(split_line[1]), float(split_line[2])]
                    poly.verts.append(vert)
                    
                    if len(poly.verts) == vert_count:
                        element_index += 1
                        
                elif element == "edge" and len(split_line) >= 2:
                    # Parse edge indices
                    property_count = int(split_line[0])
                    if len(split_line) >= 1 + property_count:
                        edge = set()
                        for i in range(1, 1 + property_count):
                            edge.add(int(split_line[i]))
                        poly.edges.append(frozenset(edge))
                    
                    if len(poly.edges) == edge_count:
                        element_index += 1
                        
                elif element == "face" and len(split_line) >= 2:
                    # Parse face edge indices
                    property_count = int(split_line[0])
                    if len(split_line) >= 1 + property_count:
                        face = set()
                        for i in range(1, 1 + property_count):
                            face.add(int(split_line[i]))
                        poly.faces.append(frozenset(face))
                    
                    if len(poly.faces) == face_count:
                        element_index += 1
        
        # If called as class method, return the new instance
        if isinstance(self, str):
            return poly
    @classmethod
    def load(cls, filename):
        """Load polyhedron from PLY file"""
        with open(filename, 'r') as f:
            poly = cls()
            poly.loads(f.read())
            return poly 
    def circuit_overlap(circuits1, circuits2):
        circuit1 = Polyhedron.circuit_cut(Polyhedron.make_clockwise(circuits1))
        circuit2 = Polyhedron.circuit_cut(Polyhedron.make_clockwise(circuits2))
        if any(Polyhedron.intersect_segments(frozenset([circuit1[i-1],x]),frozenset([circuit2[j-1],y])) for i,x in enumerate(circuit1) for j,y in enumerate(circuit2)):
            return True
        if any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit1) for y in circuit2):
            return True
        if any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit2) for y in circuit1):
            return True
        return False
    def project_on_face_plane(self, face_index, point):
        circuit = Polyhedron.circuit_cut(Polyhedron.make_clockwise(self.circuits(face_index)))
        vec1 = tuple(circuit[1][i]-circuit[0][i] for i in range(3))
        vec2 = tuple(circuit[2][i]-circuit[0][i] for i in range(3))
        vec0 = cross3D(vec1,vec2)
        vec0 = tuple(vec0[i]/distance(vec0) for i in range(3))
        ''' 
        alpha, beta, gamma = symbols("alpha beta gamma")
        eqs = [Eq(point[i]+alpha*vec0[i],beta*vec1[i]+gamma*vec2[i]+circuit[0][i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        alpha = solutions[0][alpha]
        return tuple(point[i]+alpha*vec0[i] for i in range(3))
        '''
        a = np.array([[vec0[i], -vec1[i], -vec2[i]] for i in range(3)]).astype('float')
        b = np.array([circuit[0][i]-point[i] for i in range(3)]).astype('float')
        x = np.linalg.lstsq(a,b)[0]
        return tuple(point[i]+x[0]*vec0[i] for i in range(3))
    def project_ray_on_face_plane(self, face_index, ray):
        p1 = self.project_on_face_plane(face_index, ray[0])
        p2 = self.project_on_face_plane(face_index, ray[1])
        return p1, p2
    def project_ray_on_face_plane_and_intersect(self, face_index, ray):
        output = []
        p1, p2 = self.project_ray_on_face_plane(face_index, ray)
        #print(ray, {self.verts[index] for edge_index in self.faces[face_index] for index in self.edges[edge_index]}, p1, p2)
        if p1 == p2:
            return output
        vec = tuple(p2[i]-p1[i] for i in range(3))
        for edge_index in self.faces[face_index]:
            edge = [self.verts[x] for x in self.edges[edge_index]]
            '''
            alpha, beta = symbols("alpha beta")
            eqs = [Eq(alpha*vec[i]+p1[i],beta*edge[0][i]+(1-beta)*edge[1][i]) for i in range(3)]
            solutions = solve(eqs, dict=True)
            if len(solutions):
                try:
                    alpha = solutions[0][alpha]
                except KeyError:
                    continue
                try:
                    beta = solutions[0][beta]
                except KeyError:
                    continue
                if beta >= 0 and beta <= 1 and alpha >= 0:
                    output.append((alpha,tuple(alpha*vec[i]+p1[i] for i in range(3))))
        for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(self.circuits(face_index))):
            intersection = Polyhedron.intersect_triangle(ray, triangle)
            if intersection is not None:
                output.append(intersection)
            '''
            a = np.array([[-vec[i],edge[0][i],edge[1][i]] for i in range(3)]+[[0.0,1.0,1.0]]).astype('float')
            b = np.array([p1[0],p1[1],p1[2],1.0]).astype('float')
            x = np.linalg.lstsq(a,b)[0]
            dist = distance(tuple(x[0]*vec[i]+p1[i] for i in range(3)),tuple(x[1]*edge[0][i]+x[2]*edge[1][i] for i in range(3)))
            #print(dist, round_float(x[0]), round_float(x[1]))
            if dist <= 10**-5 and round_float(x[1]) >= 0 and round_float(x[1]) <= 1 and round_float(x[0]) >= 0:
                #output.append((x[0],tuple(x[0]*vec[i]+p1[i] for i in range(3))))
                output.append((x[0],tuple(x[1]*edge[0][i]+x[2]*edge[1][i] for i in range(3))))
        return output
    # For noncoplanar circuits
    def intersect_circuits(circuit1, circuit2):
        intersections = set()
        for i, p2 in enumerate(circuit1):
            p1 = circuit1[i-1]
            for triangle in Polyhedron.triangulate(circuit2):
                intersect = Polyhedron.intersect_triangle((p1,p2), triangle)
                if intersect:
                    intersections.add(intersect[1])
        for i, p2 in enumerate(circuit2):
            p1 = circuit1[i-1]
            for triangle in Polyhedron.triangulate(circuit1):
                intersect = Polyhedron.intersect_triangle((p1,p2), triangle)
                if intersect:
                    intersections.add(intersect[1])
        sorted_intersections = []
        if len(intersections):
            vec = tuple(intersections[1][i]-intersections[0][i] for i in range(3))
            for intersection in intersections:
                alpha = symbols("alpha")
                eqs = [Eq(alpha*vec[i]+intersections[0][i],intersection[i]) for i in range(3)]
                solutions = solve(eqs, dict=True)
                alpha = solutions[0][alpha]
                sorted_intersections.append((alpha,intersection))
        sorted_intersections = [x[1] for x in sorted(sorted_intersections)]
        edges = []
        for index, p1 in enumerate(sorted_intersections[:-1]):
            p2 = sorted_intersections[index+1]
            p = tuple((p1[i]+p2[i])/2 for i in range(3))
            if any(Polyhedron.inside_triangle(triangle, p) for triangle in circuit1) and any(Polyhedron.inside_triangle(triangle, p) for triangle in circuit2):
                edges.append(frozenset([p1, p2]))
        return sorted_intersections, edges
    def make_from_hull(hull, points):
        #print(hull.vertices)
        #print(hull.simplices)
        #print("points", points)
        hull_edges = set()
        simplices = [Polyhedron.make_clockwise([tuple(points[index] for index in simplex)])[0] for simplex in hull.simplices]
        #print(simplices)
        update = True
        while update:
            #print([len(x) for x in simplices])
            update = False
            for simplex1 in list(simplices):
                for simplex2 in list(simplices):
                    #if simplex1 != simplex2:
                    #    print(simplex1, simplex2, Polyhedron.coplanar(set(simplex1+simplex2)),len(set(simplex1)&set(simplex2)))
                    if simplex1 != simplex2 and Polyhedron.coplanar(set(simplex1+simplex2), rtol=0.01):
                        if len(set(simplex1)&set(simplex2)) > 1:
                            simplex3 = []
                            for i,x in enumerate(simplex1):
                                if simplex1[i-1] in simplex2 and x in simplex2:
                                    index1 = simplex2.index(simplex1[i-1])
                                    index2 = simplex2.index(x)
                                    #print(index1, index2)
                                    simplex3.extend(tuple(reversed([y for j,y in enumerate(simplex2) if j not in [index1, index2]])) + tuple([x]))
                                else:
                                    simplex3.append(x)
                            simplices = [simplex for simplex in simplices if simplex != simplex1 and simplex != simplex2] + [tuple(simplex3)]
                            update = True
                            break
                if update:
                    break
        faces = []
        for simplex in simplices:
            faces.append(set())
            for i, x in enumerate(simplex):
                edge = frozenset([simplex[i-1], x])
                faces[-1].add(edge)
                hull_edges.add(edge)
        double_break = True
        while double_break:
            double_break = False
            for edge1 in hull_edges:
                for edge2 in hull_edges:
                    #print(edge1, edge2, len(edge1&edge2) == 1 and Polyhedron.colinear(edge1|edge2))
                    if len(edge1&edge2) == 1 and Polyhedron.colinear(edge1|edge2, rtol=0.00001):
                        for point in edge1:
                            if not Polyhedron.point_on_segment(edge2, point):
                                break
                        else:
                            point = list(edge1-edge2)[0]
                            edge3 = None
                            double_break = True
                            break
                        for point in edge2:
                            if not Polyhedron.point_on_segment(edge1, point):
                                break
                        else:
                            point = list(edge2-edge1)[0]
                            edge3 = None
                            double_break = True
                            break
                        point = list(edge1&edge2)[0]
                        edge3 = edge1^edge2
                        double_break = True
                        break
                if double_break:
                    print('double break', edge1, edge2, edge3)
                    break
            if double_break:
                points.remove(point)
                if edge3 is not None:
                    hull_edges.add(edge3)
                    hull_edges.remove(edge1)
                    hull_edges.remove(edge2)
                    for face in faces:
                        if edge1 in face and edge2 in face:
                            face.add(edge3)
                        if edge1 in face:
                            face.remove(edge1)
                        if edge2 in face:
                            face.remove(edge2)
                to_delete = {edge for edge in hull_edges if point in edge}
                #print("to delete", to_delete)
                hull_edges = {edge for edge in hull_edges if point not in edge}
                for edge in to_delete:
                    to_merge = [face for face in faces if edge in face]
                    #print(to_merge)
                    faces.append(set())
                    for face in to_merge:
                        faces.remove(face)
                        faces[-1].update(face)
                        faces[-1].remove(edge)
                    if edge3 is not None:
                        faces[-1].add(edge3)
        '''
            if not double_break:
                crossing = set()
                for edge1 in hull_edges:
                    for edge2 in hull_edges:
                        intersect = Polyhedron.intersect_segments(edge1,edge2)
                        if intersect is not None:
                            if round_point(intersect) not in round_edge(edge1):
                                crossing.add(edge2)
                            if round_point(intersect) not in round_edge(edge2):
                                crossing.add(edge1)
                for edge in crossing:
                    hull_edges.remove(edge)
                    to_merge = [face for face in faces if edge in face]
                    faces.append(set())
                    for face in to_merge:
                        faces.remove(face)
                        faces[-1].update(face)
                    faces[-1].remove(edge)
        '''
                    
        #print("hull edges", hull_edges)
        poly = Polyhedron()
        poly.verts = [x for i,x in enumerate(points) if any(x in edge for edge in hull_edges)]
        poly.edges = [frozenset([round_point(vert) for vert in poly.verts].index(round_point(point)) for point in edge) for edge in hull_edges]
        poly.faces = [frozenset(poly.edges.index(frozenset([round_point(vert) for vert in poly.verts].index(round_point(point)) for point in edge)) for edge in face) for face in faces]
        print([len(face) for face in poly.faces])
        return poly
    def is_enclosed(self):
        for vert_index, vert in enumerate(self.verts):
            neighbors = set()
            for edge in self.edges:
                if vert_index in edge:
                    neighbors |= edge
            neighbors -= {vert_index}
            try:
                path = [list(neighbors)[0]]
            except IndexError:
                return True
            used_faces = set()
            while True:
                for face in self.faces:
                    if face not in used_faces and any(path[-1] in self.edges[edge_index] for edge_index in face):
                        triple_break = False
                        used_faces.add(face)
                        for edge_index in face:
                            for index in self.edges[edge_index]:
                                if index != path[-1] and index in neighbors:
                                    path.append(index)
                                    triple_break = True
                                    break
                            if triple_break:
                                break
                        if triple_break:
                            break
                else:
                    break
            if len(path) != len(neighbors)+1:
                return False
        return True
    def intersecting(self, other):
        polys = [self, other]
        for poly_index, poly in enumerate(polys):
            new_edges = set()
            other_poly = polys[poly_index-1]
            seen_verts = set()
            for vert_index, vert in enumerate(poly.verts):
                if vert_index in seen_verts:
                    continue
                if len(other_poly.in_faces(vert)):
                    continue
                leaves = []
                nonleaves = []
                root_in_poly = other_poly.is_inside(vert)
                if root_in_poly:
                    return True
                q = [vert_index]
                face_lookup = dict()
                while len(q):
                    current = q.pop()
                    if current in seen_verts:
                        continue
                    seen_verts.add(current)
                    edges = [edge for edge in poly.edges if current in edge]
                    for edge in edges:
                        v_i1 = current
                        v_i2 = list(edge-set([current]))[0]
                        segment = [poly.verts[v_i1],poly.verts[v_i2]]
                        intersects = sorted(other_poly.face_intersect(segment))
                        print("intersects",intersects, segment)
                        if len(set(x[1] for x in intersects)) > 1 or any(round_float(x[0]) != 0 and round_float(x[0]) != 1 for x in intersects):
                            return True
                        else:
                            q.append(v_i2)
        return False
    def circuit_helper(edges, start=None, previous=None, current=None, path=None, old_circuits=None):
        #print(path, previous, current)
        '''
        if path is not None:
            for x in old_circuits:
                if not len(set(path)-set(x)):
                    return set()
        '''
        edge_lookup = dict()
        for edge in edges:
            for point in edge:
                if point not in edge_lookup:
                    edge_lookup[point] = set()
                edge_lookup[point].add(edge)
        circuits = set()
        if start is None:
            seen = set()
            for start in edges:
                path = []
                point = list(start)[0]
                temp = edge_lookup[point] - set([start])
                for current in temp:
                    if Polyhedron.coplanar(set(path + [point])):
                        output = Polyhedron.circuit_helper(edges, start, start, current, path + [point], circuits)
                        circuits.update(output)
        else:
            if current == start:
                return {tuple(path)}
            point = list(current - previous)[0]
            if point in path:
                return set()
            if len(path) > 2:
                for i, p2 in enumerate(path[1:-2]):
                    p1 = path[i-1]
                    intersect = Polyhedron.intersect_segments([p1,p2],[path[-2],path[-1]])
                    if intersect is not None and distance(intersect, path[-2]) > 0.001 and distance(intersect, path[-1]) > 0.001:
                        return set()
            previous = current
            temp = list(edge_lookup[point]-set([previous]))
            for current in temp:
                if Polyhedron.coplanar(set(path + [point])):
                    circuits.update(Polyhedron.circuit_helper(edges, start, previous, current, path + [point], old_circuits|circuits))
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            for j,y in enumerate(circuits_list[i+1:]):
                if not len(set(x)-set(y)):
                    y_r = tuple(reversed(y))
                    if (y[y.index(x[0]):]+y[:y.index(x[0])])[:len(x)] == x or (y_r[y_r.index(x[0]):]+y_r[:y_r.index(x[0])])[:len(x)] == x:
                        if y in circuits:
                            circuits.remove(y)
        return circuits
    def intersect(self, other, use_cpp=True):
        def path_along_faces(p1, p2, face_path): 
            ray = other_poly.project_ray_on_face_plane(face_path[-1], (p1,p2))
            #print('ray', ray, (p1,p2), round_float((distance(ray[1],p2))), face_path)
            #if round_float(distance(ray[0],p1)) == 0 and round_float(distance(ray[1],p2)) == 0:
            if round_float(distance(ray[0],p1)) == 0 and round_float(distance(ray[1],p2)) == 0:
                #print('hey', [p1,p2])
                return [([p1,p2], face_path)]
            if ray[0] == ray[1]:
                return []
            #print(ray, current_face_index)
            #print("ray:", ray, "face:", {other_poly.verts[index] for edge_index in other_poly.faces[current_face_index] for index in other_poly.edges[edge_index]})
            intersections = sorted([x for x in other_poly.project_ray_on_face_plane_and_intersect(face_path[-1], (p1,p2)) if x[0] >= 0])
            #print('intersections',intersections)
            if len(intersections):
                #print('intersection', intersections[0])
                p = intersections[0][1]
                output = []
                if not poly.is_inside(p):
                    return []
                if p != p1:
                    for face_index2, face in enumerate(other_poly.faces):
                        print(face_index2,  any(Polyhedron.inside_triangle(triangle, p) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(Polyhedron.make_clockwise(other_poly.circuits(face_index2))))))
                        if face_index2 not in face_path and any(Polyhedron.inside_triangle(triangle, p) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(Polyhedron.make_clockwise(other_poly.circuits(face_index2))))):
                            print(path_along_faces(p, p2, face_path+[face_index2]))
                            output.extend([([p1]+x[0],x[1]) for x in path_along_faces(p, p2, face_path+[face_index2])])
                return output
            return []
        def path_around(points, polys=None):
            points = np.array(list(points), dtype=np.float64)
            print(points)
            n = len(points)
            if polys is None:
                polys = []
            # Compute distance matrix
            def create_distance_matrix(points):
                dist_matrix = {}
                for i in range(n):
                    dist_matrix[i] = {}
                    for j in range(n):
                        if i == j:
                            dist_matrix[i][j] = 0
                        else:
                            dist_matrix[i][j] = int(1000 * np.linalg.norm(points[i] - points[j]))  # scale to int
                            p1 = tuple(float(x) for x in points[i])
                            p2 = tuple(float(x) for x in points[j])
                            for poly in polys:
                                intersects = poly.face_intersect([p1,p2])
                                ps = []
                                #print(p1, p2, intersects, point, "path_around")
                                for intersect in sorted(intersects):
                                    if not len(ps) or intersect[1] != ps[-1]:
                                        ps.append(intersect[1])
                                ps = [p1] + ps + [p2]
                                for k, p in enumerate(ps[:-1]):
                                    midpoint = tuple(0.5*(p[l]+ps[k+1][l]) for l in range(3))        
                                    if not poly.is_inside(midpoint):
                                        dist_matrix[i][j] = sys.maxsize
                                        break
                                if dist_matrix[i][j] == sys.maxsize:
                                    break
                            ps = [tuple(float(x) for x in k/10*points[i] + (10-k)/10*points[j]) for k in range(1, 10)]
                            #print(points[i], points[j])
                            #midpoint = tuple(float(x) for x in 0.5*(points[i] + points[j]))
                            #if all(poly.is_inside(midpoint) for poly in polys):
                            if any(not poly.is_inside(point) for poly in polys for point in ps):
                                #print(points[i], points[j], [x for x in ps if any(not poly.is_inside(x) for poly in polys)], "path_around")
                                dist_matrix[i][j] = sys.maxsize
                return dist_matrix

            distance_matrix = create_distance_matrix(points)

            # Create the routing index manager and model
            manager = pywrapcp.RoutingIndexManager(n, 1, 0)  # 1 vehicle, depot = 0
            routing = pywrapcp.RoutingModel(manager)

            # Create a callback for the distance matrix
            def distance_callback(from_index, to_index):
                from_node = manager.IndexToNode(from_index)
                to_node = manager.IndexToNode(to_index)
                return distance_matrix[from_node][to_node]

            transit_callback_index = routing.RegisterTransitCallback(distance_callback)
            routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)
            #for i in range(n):
            #    for j in range(n):
            #        if polys is not None:
            #            midpoint = tuple(float(x) for x in 0.5*(points[i] + points[j]))
            #            if not all(poly.is_inside(midpoint) for poly in polys):
            #                a = manager.NodeToIndex(i)
            #                b = manager.NodeToIndex(j)
            #                routing.NextVar(a).RemoveValue(b)
            #                routing.NextVar(b).RemoveValue(a)
            # Force a closed tour (cycle)
            routing.AddDimension(
                transit_callback_index,
                0,  # no slack
                1_000_000,  # max distance
                True,  # start cumul to zero
                'Distance'
            )

            # Solve the problem
            search_parameters = pywrapcp.DefaultRoutingSearchParameters()
            search_parameters.first_solution_strategy = (
                routing_enums_pb2.FirstSolutionStrategy.SEQUENTIAL_CHEAPEST_INSERTION)
            #search_parameters.solution_limit = 1
            #search_parameters.time_limit.FromMilliseconds(1)
            #search_parameters.first_solution_strategy = (
            #    routing_enums_pb2.FirstSolutionStrategy.LOCAL_CHEAPEST_INSERTION)
            #search_parameters = pywrapcp.DefaultRoutingSearchParameters()
            #search_parameters.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.CHRISTOFIDES
            #search_parameters.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
            #search_parameters.time_limit.FromMilliseconds(10)
            #search_parameters.time_limit.FromMilliseconds(1000)
            #search_parameters.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.TABU_SEARCH
            start_time = time.perf_counter()
            solution = routing.SolveWithParameters(search_parameters)
            print('TSP time:', time.perf_counter()-start_time)
            # Extract solution
            def get_ordered_points(solution, manager, routing, points):
                index = routing.Start(0)
                ordered_indices = []
                while not routing.IsEnd(index):
                    ordered_indices.append(manager.IndexToNode(index))
                    try:
                        index = solution.Value(routing.NextVar(index))
                    except AttributeError:
                        return None, None
                ordered_points = points[ordered_indices]
                return ordered_points, ordered_indices

            ordered_points, ordered_indices = get_ordered_points(solution, manager, routing, points)
            print(ordered_indices)
            if ordered_points is None:
                return None
            return [tuple(float(y) for y in x) for x in ordered_points]
        '''
        def path_around(points):
            points = tuple(points)
            def write_tsplib_2d(filename, x, y, name="problem"):
                with open(filename, "w") as f:
                    f.write(f"NAME: {name}\n")
                    f.write("TYPE: TSP\n")
                    f.write(f"DIMENSION: {len(x)}\n")
                    f.write("EDGE_WEIGHT_TYPE: EUC_2D\n")
                    f.write("NODE_COORD_SECTION\n")
                    for i, (xi, yi) in enumerate(zip(x, y), start=1):
                        f.write(f"{i} {xi} {yi}\n")
                    f.write("EOF\n")
            def read_tsplib_2d(filename):
                with open(filename) as f:
                    f.readline()
                    return [int(x) for x in f.readline().split()]
            planar = Polyhedron.make_planar(points)
            X = [x[0] for x in points]
            Y = [x[1] for x in points]
            write_tsplib_2d("problem.tsp", x, y)
            subprocess.run(["concorde", "-x", "problem.tsp"])
            indices = read_tsplib_2d("problem.sol")
            return [points[x] for x in indices]
        '''
        def path_around2(points, vert):
            centroid = tuple(sum(x[i] for x in points)/len(points) for i in range(3))
            vec0 = tuple(centroid[i]-vert[i] for i in range(3))
            gamma = symbols("gamma")
            eqs = [Eq(vec0[0]*vec0[1]+vec0[1]*-vec0[0]+gamma*vec0[2],0)]
            solutions = solve(eqs, dict=True)
            vec1 = (vec0[1], -vec0[0], solutions[0][gamma])
            vec2 = cross3D(vec0, vec1)
            vec1 = tuple(vec1[i]/distance(vec1) for i in range(3))
            vec2 = tuple(vec2[i]/distance(vec2) for i in range(3))
            projection = []
            for point in points:
                alpha, beta, gamma = symbols("alpha beta gamma")
                eqs = [Eq(alpha*vec0[i]+point[i],beta*vec1[i]+gamma*vec2[i]) for i in range(3)]
                solutions = solve(eqs, dict=True)
                alpha = solutions[0][alpha]
                projection.append(tuple(alpha*vec0[i]+point[i] for i in range(3)))
            planar_projection = Polyhedron.make_planar(projection)
            reverse_projection = dict(zip(planar_projection, points))
            point = min(planar_projection)
            path = []
            rotated_mapping = {x:x for x in planar_projection}
            reverse_mapping = {x:x for x in planar_projection}
            angle = (0,0,0)
            while (not len(path) or point != path[0]) and point is not None:
                path.append(point)
                keys = []
                temp_keys = []
                for key in rotated_mapping:
                    keys.append(key)
                    temp_keys.append(tuple(key[i]-reverse_mapping[path[-1]][i] for i in range(3)))
                if len(path) > 1:
                    angle = (0,0,angle[2]+atan2(reverse_mapping[path[-1]][1],reverse_mapping[path[-1]][0]))
                    temp_keys = rotate(temp_keys, angle)
                rotated_mapping = {temp_keys[i]:rotated_mapping[x] for i, x in enumerate(keys)}
                reverse_mapping = {rotated_mapping[key]:key for key in rotated_mapping}
                gift_wrap_angle_distance = []
                gift_wrap = []
                for key in rotated_mapping:
                    print('key', key)
                    if key == (0,0,0):
                        continue
                    if np.allclose(atan2(key[1],key[0]),pi, rtol=0.00001):
                        continue
                    gift_wrap_angle_distance.append((atan2(key[1],key[0]), distance(key)))
                    gift_wrap.append(rotated_mapping[key])
                maxi = ((-float("inf"),-float("inf")), None)
                for i,x in enumerate(gift_wrap_angle_distance):
                    print(path)
                    if x > maxi[0] and gift_wrap[i] not in path:# and (maxi[0][0] == -float("inf") or round_float(maxi[0][0] > -round_float(pi/2))):
                        maxi = (x, gift_wrap[i])
                print('maxi', maxi)
                point = maxi[1]
            return [reverse_projection[x] for x in path]

        def path_around3(points, polys = None):
            mini = (float('inf'), None)
            for ordering in permutations(points):
                break_continue = False

                for i,p2 in enumerate(ordering):
                    midpoint = tuple(ordering[i-1][j]+p2[j] for j in range(3))
                    if not all(poly.is_inside(midpoint) for poly in polys):
                        break_continue = True
                        break
                if break_continue:
                    continue
                perimeter = sum(distance(ordering[i-1],x) for i,x in enumerate(ordering))
                if perimeter < mini[0]:
                    mini = (perimeter, ordering)
                '''
                double_break = False
                for i,p2 in enumerate(ordering):
                    p1 = ordering[i-1]
                    for j in range(i+1, len(ordering)):
                        p4 = ordering[j]
                        p3 = ordering[j-1]
                        intersect = Polyhedron.intersect_segments(frozenset([p1,p2]), frozenset([p3,p4]))
                        if intersect is not None and intersect not in {p1,p2,p3,p4}:
                            double_break = True
                    if double_break:
                        break
                else:
                    return ordering
                '''
            return mini[1]
        edge_sets_per_poly = [[],[]]
        polys = [self, other]
        if use_cpp:
            self.dump("poly1.ply")
            other.dump("poly2.ply")
            subprocess.run(["./intersect_polyhedron"])
            blank = False
            new_poly = False
            poly_index = 0
            index = -1
            with open("out.txt") as f:
                for line in f:
                    if not(len(line.strip())):
                        blank = True
                    elif line.strip() == "}":
                        new_poly = True
                    else:
                        if new_poly:
                            poly_index += 1
                            if poly_index > 1:
                                break
                            index = -1
                            new_poly = False
                        if index == -1:
                            edge_sets_per_poly[poly_index].append(set())
                            index += 1
                        elif blank:
                            edge_sets_per_poly[poly_index].append(set())
                            blank = False
                        try:
                            p1_text, p2_text = line.strip().split(" ")
                            p1 = tuple(float(x) for x in p1_text.strip('[').strip(']').split(','))
                            p2 = tuple(float(x) for x in p2_text.strip('[').strip(']').split(','))
                            edge_sets_per_poly[poly_index][-1].add(frozenset([p1,p2]))
                        except ValueError:
                            pass
            print(edge_sets_per_poly)
        else:
            for poly_index, poly in enumerate(polys):
                other_poly = polys[poly_index-1]
                edges_are_subset = True
                for edge1 in poly.edges:
                    edge_copresent = False
                    for edge2 in other_poly.edges:
                        edge2 = frozenset(other_poly.verts[index] for index in edge2)
                        if all(Polyhedron.point_on_segment(edge2, poly.vert[index]) for index in edge1):
                            edge_copresent = True
                            break
                    if not edge_copresent:
                        edge_are_subset = False
                        break
                if edges_are_subset:
                    return Polyhedron()
            for poly_index, poly in enumerate(polys):
                new_edges = set()
                other_poly = polys[poly_index-1]
                seen_verts = set()
                seen_leaves = set()
                for vert_index, vert in enumerate(poly.verts):
                    #if vert_index in seen_verts:
                    #    continue
                    if len(other_poly.in_faces(vert)):
                        continue
                    leaves = []
                    nonleaves = []
                    root_in_poly = other_poly.is_inside(vert)
                    q = [vert_index]
                    face_lookup = dict()
                    print('vert',vert)
                    while len(q):
                        current = q.pop()
                        if current in seen_verts:
                            continue
                        seen_verts.add(current)
                        nonleaves.append(current)
                        edges = [edge for edge in poly.edges if current in edge]
                        for edge in edges:
                            v_i1 = current
                            v_i2 = list(edge-set([current]))[0]
                            segment = [poly.verts[v_i1],poly.verts[v_i2]]
                            intersects = sorted(other_poly.face_intersect(segment))
                            print("intersects",intersects, segment)
                            in_faces = other_poly.in_faces(poly.verts[v_i2])
                            if len(in_faces):
                                leaves.append(poly.verts[v_i2])
                                if poly.verts[v_i2] not in face_lookup:
                                    face_lookup[poly.verts[v_i2]] = []
                                face_lookup[poly.verts[v_i2]].extend(in_faces)
                                q.append(v_i2)
                            elif len(intersects):
                                print("alpha", intersects[0][0])
                                leaves.append(intersects[0][1])
                                for index, intersect in enumerate(intersects):
                                    if intersect[0] == intersects[0][0]:
                                        if intersect[1] not in face_lookup:
                                            face_lookup[intersect[1]] = []
                                        face_lookup[intersect[1]].append(intersect[2])
                                    else:
                                        break
                                if root_in_poly:
                                    new_edges.add(frozenset([poly.verts[v_i1], intersects[0][1]]))
                                else:
                                    for intersect in intersects[1:]:
                                        if intersect[0] > intersects[0][0]:
                                            new_edges.add(frozenset([intersects[0][1], intersect[1]]))
                                            break
                            else:
                                q.append(v_i2)
                                if root_in_poly:
                                    new_edges.add(frozenset([poly.verts[v_i1], poly.verts[v_i2]]))
                    print('leaves',leaves)
                    print(face_lookup)
                    '''
                    if not len(leaves) and root_in_poly:
                        output = []
                        new_poly = Polyhedron
                        new_poly.verts = list(poly.verts)
                        new_poly.edges = list(poly.edges)
                        new_poly.faces = list(poly.faces)
                        output.append(new_poly)
                        #return output
                    '''
                    #print(vert, nonleaves, leaves)
                    if len(leaves) > 2 and frozenset([round_point(x) for x in leaves]) not in seen_leaves:
                        seen_leaves.add(frozenset([round_point(x) for x in leaves]))
                        print(leaves)
                        path = path_around2(leaves, vert)
                        print('path',path)
                        for index, p2 in enumerate(path):
                            p1 = path[index-1]
                            #print("p1 p2", p1, p2)
                            starting_faces = face_lookup[p1]
                            print(p1, starting_faces)
                            #print("starting_faces",starting_faces)
                            for face_index1 in starting_faces:
                                paths = path_along_faces(p1,p2,[face_index1])
                                print(len(paths))
                                for point_path, face_path in paths:
                                    print('paths:', point_path, face_path)
                                    if len(point_path):
                                        for i,x in enumerate(point_path[:-1]):
                                            new_edges.add(frozenset([x,point_path[i+1]]))
                                if not len(paths):
                                    new_edges.add(frozenset([p1,p2]))
                #print('leaves', seen_leaves)
                '''
                if len(new_edges) == 0 and all(len(other_poly.in_faces(vert)) for vert in poly.verts):
                    output = []
                    new_poly = Polyhedron
                    new_poly.verts = list(poly.verts)
                    new_poly.edges = list(poly.edges)
                    new_poly.faces = list(poly.faces)
                    output.append(new_poly)
                    #return output
                '''
                edge_sets_per_poly[poly_index].append(new_edges)
        print('edge sets per poly:', edge_sets_per_poly)
        #polyhedrons_per_poly = [[],[]]
        #first = True
        edges_per_poly = [set(), set()]
        for poly_index, edge_sets in enumerate(edge_sets_per_poly):
            print(edge_sets)
            for edges in edge_sets:
                edges_per_poly[poly_index] |= edges
        points_per_poly = [{point for edge in edges for point in edge} for edges in edges_per_poly]
        points = points_per_poly[0] | points_per_poly[1]
        print('POINTS', points)
        point_map = dict()
        for point in points:
            if round_point(point) in point_map:
                if distance(point, round_point(point)) < distance(point_map[round_point(point)], round_point(point)):
                    point_map[round_point(point)] = point
            else:
                point_map[round_point(point)] = point
        points = set(point_map.values())
        print(points)
        poly = Polyhedron()
        poly.verts = list(points)
        seen_components = set()
        for poly_index in range(2):
            for face_index, face in enumerate(polys[poly_index].faces):
                #if any(all(any(Polyhedron.inside_triangle(triangle,polys[(poly_index+1)%2].verts[index]) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(polys[poly_index].circuits(face_index)))) for edge_index in polys[(poly_index+1)%2].faces[face_index2] for index in polys[(poly_index+1)%2].edges[edge_index]) for face_index2 in range(len(polys[(poly_index+1)%2].faces))):
                #    continue
                cofacial_points = set()
                for point in points:
                    if any(Polyhedron.inside_triangle(triangle,point) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(Polyhedron.make_clockwise(polys[poly_index].circuits(face_index))))):
                        cofacial_points.add(point)
                cofacial_points = list(cofacial_points)
                print("cofacial_points", cofacial_points)
                if len(cofacial_points) > 2:
                    edges = set()
                    print(len(cofacial_points))
                    for index1, point1 in enumerate(cofacial_points):
                        for index2, point2 in enumerate(cofacial_points[index1+1:]):
                            double_break = False
                            for p in polys:
                                intersects = p.face_intersect([point1,point2])
                                ps = []
                                for intersect in sorted(intersects):
                                    if not len(ps) or intersect[1] != ps[-1]:
                                        ps.append(intersect[1])
                                for k, point in enumerate(ps):
                                    if distance(point, point1) < 0.001:
                                        ps[k] = point1
                                    if distance(point, point2) < 0.001:
                                        ps[k] = point2
                                ps = [point1] + ps + [point2]
                                for k, point in enumerate(ps[:-1]):
                                    midpoint = tuple(0.5*(point[l]+ps[k+1][l]) for l in range(3))
                                    if not p.is_inside(midpoint):
                                        print(ps)
                                        double_break = True
                                        break
                                if double_break:
                                    break
                            if double_break:
                                continue
                            for p in polys:
                                print(point1, p.in_faces(point1), point2, p.in_faces(point2))
                                if any(not Polyhedron.coplanar(set(cofacial_points)|{p.verts[index] for edge_index in p.faces[i] for index in p.edges[edge_index]}) for i in set(p.in_faces(point1))&set(p.in_faces(point2))):
                                    edge = frozenset([point1,point2])
                                    if not any(point != point1 and point != point2 and Polyhedron.point_on_segment(edge, point) for point in cofacial_points):
                                        edges.add(edge)
                                        break
                    edges = list(edges)
                    '''
                    updated = True
                    while updated:
                        updated = False
                        for edge_index1, edge1 in enumerate(edges):
                            for edge_index2, edge2 in enumerate(edges):
                                if len(edge1&edge2) == 1 and Polyhedron.colinear(edge1|edge2):
                                    edges.remove(edge1)
                                    edges.remove(edge2)
                                    edges.append(edge1^edge2)
                                    updated = True
                                    break
                            if updated:
                                break
                    '''
                    print(edges)
                    for circuit in Polyhedron.circuit_helper(edges):
                        new_face = set()
                        #circuit = tuple(x for i,x in enumerate(circuit) if not Polyhedron.colinear([circuit[i-1],x,circuit[(i+1)%len(circuit)]]))
                        if not len(circuit):
                            continue
                        for i,p2 in enumerate(circuit):
                            p1 = circuit[i-1]
                            print(p1, p2)
                            p1 = poly.verts.index(p1)
                            p2 = poly.verts.index(p2)
                            edge = frozenset([p1,p2])
                            try:
                                edge_index = poly.edges.index(edge)
                                new_face.add(edge_index)
                            except ValueError:
                                poly.edges.append(edge)
                                new_face.add(len(poly.edges)-1)
                        print(circuit)
                        #for face_index, face in enumerate(poly.faces):
                        #    pass
                            #if Polyhedron.find_exterior_circuit(poly.circuits(face_index)|{circuit}):
                            #    poly.faces[face_index] |= new_face
                            #    break
                            
                        #else:
                        poly.faces.append(frozenset(new_face))
        print('after cofacial', len(poly.verts), len(poly.edges), len(poly.faces))
        edge_index = 0
        while edge_index < len(poly.edges):
            if sum(edge_index in face for face in poly.faces) < 2:
                for face_index, face in enumerate(poly.faces):
                    poly.faces[face_index] = frozenset(x-1 if x > edge_index else x for x in face if x != edge_index)
                del poly.edges[edge_index]
            else:
                edge_index += 1
        updated = True
        while updated:
            updated = False
            for face_index1, face1 in enumerate(poly.faces):
                face1 = frozenset(frozenset(poly.verts[index] for index in poly.edges[edge_index]) for edge_index in poly.faces[face_index1])
                for face_index2, face2 in enumerate(poly.faces):
                    face2 = frozenset(frozenset(poly.verts[index] for index in poly.edges[edge_index]) for edge_index in poly.faces[face_index2])
                    double_break = False
                    if face_index1 < face_index2 and Polyhedron.coplanar({point for edge in face1 for point in edge}|{point for edge in face2 for point in edge}):
                        for edge1 in face1:
                            for edge2 in face2:
                                if sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 2:
                                    poly.faces[face_index1] |= poly.faces[face_index2]
                                    del poly.faces[face_index2]
                                    updated = True
                                    break
                                elif sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 2:
                                    poly.faces[face_index1] |= poly.faces[face_index2]
                                    del poly.faces[face_index2]
                                    updated = True
                                    break
                                elif sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 1 and sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 1 and len(edge1&edge2) == 0:
                                    poly.faces[face_index1] |= poly.faces[face_index2]
                                    del poly.faces[face_index2]
                                    updated = True
                                    break
                                elif len(edge1&edge2) == 1 and Polyhedron.colinear(edge1|edge2):
                                    poly.faces[face_index1] |= poly.faces[face_index2]
                                    del poly.faces[face_index2]
                                    updated = True
                                    break
                            if updated:
                                print(face1, face2)
                                break
                    if updated:
                        break
                if updated:
                    break
        edge_index = 0
        while edge_index < len(poly.edges):
            if sum(edge_index in face for face in poly.faces) < 2:
                for face_index, face in enumerate(poly.faces):
                    poly.faces[face_index] = frozenset(x-1 if x > edge_index else x for x in face if x != edge_index)
                del poly.edges[edge_index]
            else:
                edge_index += 1
        print('after cofacial2', len(poly.verts), len(poly.edges), len(poly.faces))
        faces = []
        for face_index, face in enumerate(poly.faces):
            print(frozenset(frozenset(poly.verts[index] for index in poly.edges[edge_index]) for edge_index in face))
            circuits = poly.circuits(face_index)
            circuit_sets = [{circuit} for circuit in circuits]
            index1 = 0
            print('CIRCUITS', circuits)
            while index1 < len(circuit_sets):
                index2 = index1 + 1
                while index2 < len(circuit_sets):
                    print(circuit_sets[index1])
                    print(circuit_sets[index2])
                    if Polyhedron.find_exterior_circuit({tuple(x for i,x in enumerate(circuit) if not Polyhedron.colinear([circuit[i-1],x,circuit[(i+1)%len(circuit)]])) for circuit in circuit_sets[index1]}|{tuple(x for i,x in enumerate(circuit) if not Polyhedron.colinear([circuit[i-1],x,circuit[(i+1)%len(circuit)]])) for circuit in circuit_sets[index2]}) is not None:
                        circuit_sets[index1] |= circuit_sets[index2]
                        del circuit_sets[index2]
                    else:
                        index2 += 1
                index1 += 1
            print("CIRCUIT SETS", circuit_sets)
            new_faces = []
            for circuit_set in circuit_sets:
                face = set()
                for circuit in circuit_set:
                    for i, p2 in enumerate(circuit):
                        p1 = circuit[i-1]
                        face.add(frozenset([p1,p2]))
                new_faces.append(list(face))
            index1 = 0
            while index1 < len(new_faces):
                index2 = index1 + 1
                while index2 < len(new_faces):
                    combine = False
                    updated = True
                    while updated:
                        updated = False
                        edge_index1 = 0
                        while edge_index1 < len(new_faces[index1]):
                            edge_index2 = 0
                            while edge_index2 < len(new_faces[index2]):
                                edge1 = new_faces[index1][edge_index1]
                                edge2 = new_faces[index2][edge_index2]
                                if sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 2:
                                    print(1, edge2, edge1)
                                    del new_faces[index1][edge_index1]
                                    del new_faces[index2][edge_index2]
                                    p1, p2 = edge1
                                    edge3 = [p1, sorted([(distance(p1,x),x) for x in edge2])[0][1]]
                                    if round_float(distance(edge3[0],edge3[1])) > 0:
                                        edge3 = frozenset(edge3)
                                        new_faces[index1].append(edge3)
                                        edge3 = frozenset(poly.verts.index(point) for point in edge3)
                                        if edge3 not in poly.edges:
                                            poly.edges.append(edge3)
                                    edge3 = [p2, sorted([(distance(p2,x),x) for x in edge2])[0][1]]
                                    if round_float(distance(edge3[0],edge3[1])) > 0:
                                        edge3 = frozenset(edge3)
                                        new_faces[index1].append(edge3)
                                        edge3 = frozenset(poly.verts.index(point) for point in edge3)
                                        if edge3 not in poly.edges:
                                            poly.edges.append(edge3)
                                    combine = True
                                    updated = True
                                    break
                                elif sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 2:
                                    print(2, edge2, edge1)
                                    del new_faces[index1][edge_index1]
                                    del new_faces[index2][edge_index2]
                                    p1, p2 = edge2
                                    edge3 = [p1, sorted([(distance(p1,x),x) for x in edge1])[0][1]]
                                    if round_float(distance(edge3[0],edge3[1])) > 0:
                                        edge3 = frozenset(edge3)
                                        new_faces[index2].append(edge3)
                                        edge3 = frozenset(poly.verts.index(point) for point in edge3)
                                        if edge3 not in poly.edges:
                                            poly.edges.append(edge3)
                                    edge3 = [p2, sorted([(distance(p2,x),x) for x in edge1])[0][1]]
                                    if round_float(distance(edge3[0],edge3[1])) > 0:
                                        edge3 = frozenset(edge3)
                                        new_faces[index2].append(edge3)
                                        edge3 = frozenset(poly.verts.index(point) for point in edge3)
                                        if edge3 not in poly.edges:
                                            poly.edges.append(edge3)
                                    combine = True
                                    updated = True
                                    break
                                else:
                                    edge_index2 += 1
                            if updated:
                                break
                            edge_index1 += 1
                    if combine:
                        new_faces[index1] += new_faces[index2]
                        del new_faces[index2]
                    else:
                        index2 += 1
                index1 += 1
            index1 = 0
            while index1 < len(new_faces):
                index2 = index1 + 1
                while index2 < len(new_faces):
                    combine = False
                    updated = True
                    while updated:
                        updated = False
                        edge_index1 = 0
                        while edge_index1 < len(new_faces[index1]):
                            edge_index2 = 0
                            while edge_index2 < len(new_faces[index2]):
                                edge1 = new_faces[index1][edge_index1]
                                edge2 = new_faces[index2][edge_index2]
                                if sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 1 and sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 1 and len(edge1&edge2) == 0:
                                    print(3, edge2, edge1)
                                    del new_faces[index1][edge_index1]
                                    del new_faces[index2][edge_index2]
                                    p1, p2 = edge1
                                    p3, p4 = edge2
                                    if Polyhedron.point_on_segment(edge2, p1):
                                        p1, p2 = p2, p1
                                    if Polyhedron.point_on_segment(edge1, p4):
                                        p3, p4 = p4, p3
                                    edge3 = [p1, p3]
                                    if round_float(distance(edge3[0],edge3[1])) > 0:
                                        edge3 = frozenset(edge3)
                                        new_faces[index1].append(edge3)
                                        edge3 = frozenset(poly.verts.index(point) for point in edge3)
                                        if edge3 not in poly.edges:
                                            poly.edges.append(edge3)
                                    edge3 = [p2, p4]
                                    if round_float(distance(edge3[0],edge3[1])) > 0:
                                        edge3 = frozenset(edge3)
                                        new_faces[index2].append(edge3)
                                        edge3 = frozenset(poly.verts.index(point) for point in edge3)
                                        if edge3 not in poly.edges:
                                            poly.edges.append(edge3)
                                    combine = True
                                    updated = True
                                    break
                                else:
                                    edge_index2 += 1
                            if updated:
                                break
                            edge_index1 += 1
                    if combine:
                        new_faces[index1] += new_faces[index2]
                        del new_faces[index2]
                    else:
                        index2 += 1
                index1 += 1
            faces.extend([frozenset(poly.edges.index(frozenset(poly.verts.index(point) for point in edge)) for edge in face) for face in new_faces])
        poly.faces = faces
        print('FACES', faces)
        for face_index, face in enumerate(poly.faces):
            face = list(face)
            edge_index1 = 0
            while edge_index1 < len(face):
                edge_index2 = edge_index1 + 1
                while edge_index2 < len(face):
                    edge1 = frozenset(poly.verts[index] for index in poly.edges[face[edge_index1]])
                    edge2 = frozenset(poly.verts[index] for index in poly.edges[face[edge_index2]])
                    if len(edge1&edge2) == 1 and Polyhedron.colinear(edge1|edge2):
                        edge3 = poly.edges[face[edge_index1]]^poly.edges[face[edge_index2]]
                        if edge3 not in poly.edges:
                            poly.edges.append(edge3)
                        face[edge_index1] = poly.edges.index(edge3)
                        del face[edge_index2]
                    else:
                        edge_index2 += 1
                edge_index1 += 1
            poly.faces[face_index] = frozenset(face)
        edge_index = 0
        while edge_index < len(poly.edges):
            if sum(edge_index in face for face in poly.faces) < 2:
                for face_index, face in enumerate(poly.faces):
                    poly.faces[face_index] = frozenset(x-1 if x > edge_index else x for x in face if x != edge_index)
                del poly.edges[edge_index]
            else:
                edge_index += 1
        vert_index = 0
        while vert_index < len(poly.verts):
            if all(vert_index not in edge for edge in poly.edges):
                del poly.verts[vert_index]
                for edge_index, edge in enumerate(poly.edges):
                    poly.edges[edge_index] = frozenset(x-1 if x > vert_index else x for x in edge)
            else:
                vert_index += 1
        poly.faces = [face for face in poly.faces if len(face)]
        return poly

    def add_subtract_helper(poly1, poly2):
        pairings = []
        faces1 = [frozenset(frozenset(poly1.verts[index] for index in poly1.edges[edge_index]) for edge_index in face) for face in poly1.faces]
        faces2 = [frozenset(frozenset(poly2.verts[index] for index in poly2.edges[edge_index]) for edge_index in face) for face in poly2.faces]
        for face_index1, face1 in enumerate(faces1):
            face1 = frozenset(frozenset(poly1.verts[index] for index in poly1.edges[edge_index]) for edge_index in poly1.faces[face_index1])
            for face_index2, face2 in enumerate(faces2):
                face2 = frozenset(frozenset(poly2.verts[index] for index in poly2.edges[edge_index]) for edge_index in poly2.faces[face_index2])
                double_break = False
                if Polyhedron.coplanar({point for edge in face1 for point in edge}|{point for edge in face2 for point in edge}):
                    for edge1 in face1:
                        for edge2 in face2:
                            if sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 2:
                                pairings.append((face_index1, face_index2))
                                double_break = True
                                break
                            elif sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 2:
                                pairings.append((face_index1, face_index2))
                                double_break = True
                                break
                            elif sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 1 and sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 1 and len(edge1&edge2) == 0:
                                pairings.append((face_index1, face_index2))
                                double_break = True
                                break
                        if double_break:
                            break
        print(pairings)
        faces = [faces1, faces2]
        maps = [{i:[] for i in range(len(x))} for x in faces]
        for v1, v2 in pairings:
            maps[0][v1].append(v2)
            maps[1][v2].append(v1)
        visited = set()
        new_faces = []
        for poly_index in range(2):
            for face_index in range(len(faces[poly_index])):
                if (poly_index, face_index) not in visited:
                    face = []
                    face.extend(faces[poly_index][face_index])
                    q = [(poly_index, face_index)]
                    while len(q):
                        indices = q.pop()
                        print(indices)
                        visited.add(indices)
                        p_i1, f_i1 = indices
                        p_i2 = (p_i1+1)%2
                        for f_i2 in maps[p_i1][f_i1]:
                            if (p_i2,f_i2) in visited:
                                continue
                            q.append((p_i2,f_i2))
                            for edge2 in faces[p_i2][f_i2]:
                                edge_index1 = 0
                                new_edges = []
                                combine = False
                                while edge_index1 < len(face):
                                    edge1 = face[edge_index1]
                                    if sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 2:
                                        print(1, edge2, edge1)
                                        del face[edge_index1]
                                        p1, p2 = edge1
                                        edge3 = [p1, sorted([(distance(p1,x),x) for x in edge2])[0][1]]
                                        if round_float(distance(edge3[0],edge3[1])) > 0:
                                            edge3 = frozenset(edge3)
                                            new_edges.append(edge3)
                                        edge3 = [p2, sorted([(distance(p2,x),x) for x in edge2])[0][1]]
                                        if round_float(distance(edge3[0],edge3[1])) > 0:
                                            edge3 = frozenset(edge3)
                                            new_edges.append(edge3)
                                        combine = True
                                    elif sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 2:
                                        print(2, edge2, edge1)
                                        del face[edge_index1]
                                        p1, p2 = edge2
                                        edge3 = [p1, sorted([(distance(p1,x),x) for x in edge1])[0][1]]
                                        if round_float(distance(edge3[0],edge3[1])) > 0:
                                            edge3 = frozenset(edge3)
                                            new_edges.append(edge3)
                                        edge3 = [p2, sorted([(distance(p2,x),x) for x in edge1])[0][1]]
                                        if round_float(distance(edge3[0],edge3[1])) > 0:
                                            edge3 = frozenset(edge3)
                                            new_edges.append(edge3)
                                        combine = True
                                    elif sum(Polyhedron.point_on_segment(edge1, x) for x in edge2) == 1 and sum(Polyhedron.point_on_segment(edge2, x) for x in edge1) == 1:
                                        if len(edge1&edge2) == 0:
                                            print(3, edge2, edge1)
                                            del face[edge_index1]
                                            p1, p2 = edge1
                                            p3, p4 = edge2
                                            if Polyhedron.point_on_segment(edge2, p1):
                                                p1, p2 = p2, p1
                                            if Polyhedron.point_on_segment(edge1, p4):
                                                p3, p4 = p4, p3
                                            edge3 = [p1, p3]
                                            if round_float(distance(edge3[0],edge3[1])) > 0:
                                                edge3 = frozenset(edge3)
                                                new_edges.append(edge3)
                                            edge3 = [p2, p4]
                                            if round_float(distance(edge3[0],edge3[1])) > 0:
                                                edge3 = frozenset(edge3)
                                                new_edges.append(edge3)
                                            combine = True
                                        elif len(edge1&edge2) == 1 and Polyhedron.colinear(edge1|edge2):
                                            print(4, edge2, edge1)
                                            del face[edge_index1]
                                            if edge2 in new_edges:
                                                new_edges.remove(edge2)
                                            new_edges.append(edge1^edge2)
                                            edge2 = edge1^edge2
                                            combine = True
                                        else:
                                            edge_index1 += 1
                                    else:
                                        edge_index1 += 1
                                if combine:
                                    face.extend(new_edges)
                                else:
                                    face.append(edge2)
                            print()
                    print()
                    face = frozenset(face)
                    circuits = Polyhedron.circuit_helper(face)
                    #circuits = {circuit for circuit in circuits if not any(Polyhedron.colinear([circuit[i-1],x,circuit[(i+1)%len(circuit)]]) and not Polyhedron.point_on_segment([circuit[i-1],circuit[(i+1)%len(circuit)]],x) for i,x in enumerate(circuit))}
                    circuits = {tuple(x for i,x in enumerate(circuit) if not Polyhedron.colinear([circuit[i-1],x,circuit[(i+1)%len(circuit)]])) for circuit in circuits}
                    circuits = {circuit for circuit in circuits if len(circuit)}
                    print(circuits)
                    if Polyhedron.find_exterior_circuit(circuits) is None:
                        for circuit in circuits:
                            face = set()
                            for i, p2 in enumerate(circuit):
                                face.add(frozenset([circuit[i-1],p2]))
                            new_faces.append(frozenset(face))
                    else:
                        new_faces.append(face)
        for face_index, face in enumerate(new_faces):
            updated = True
            while updated:
                updated = False
                for edge1 in new_faces[face_index]:
                    for edge2 in new_faces[face_index]:
                        if edge1 == edge2:
                            continue
                        intersect = Polyhedron.intersect_segments(edge1, edge2)
                        if intersect is not None and not any(distance(intersect,point) < 0.001 for point in edge1|edge2):
                            print(edge1, edge2, intersect)
                            new_faces[face_index] -= frozenset([edge1, edge2])
                            p1, p2 = edge1 
                            p3, p4 = edge2
                            new_faces[face_index] |= frozenset([frozenset([p1, intersect]),frozenset([p2,intersect]),frozenset([p3,intersect]),frozenset([p4,intersect])])
                            updated = True
                            break
                    if updated:
                        break
        updated = True
        while updated:
            updated = False
            for face_index1, face1 in enumerate(new_faces):
                for face_index2, face2 in enumerate(new_faces):
                    if face_index1 < face_index2 and Polyhedron.coplanar({point for edge in new_faces[face_index1] for point in edge}|{point for edge in new_faces[face_index2] for point in edge}):
                        circuits = Polyhedron.circuit_helper(face1|face2)
                        circuits = {tuple(x for i,x in enumerate(circuit) if not Polyhedron.colinear([circuit[i-1],x,circuit[(i+1)%len(circuit)]])) for circuit in circuits}
                        circuits = {circuit for circuit in circuits if len(circuit)}
                        if Polyhedron.find_exterior_circuit(circuits) is not None:
                            del new_faces[face_index2]
                            del new_faces[face_index1]
                            new_faces.append(face1|face2)
                            updated = True
                            break
                if updated:
                    break
        face_index = 0
        while face_index < len(new_faces):
            circuits = Polyhedron.circuit_helper(new_faces[face_index])
            circuits = {tuple(x for x in circuit if not any(Polyhedron.colinear([circuit[i-1], x, circuit[(i+1)%len(circuit)]]) for i,x in enumerate(circuit))) for circuit in circuits}
            circuits = {circuit for circuit in circuits if len(circuit)}
            print(circuits)
            if Polyhedron.find_exterior_circuit(circuits) is None:
                del new_faces[face_index]
                for circuit in circuits:
                    face = set()
                    for i, p2 in enumerate(circuit):
                        p1 = circuit[i-1]
                        face.add(frozenset([p1,p2]))
                    new_faces.append(frozenset(face))
            else:
                face_index += 1
        print("New Faces:")
        for face in new_faces:
            print(face)
        poly = Polyhedron()
        poly.verts = list({point for face in new_faces for edge in face for point in edge})
        poly.edges = list({frozenset(poly.verts.index(point) for point in edge) for face in new_faces for edge in face})
        poly.faces = [frozenset(poly.edges.index(frozenset(poly.verts.index(point) for point in edge)) for edge in face) for face in new_faces if len(face)]
        return poly
    def subtract(self, other, use_cpp=True):
        new_poly = Polyhedron()
        new_poly.verts = list(self.verts)
        new_poly.edges = list(self.edges)
        new_poly.faces = list(self.faces)
        poly = self.intersect(other, use_cpp)
        print(len(poly.verts), len(poly.edges), len(poly.faces))
        return Polyhedron.add_subtract_helper(poly,new_poly)
    def add(self, other, use_cpp=True):
        new_poly1 = Polyhedron()
        new_poly1.verts = list(other.verts)
        new_poly1.edges = list(other.edges)
        new_poly1.faces = list(other.faces)
        new_poly2 = Polyhedron()
        new_poly2.verts = list(self.verts)
        new_poly2.edges = list(self.edges)
        new_poly2.faces = list(self.faces)
        poly = self.intersect(other, use_cpp)
        new_poly1 = new_poly1.subtract(poly, use_cpp)
        return Polyhedron.add_subtract_helper(new_poly1, new_poly2)

                        
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

if __name__=="__main__":
    #poly1 = get_cube(factors=2, angles=(pi/4,pi/4,pi/4))
    #poly2 = get_cube((1,1,1), factors=2)
    #poly1 = get_cube(factors=2, angles=(pi/3,pi/3,pi/3))
    #poly2 = get_cube((1,1,1), factors=2)
    #poly1.dump('poly1.ply')
    #poly2.dump('poly2.ply')
    #poly1 = get_cube(factors=2, angles=(pi/6,pi/6,pi/6))
    #poly2 = get_cube((1,1,1), factors=2)
    poly1 = get_cube(factors=2, angles=(pi/5,pi/5,pi/5))
    poly2 = get_cube((1,1,1), factors=2)
    #poly1 = get_cube()
    #poly2 = get_cube((1.5, 1.5, 0), factors=3)
    #poly2 = get_cube(factors=(1,1,5), angles=(pi/10,pi/10,pi/10))
    #poly2 = get_cube((0,0,0), factors=3)
    #poly1 = get_cube(factors=(1,1,5))
    #poly2 = get_cube((1,0,0), factors=3)
    #poly2 = get_cube(factors=5)
    #poly2 = get_cube((1,1,2), factors=2)
    #poly2 = get_cube(factors=2)
    #poly1 = get_cube(factors=5)
    #poly1 = get_cube()
    #poly1 = poly1.add(get_cube((1,0,0)))
    #poly2 = get_cube((1,0,0),(1,3,1))
    #poly1 = poly1.add(get_cube((0,1,0)))
    #poly1 = get_cube((0,0,0),(3,3,1))
    #poly1 = poly1.subtract(get_cube((1.5,0,0),(2,1,1)))
    #poly1 = get_cube((0,0,0), (2,3,1))
    #poly1 = poly1.subtract(get_cube((0.5,0,0),(1,.33,1)))
    #poly1 = get_cube(factors=(2,2,1))
    #poly2 = get_cube((1,1,0))
    #poly3 = poly1.subtract(get_cube((0,0,1),factors=4))
    start_time = time.perf_counter()
    #poly3 = poly1.intersect(poly2, use_cpp=True)
    print("Time:", time.perf_counter()-start_time)
    poly3 = poly1.intersect(poly2, use_cpp=True)
    print([vert for vert in poly1.verts])
    print([vert for vert in poly2.verts])

    
    #intersections = poly1.intersect(poly2, True)
    #poly3 = poly1.add(poly2, True)
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
    meters = {1}
    print(poly1.faces)
    print(poly2.faces)
    print(len(poly1.edges))
    #print(poly3.faces)
    colors = {face:tuple(random.randrange(256) for i in range(3)) for face in poly3.faces}
    count = 0
    temp = 0
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
                if event.key == pygame.K_z:
                        camera.zoom *= 1.1
                if event.key == pygame.K_x:
                    camera.zoom /= 1.1
        screen.fill("black")
        for edge in poly1.edges:
            p1, p2 = tuple(poly1.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "red", p1, p2)
        for edge in poly2.edges:
            p1, p2 = tuple(poly2.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "blue", p1, p2)
        for edge in poly3.edges:
            p1, p2 = tuple(poly3.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "white", p1, p2)
        if count % 1000 == 0:
            temp += 1
            temp %= len(poly3.faces)
        for face_index, face in enumerate(poly3.faces):
            if face_index != temp:
                continue
            color = colors[face]
            for edge_index in face:
                edge = poly3.edges[edge_index]
                p1, p2 = tuple(poly3.verts[index] for index in edge)
                p1, p2 = camera.project(p1), camera.project(p2)
                if p1 and p2:
                    p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                    p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                    #print(p1,p2)
                    pygame.draw.line(screen, "green", p1, p2)
        '''
        for poly in intersections:
            for edge in poly.edges:
                p1, p2 = [poly.verts[index] for index in edge]
                p1, p2 = camera.project(p1), camera.project(p2)
                if p1 and p2:
                    p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                    p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                    #print(p1,p2)
                    pygame.draw.line(screen, "white", p1, p2)
        '''
        '''
        for edge in hull_edges:
            p1, p2 = edge
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "green", p1, p2)
        '''
        pygame.display.flip()
        count += 1
        count %= 1000
