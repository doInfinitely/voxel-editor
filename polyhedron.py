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
                if Polyhedron.coplanar(set([round_point(x) for x in path])|{round_point(self.verts[x]) for x in y}):
                    current = self.edges.index(y)
                    circuits.update(self.circuits(face_index, start, previous, current, path, old_circuits|circuits))
        circuits_list = list(circuits)
        for i,x in enumerate(circuits_list):
            for j,y in enumerate(circuits_list[i+1:]):
                if not len(set(x)-set(y)):
                    y_r = tuple(reversed(y))
                    if (y[y.index(x[0]):]+y[:y.index(x[0])])[:len(x)] == x or (y_r[y_r.index(x[0]):]+y_r[:y_r.index(x[0])])[:len(x)] == x:
                        if y in circuits:
                            circuits.remove(y)
        return circuits
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
        if Polyhedron.colinear(points):
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

    @timeout(1)
    def timed_solve(*args, **kwargs):
        return solve(*args, **kwargs)
    def inside_triangle(triangle, point):
        if point in triangle:
            return True
        triangle = list(triangle)
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
                    print('triangle intersect', tuple(alpha*p2[i]+(1-alpha)*p1[i] for i in range(3)))
                    return (alpha,tuple(alpha*p2[i]+(1-alpha)*p1[i] for i in range(3)))
            except TypeError:
                pass
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
        return tuple((y[0],y[1],0) for y in output)
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
                    remainder = tuple([x for k,x in enumerate(remainder) if x != remainder[k-1]])
                    return tuple([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]]), remainder
        sys.exit(1)
    def triangulate(circuit):
        #print(circuit)
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
                            output.append((alpha,tuple(alpha*p2[i]+(1-alpha)*p1[i] for i in range(3)),(p3,p4),circuit))
                    except TypeError:
                        pass
        return output
    def face_intersect(self, segment):
        output = []
        for face_index, face in enumerate(self.faces):
            for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(self.circuits(face_index))):
                intersect = Polyhedron.intersect_triangle(segment, triangle)
                if intersect:
                    print('inter', segment, triangle, intersect)
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
            output[-1] = output[-1][:i] + tuple([first_exterior_intersection[1], last_interior_intersection[1]]) + last_interior_intersection[3][j:] + last_interior_intersection[3][:j] + tuple([last_interior_intersection[1], first_exterior_intersection[1]]) + output[-1][i:]
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
        return output[0]
    def is_inface(self, point):
        return any(Polyhedron.inside_triangle(x,point) for face_index,face in enumerate(self.faces) for x in Polyhedron.triangulate(Polyhedron.circuit_cut(self.circuits(face_index))))
    def is_inside(self, point):
        #print('is_inside', point)
        if self.is_inface(point):
            return True
        vec = (random.random()*2-1,random.random()*2-1,random.random()*2-1)
        output = []
        for face_index,face in enumerate(self.faces):
            circuit = Polyhedron.circuit_cut(self.circuits(face_index))
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
    def loads(string):
        pass
    def load(filename):
        with open(filename, 'r') as f:
            return Polyhedron.loads(f.read())
    def circuit_overlap(circuits1, circuits2):
        circuit1 = Polyhedron.circuit_cut(circuits1)
        circuit2 = Polyhedron.circuit_cut(circuits2)
        if any(Polyhedron.intersect_segments(frozenset([circuit1[i-1],x]),frozenset([circuit2[j-1],y])) for i,x in enumerate(circuit1) for j,y in enumerate(circuit2)):
            return True
        if any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit1) for y in circuit2):
            return True
        if any(Polyhedron.inside_triangle(x,y) for x in Polyhedron.triangulate(circuit2) for y in circuit1):
            return True
        return False
    def project_on_face_plane(self, face_index, point):
        circuit = Polyhedron.circuit_cut(self.circuits(face_index))
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
    def intersect(self, other):
        polys = [self, other]
        old_edges = set()
        new_edges = set()
        new_edges2 = set()
        def path_along_faces(p1, p2, face_path): 
            ray = other_poly.project_ray_on_face_plane(face_path[-1], (p1,p2))
            #print('ray', ray, (p1,p2), round_float((distance(ray[1],p2))), face_path)
            #if round_float(distance(ray[0],p1)) == 0 and round_float(distance(ray[1],p2)) == 0:
            if round_float(distance(ray[0],p1)) == 0 and round_float(distance(ray[1],p2)) == 0:
                #print('hey', [p1,p2])
                return [([p1,p2], face_path)]
            if ray[0] == ray[1]:
                faces = [face_index for face_index, face in enumerate(poly.faces) if p1 in {poly.verts[index] for edge_index in poly.faces[face_index] for index in poly.edges[edge_index]}]
                output = []
                for face_index in faces:
                    intersections, edges = Polyhedron.intersect_circuits(Polyhedron.circuit_cut(poly.circuits(face_index)), Polyhedron.circuit_cut(other_poly.circuits(face_path[-1])))
                    for edge in edges:
                        if p1 in edge:
                            p = list(edge-frozenset([p1]))[0]
                            for face_index2, face2 in enumerate(other_poly.faces):
                                if face_index2 not in face_path and any(Polyhedron.inside_triangle(triangle, p) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(other_poly.circuits(face_index2)))):
                                     output.extend([([p1]+x[0],x[1]) for x in path_along_faces(p, p2, face_path+[face_index_2])])
                return output
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
                        print(face_index2,  any(Polyhedron.inside_triangle(triangle, p) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(other_poly.circuits(face_index2)))))
                        if face_index2 not in face_path and any(Polyhedron.inside_triangle(triangle, p) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(other_poly.circuits(face_index2)))):
                            print(path_along_faces(p, p2, face_path+[face_index2]))
                            output.extend([([p1]+x[0],x[1]) for x in path_along_faces(p, p2, face_path+[face_index2])])
                return output
            return [([],[])]
        for poly_index, poly in enumerate(polys):
            other_poly = polys[poly_index-1]
            seen_verts = set()
            seen_leaves = set()
            for vert_index, vert in enumerate(poly.verts):
                if vert_index in seen_verts:
                    continue
                if other_poly.is_inface(vert):
                    print('hey')
                    continue
                leaves = []
                nonleaves = []
                root_in_poly = other_poly.is_inside(vert)
                q = [vert_index]
                face_lookup = dict()
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
                        if len(intersects):
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
                                old_edges.add(frozenset([poly.verts[v_i1], intersects[0][1]]))
                            else:
                                for intersect in intersects[1:]:
                                    if intersect[0] > intersects[0][0]:
                                        old_edges.add(frozenset([intersects[0][1], intersect[1]]))
                                        break
                        else:
                            q.append(v_i2)
                            if root_in_poly:
                                old_edges.add(frozenset([poly.verts[v_i1], poly.verts[v_i2]]))
                print('leaves',leaves)
                print(face_lookup)
                #print(vert, nonleaves, leaves)
                if len(leaves) > 2 and frozenset([round_point(x) for x in leaves]) not in seen_leaves:
                    seen_leaves.add(frozenset([round_point(x) for x in leaves]))
                    leaf_centroid = tuple(sum(x[i] for x in leaves)/len(leaves) for i in range(3))
                    vec0 = tuple(leaf_centroid[i]-vert[i] for i in range(3))
                    gamma = symbols("gamma")
                    eqs = [Eq(vec0[0]*vec0[1]+vec0[1]*-vec0[0]+gamma*vec0[2],0)]
                    solutions = solve(eqs, dict=True)
                    vec1 = (vec0[1], -vec0[0], solutions[0][gamma])
                    vec2 = cross3D(vec0, vec1)
                    vec1 = tuple(vec1[i]/distance(vec1) for i in range(3))
                    vec2 = tuple(vec2[i]/distance(vec2) for i in range(3))
                    #print(vec0, vec1, vec2)
                    leaf_projection = []
                    for leaf in leaves:
                        alpha, beta, gamma = symbols("alpha beta gamma")
                        eqs = [Eq(alpha*vec0[i]+leaf[i],beta*vec1[i]+gamma*vec2[i]) for i in range(3)]
                        solutions = solve(eqs, dict=True)
                        alpha = solutions[0][alpha]
                        leaf_projection.append(tuple(alpha*vec0[i]+leaf[i] for i in range(3)))
                    planar_leaf_projection = Polyhedron.make_planar(leaf_projection)
                    reverse_projection = dict(zip(planar_leaf_projection, leaves))
                    point = min(planar_leaf_projection)
                    path = []
                    rotated_mapping = {x:x for x in planar_leaf_projection}
                    reverse_mapping = {x:x for x in planar_leaf_projection}
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
                    print('path',path)
                    path = [reverse_projection[x] for x in path]
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
                                        #print('blah',[x,point_path[i+1]], sum([abs(x[j]-point_path[i+1][j]) for j in range(3)]))
                                        new_edges.add(frozenset([x,point_path[i+1]]))
                            if not len(paths):
                                new_edges.add(frozenset([p1,p2]))
            print('leaves', seen_leaves)
        print(len(new_edges|old_edges))        
        '''
        components = {frozenset([edge1,edge2]) for edge1 in new_edges|old_edges for edge2 in new_edges|old_edges if len(edge1 & edge2) == 1 and Polyhedron.coplanar(edge1|edge2)}
        components = [set(x) for x in components]
        update = True
        while update:
            print([len(x) for x in components])
            update = False
            for edge in new_edges|old_edges:
                for index1, component1 in enumerate(list(components)):
                    points1 = {point for point in edge}
                    points2 = {point for edge1 in component1 for point in edge1}
                    if len(points1 & points2) == 1 and Polyhedron.coplanar(points1|points2) and edge not in component1:
                        components[index1].add(edge)
                        update = True
        for component in components:
            #print(component)
            endpoints = set()
            for edge in component:
                for point in edge:
                    if point not in endpoints:
                        endpoints.add(point)
                    else:
                        endpoints.remove(point)
            try:
                p1, p2 = endpoints
                new_edges2.add(frozenset([p1,p2]))
            except ValueError:
                pass
            starting_faces = []
            for face_index, face in enumerate(other_poly.faces):
                if any(Polyhedron.inside_triangle(triangle, p1) for triangle in Polyhedron.triangulate(Polyhedron.circuit_cut(other_poly.circuits(face_index)))):
                    starting_faces.append(face_index)
            print(p1, "starting_faces",starting_faces)
            for face_index1 in starting_faces:
                for point_path, face_path in paths(p1,p2,[face_index1]):
                    if point_path is None:
                        new_edges.add(frozenset([p1,p2]))
                        break
                    #print(point_path, face_path)
                    if len(point_path):
                        for i,x in enumerate(point_path[:-1]):
                            #print('blah',[x,point_path[i+1]], sum([abs(x[j]-point_path[i+1][j]) for j in range(3)]))
                            new_edges.add(frozenset([x,point_path[i+1]]))
        '''
        #print('old',old_edges)
        #print('new',new_edges)
        points = {point for edge in new_edges|old_edges for point in edge}
        points = np.array(list(points))
        hull = ConvexHull(points)
        print(hull.vertices)
        print(hull.simplices)
        print(points)
        simplices = [Polyhedron.make_clockwise([tuple(round_point(tuple(points[index])) for index in simplex)])[0] for simplex in hull.simplices]
        hull_edges = set()
        update = True
        while update:
            print([len(x) for x in simplices])
            update = False
            for simplex1 in list(simplices):
                for simplex2 in list(simplices):
                    if simplex1 != simplex2:
                        print(simplex1, simplex2, Polyhedron.coplanar(set(simplex1+simplex2)))
                    if simplex1 != simplex2 and Polyhedron.coplanar(set(simplex1+simplex2)) and len(set(simplex1)|set(simplex2)) > 1:
                        simplex3 = []
                        for i,x in enumerate(simplex1):
                            if simplex1[i-1] in simplex2 and x in simplex2:
                                index1 = simplex2.index(simplex1[i-1])
                                index2 = simplex2.index(x)
                                print(index1, index2)
                                simplex3.extend(tuple(y for j,y in enumerate(simplex2) if j not in [index1, index2]) + tuple([x]))
                            else:
                                simplex3.append(x)    
                        simplices = [simplex for simplex in simplices if simplex != simplex1 and simplex != simplex2] + [tuple(simplex3)]
                        update = True
                        break
                if update:
                    break
        for simplex in simplices:
            for i, x in enumerate(simplex):
                edge = frozenset([simplex[i-1], x])
                if edge in new_edges|old_edges:
                    hull_edges.add(edge)
        return new_edges|old_edges, hull_edges

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
    poly1 = get_cube(factors=2, angles=(pi/4,pi/4,pi/4))
    poly2 = get_cube((1,1,1), factors=2)

    poly1.dump('poly1.ply')
    poly2.dump('poly2.ply')
    #poly1 = get_cube(factors=2, angles=(pi/6,pi/6,pi/6))
    #poly2 = get_cube((1,1,1), factors=2)
    #poly1 = get_cube()
    #poly2 = get_cube((1.5, 1.5, 0), factors=3)
    #poly1 = get_cube(factors=(1,1,5))
    #poly2 = get_cube((0,0,0), factors=3)
    #poly1 = get_cube(factors=5)
    #poly2 = get_cube((1,1,2), factors=2)
    print([vert for vert in poly1.verts])
    print([vert for vert in poly2.verts])

    print(poly1.circuits(0))
    edges, hull_edges = poly1.intersect(poly2)
    for edge in edges:
        print(edge)
    #edges2 = poly2.union(poly1)
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
        for edge in edges:
            p1, p2 = edge
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "white", p1, p2)
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
