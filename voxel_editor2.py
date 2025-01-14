from sympy import nsolve, solve, symbols, Eq, diff
import sympy
from math import sin, cos, pi, sqrt, copysign
import numpy as np
import subprocess
import random
import time

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

class Polyhedron:
    def __init__(self):
        self.verts = []
        self.edges = []
        self.faces = []
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
    def circuits(self, face_index):
        face = self.faces[face_index]
        if not len(face):
            return []
        edge_lookup = dict()
        for edge_index in face:
            edge = self.edges[edge_index]
            for index in self.edges[edge_index]:
                point = self.verts[index]
                #print(point, edge)
                if point not in edge_lookup:
                    edge_lookup[point] = set()
                edge_lookup[point].add(edge)
        circuits = set()
        edges = {self.edges[x] for x in self.faces[face_index]}
        #print(self.faces[face_index])
        while len(edges):
            start = list(edges)[0]
            edges.remove(start)
            previous = start
            point = self.verts[list(start)[0]]
            circuit = []
            circuit.append(point)
            #print([[self.verts[index] for index in self.edges[edge_index]] for edge_index in face])
            #print(point, edge_lookup[point])
            current = list(edge_lookup[point] - set([start]))[0]
            #print([[self.verts[index] for index in edge] for edge in self.edges])
            while current != start:
                edges.remove(current)
                point = self.verts[list(current - previous)[0]]
                circuit.append(point)
                previous = current
                #print(point, [[self.verts[index] for index in edge] for edge in edge_lookup[point]])
                #print('circuit', circuit)
                for p in edge_lookup:
                    if current in edge_lookup[p]:
                        edge_lookup[p].remove(current)
                current = list(edge_lookup[point])[0]
            circuits.add(tuple(circuit))
        return circuits
    def intersect_segments(edge1, edge2):
        output = []
        p1, p2 = edge1
        p3, p4 = edge2
        alpha, beta = symbols("alpha beta")
        eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], beta*p3[i]+(1-beta)*p4[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            alpha, beta = solutions[0][alpha], solutions[0][beta]
            if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1:
                output.append(tuple(alpha*p1[i]+(1-alpha)*p2[i] for i in range(3)))
        return output
    # returns any edge that intersects the line segment between p1 and p2
    def intersect(self, p1, p2):
        output = []
        for edge in self.edges:
            if len(Polyhedron.intersect_segments(frozenset([p1,p2]), tuple(self.verts[index] for index in edge))):
                output.append(edge)
        return output
    def inside_triangle(triangle, point):
        triangle = list(triangle)
        alpha, beta = symbols("alpha beta")
        #print(triangle)
        exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
        exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
        eqs = [Eq(x,point[i]) for i,x in enumerate(exprs)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            alpha, beta = solutions[0][alpha], solutions[0][beta]
            if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1:
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
        return [x==signs[i-1] and x == signs[(i+1)%len(signs)] for i,x in enumerate(signs)]

    # A circuit is a representation of a Jordan Polygon
    def clip_ear(circuit):
        is_convex = Polyhedron.convex_angles(circuit)
        for i in range(len(circuit)):
            if is_convex[i]:
                for j,y in enumerate(circuit):
                    if j != i and j != (i+1)%len(circuit) and j != (i-1)%len(circuit):
                        if Polyhedron.inside_triangle([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]],y):
                            break
                else:
                    return tuple([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]]), tuple([x for k,x in enumerate(circuit) if k != i])
    def triangulate(circuit):
        output = []
        remainder = circuit
        while len(remainder) > 3:
            ear, remainder = Polyhedron.clip_ear(circuit)
            output.append(ear)
        output.append(remainder)
        return output

    # Given a set of circuits where one circuit surrounds the others, find that circuit
    def find_exterior_circuit(circuits):
        for i,x in enumerate(circuits):
            triangulation = Polyhedron.triangulate(x) 
            for j,y in enumerate(circuits):
                if j != i:
                    doublebreak = False
                    for point in y:
                        if not any(Polyhedron.inside_triangle(z,point) for z in triangulation):
                            doublebreak = True
                            break
                    if doublebreak:
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
                    alpha, beta = solutions[0][alpha], solutions[0][beta]
                    if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1:
                        output.append((alpha,tuple(alpha*p2[i]+(1-alpha)*p1[i] for i in range(3)),(p3,p4),circuit))
        return output
    # Converts a polygon with one or more holes into a single circuit by forming cuts from the holes to the exterior
    def circuit_cut(circuits):
        exterior = Polyhedron.find_exterior_circuit(circuits)
        if exterior is None:
            return None
        output = [x for x in circuits if x != exterior]
        output.append(exterior)
        while len(output) > 1:
            interior = output[0]
            segment = (interior[0],output[-1][0])
            intersections = Polyhedron.circuit_intersect(segment, output)
            intersections.sort()
            for intersection in intersections:
                if intersection[3] == output[-1]:
                    first_exterior_intersection = intersection
                    break
            for intersection in intersections:
                if intersection[3] == output[-1]:
                    continue
                if intersection[0] > first_exterior_intersection[0]:
                    break
                last_interior_intersection = intersection
            for i,x in enumerate(output[-1]):
                if x in first_exterior_intersection[2] and output[-1][i-1] in first_exterior_intersection[2]:
                    break
            for j,y in enumerate(last_interior_intersection[3]):
                if y in last_interior_intersection[2] and last_interior_intersection[3][j-1] in first_exterior_intersection[2]:
                    break
            output[-1] = output[-1][:i] + tuple([first_exterior_intersection[1], last_interior_intersection[1]]) + tuple(reversed(last_interior_intersection[3][:j])) + tuple(reversed(last_interior_intersection[3][j:])) + tuple([last_interior_intersection[1], first_exterior_intersection[1]]) + output[-1][i:]
            output.remove(last_interior_intersection[3])
            output[-1] = list(output[-1])
            i = 1
            while i < len(output[-1]):
                if output[-1][i] == output[-1][i-1]:
                    del output[-1][i]
                else:
                    i += 1
            output[-1] = tuple(output[-1])
        return output[0]

    def write_poly_file(self, filename='polyhedron.poly'):
        with open(filename, 'w') as f:
            f.write('# part 1 - node list\n')
            f.write('{0} 3 0 0\n\n'.format(len(self.verts)))
            for i,x in enumerate(self.verts):
                f.write('{0} {1} {2}, {3}\n'.format(i, x[0], x[1], x[2]))
            f.write('\n# part 2 - facet list\n')
            f.write('{0} 0\n\n'.format(len(self.faces)))
            hole_index = 0
            for i,x in enumerate(self.faces):
                holes = []
                circuits = self.circuits(i)
                exterior = Polyhedron.find_exterior_circuit(circuits)
                f.write('{0} {1}\n'.format(len(circuits), len(circuits)-1))
                for circuit in circuits:
                    line = [len(circuit)] + [self.verts.index(y) for y in circuit]
                    f.write(' '.join([str(y) for y in line])+'\n')
                    if circuit != exterior:
                        triangles = Polyhedron.triangulate(circuit)
                        holes.append(tuple(0.5*(triangles[0][0][i]+0.5*(triangles[0][1][i]+triangles[0][2][i])) for i in range(3)))
                for hole in holes:
                    f.write("{0} {1} {2} {3}\n".format(hole_index, hole[0], hole[1], hole[2]))
                    hole_index += 1
                f.write('\n')
    def tetrahedralize(self):
        self.write_poly_file()
        out = subprocess.run(["./tetgen1.6.0/tetgen", "-p", "polyhedron.poly"], capture_output=True)
        tetrahedrons = set()
        with open('polyhedron.1.ele') as f:
            first = True
            for line in f:
                if first:
                    first = False
                    continue
                if line.startswith('#'):
                    continue
                splitline = line.split()
                tetrahedrons.add(frozenset([int(x) for x in splitline[1:5]]))
        points = dict()
        with open('polyhedron.1.node') as f:
            first = True
            for line in f:
                if first:
                    first = False
                    continue
                if line.startswith('#'):
                    continue
                splitline = line.split()
                points[int(splitline[0])] = tuple([int(x) for x in splitline[1:4]])
        tetrahedrons = {frozenset([points[y] for y in x]) for x in tetrahedrons}
        '''
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for tetra in tetrahedrons:
            for p1 in tetra:
                for p2 in tetra:
                    if p1 != p2:
                        xline = np.linspace(p1[0], p2[0], 1000)
                        yline = np.linspace(p1[1], p2[1], 1000)
                        zline = np.linspace(p1[2], p2[2], 1000)
                        ax.plot3D(xline, yline, zline, 'gray')
        for edge in self.edges:
            p1, p2 = [self.verts[x] for x in edge]
            xline = np.linspace(p1[0], p2[0], 1000)
            yline = np.linspace(p1[1], p2[1], 1000)
            zline = np.linspace(p1[2], p2[2], 1000)
            ax.plot3D(xline, yline, zline, 'black')
        plt.show()
        '''
        return tetrahedrons
    def triangle_noncoplanar_intersect(triangle, segment):
        triangle = list(triangle)
        alpha, beta, gamma = symbols("alpha beta gamma")
        p1, p2 = segment
        #print(triangle)
        exprs = [alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3)]
        exprs = [beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(exprs)]
        eqs = [Eq(x,gamma*p1[i]+(1-gamma)*p2[i]) for i,x in enumerate(exprs)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            for x in solutions:
                try:
                    alpha, beta, gamma = x[alpha], x[beta], x[gamma]
                    break
                except KeyError:
                    pass
            else:
                raise KeyError
            if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1 and gamma >= 0 and gamma <= 1:
                return tuple(gamma*p1[i]+(1-gamma)*p2[i] for i in range(3))
        return None
    def interval_interesect(segment1, segment2):
        output = []
        p1, p2 = segment1
        p3, p4 = segment2
        alpha = symbols("alpha")
        eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], p3[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        alpha1 = None
        if len(solutions):
            alpha1 = solutions[0][alpha]
        alpha = symbols("alpha")
        eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], p4[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        alpha2 = None
        if len(solutions):
            alpha2 = solutions[0][alpha]
        if alpha1 is not None and alpha1 >= 0 and alpha1 <= 1:
            if alpha2 is not None and alpha2 >= 0 and alpha2 <= 1:
                return [p3, p4]
        alpha = symbols("alpha")
        eqs = [Eq(alpha*p3[i]+(1-alpha)*p4[i], p1[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        alpha3 = None
        if len(solutions):
            alpha3 = solutions[0][alpha]
        alpha = symbols("alpha")
        eqs = [Eq(alpha*p3[i]+(1-alpha)*p4[i], p2[i]) for i in range(3)]
        solutions = solve(eqs, dict=True)
        alpha4 = None
        if len(solutions):
            alpha4 = solutions[0][alpha]
        if alpha3 is not None and alpha3 >= 0 and alpha3 <= 1:
            if alpha4 is not None and alpha4 >= 0 and alpha4 <= 1:
                return [p1, p2]
        if alpha1 is not None and alpha1 >= 0 and alpha1 <= 1:
            if alpha3 is not None and alpha3 >= 0 and alpha3 <= 1:
                return [p3, p1]
            if alpha4 is not None and alpha4 >= 0 and alpha4 <= 1:
                return [p3, p2]
        if alpha2 is not None and alpha2 >= 0 and alpha2 <= 1:
            if alpha3 is not None and alpha3 >= 0 and alpha3 <= 1:
                return [p4, p1]
            if alpha4 is not None and alpha4 >= 0 and alpha4 <= 1:
                return [p4, p2]
        return []
    def tetrahedron_surface_intersect(tetrahedron, segment):
        start = time.time()
        output = set()
        for excluded_point in tetrahedron:
            triangle = {x for x in tetrahedron if x != excluded_point}
            if Polyhedron.coplanar(triangle|set(segment)):
                for excluded_point in triangle:
                    edge = {x for x in triangle if x != excluded_point}
                    if Polyhedron.colinear(edge|set(segment)):
                        intersections = Polyhedron.interval_interesect(edge, segment)
                        output.update(intersections)
                        break
                else:
                    intersections = [x[1] for x in Polyhedron.circuit_intersect(segment, [list(triangle)])]
                    output.update(intersections)
                    break
            else:
                print(triangle, segment)
                intersection = Polyhedron.triangle_noncoplanar_intersect(triangle, segment)
                if intersection is not None:
                    output.add(intersection)
        output = list(output)
        for i,x in enumerate(output):
            remove = False
            for j,y in enumerate(output):
                if i < j and x == y:
                    remove = True
                    break
            if remove:
                del output[j]
        print(time.time()-start)
        return set(output)
    def inside_tetrahedron(tetrahedron, point):
        t = list(tetrahedron)
        alpha, beta, gamma = symbols("alpha beta gamma")
        exprs = [alpha*t[0][i]+(1-alpha)*t[1][i] for i in range(3)]
        exprs = [beta*x+(1-beta)*t[2][i] for i,x in enumerate(exprs)]
        exprs = [gamma*x+(1-gamma)*t[3][i] for i,x in enumerate(exprs)]
        eqs = [Eq(x,point[i]) for i,x in enumerate(exprs)]
        solutions = solve(eqs, dict=True)
        if len(solutions):
            alpha, beta, gamma = solutions[0][alpha], solutions[0][beta], solutions[0][gamma]
            if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1 and gamma >= 0 and gamma <= 1:
                return True
        return False

    def tetrahedron_get_edges(tetrahedron):
        output = set()
        for i,point1 in enumerate(tetrahedron):
            for j,point2 in enumerate(tetrahedron[i+1:]):
                output.add(frozenset([point1,point2]))
        return output
    def tetrahedron_subtract(tetra1, tetra2):
        intersect_lookup1 = dict()
        intersect_lookup2 = dict()
        edges = []
        new_points = set()
        for edge in Polyhedron.tetrahedron_get_edges(tetra1):
            intersect_lookup1[edge] = list(Polyhedron.tetrahedron_surface_intersect(tetra2, edge))
        for edge in Polyhedron.tetrahedron_get_edges(tetra2):
            intersect_lookup2[edge] = list(Polyhedron.tetrahedron_surface_intersect(tetra1, edge))
        for edge in intersect_lookup1:
            print('blah', edge, intersect_lookup1[edge])
            point1, point2 = edge
            if len(intersect_lookup1[edge]) == 1:
                is_outside1 = not Polyhedron.inside_tetrahedron(tetra2, point1)
                is_outside2 = not Polyhedron.inside_tetrahedron(tetra2, point2)
                if is_outside1 and is_outside2:
                    edges.append(frozenset([point1,point2]))
                elif is_outside1:
                    edges.append(frozenset([point1,intersect_lookup1[edge][0]]))
                    new_points.add(intersect_lookup1[edge][0])
                elif is_outside2:
                    edges.append(frozenset([point2,intersect_lookup1[edge][0]]))
                    new_points.add(intersect_lookup1[edge][0])
            elif len(intersect_lookup1[edge]) == 2:
                if distance(point1, intersect_lookup1[edge][0]) + distance(point2, intersect_lookup1[edge][1]) < distance(point2, intersect_lookup1[edge][0]) + distance(point1, intersect_lookup1[edge][1]):
                    edges.append(frozenset([point1,intersect_lookup1[edge][0]]))
                    edges.append(frozenset([point2,intersect_lookup1[edge][1]]))
                else:
                    edges.append(frozenset([point1,intersect_lookup1[edge][1]]))
                    edges.append(frozenset([point2,intersect_lookup1[edge][0]]))
                new_points.add(intersect_lookup1[edge][0])
                new_points.add(intersect_lookup1[edge][1])
            else:
                edges.append(frozenset([point1,point2]))
        for edge in intersect_lookup2:
            print('bleh', edge, intersect_lookup2[edge])
            point1, point2 = edge
            is_inside1 = Polyhedron.inside_tetrahedron(tetra1, point1)
            is_inside2 = Polyhedron.inside_tetrahedron(tetra1, point2)
            if len(intersect_lookup2[edge]) == 1:
                if is_inside1 and is_inside2:
                    edges.append(frozenset([point1,point2]))
                elif is_inside1:
                    edges.append(frozenset([point1,intersect_lookup2[edge][0]]))
                    new_points.add(intersect_lookup2[edge][0])
                elif is_inside2:
                    edges.append(frozenset([point2,intersect_lookup2[edge][0]]))
                    new_points.add(intersect_lookup2[edge][0])
            elif len(intersect_lookup2[edge]) == 2:
                edges.append(frozenset(intersect_lookup2[edge])-frozenset(tetra1))
                new_points.add(intersect_lookup2[edge][0])
                new_points.add(intersect_lookup2[edge][1])
            else:
                if is_inside1 and is_inside2:
                    edges.append(frozenset([point1,point2]))
        for point1 in new_points:
            for point2 in new_points:
                if point1 < point2:
                    if frozenset([point1,point2]) not in edges:
                        for edge in Polyhedron.tetrahedron_get_edges(tetra1):
                            if Polyhedron.colinear(edge|frozenset([point1,point2])):
                                break
                        else:
                            edges.append(frozenset([point1,point2]))
        for edge1 in edges:
            double_break = False
            for edge2 in edges:
                if edge1 != edge2 and len(Polyhedron.intersect_segments(edge1, edge2)):
                    double_break = True
                    break
            if double_break:
                break
        if double_break:
            edges.remove(edge1)
            edges.remove(edge2)

        print([x for x in edges if len(x) == 2])
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for edge in edges:
            try:
                p1, p2 = edge
            except ValueError:
                print('edge', edge)
                continue
            print(p1,p2)
            xline = np.linspace(float(p1[0]), float(p2[0]), 1000)
            yline = np.linspace(float(p1[1]), float(p2[1]), 1000)
            zline = np.linspace(float(p1[2]), float(p2[2]), 1000)
            ax.plot3D(xline, yline, zline, 'black')
        plt.show()
        
    '''
    def subtract(self, other):
        self_tetra = self.tetrahedralize()
        other_tetra = other.tetrahedralize()
        for tetra1 in self_tetra:
            for tetra2 in other_tetra:
                for p1 in tetra2:
                    for p2 in tetra2:
                        if p1 != p2:
    '''                     

def polygon_with_hole():
    points = [(0,0,0),(3,0,0),(3,3,0),(0,3,0)]
    points.extend([(1,1,0),(2,1,0),(2,2,0),(1,2,0)])
    edges = [(i,(i+1)%4) for i in range(4)]
    edges.extend([(4+i,4+(i+1)%4) for i in range(4)])
    faces = [{i for i,x in enumerate(edges)}]
    polyhedron = Polyhedron()
    polyhedron.verts = points
    polyhedron.edges = [frozenset(x) for x in edges]
    polyhedron.faces = faces
    return polyhedron

polygon = polygon_with_hole()
circuits = polygon.circuits(0)
print(circuits)
print(Polyhedron.find_exterior_circuit(circuits))
print(Polyhedron.circuit_cut(circuits))
    
def cube_with_hole():
    points = [(0,0,0),(3,0,0),(3,3,0),(0,3,0)]
    points.extend([(1,1,0),(2,1,0),(2,2,0),(1,2,0)])
    points.extend([(0,0,3),(3,0,3),(3,3,3),(0,3,3)])
    points.extend([(1,1,3),(2,1,3),(2,2,3),(1,2,3)])
    edges = [(i,(i+1)%4) for i in range(4)]
    edges.extend([(4+i,4+(i+1)%4) for i in range(4)])
    edges.extend([(8+i,8+(i+1)%4) for i in range(4)])
    edges.extend([(12+i,12+(i+1)%4) for i in range(4)])
    edges.extend([(i,8+i) for i in range(4)])
    edges.extend([(4+i,12+i) for i in range(4)])
    faces = [{i for i in range(4)}|{4+i for i in range(4)},{8+i for i in range(4)}|{12+i for i in range(4)}]
    faces.extend([{0,16,17,8},{1,17,18,9},{2,18,19,10},{3,19,16,11}])
    faces.extend([{4,20,21,12},{5,21,22,13},{6,22,23,14},{7,23,20,15}])
    polyhedron = Polyhedron()
    polyhedron.verts = points
    polyhedron.edges = [frozenset(x) for x in edges]
    polyhedron.faces = faces
    return polyhedron

cube_with_hole = cube_with_hole()
#cube_with_hole.write_poly_file()
cube_with_hole.tetrahedralize()

def tetrahedron():
    points = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    faces = [{0,1,3},{1,2,5},{0,2,4},{3,4,5}]
    polyhedron = Polyhedron()
    polyhedron.verts = points
    polyhedron.edges = [frozenset(x) for x in edges]
    polyhedron.faces = faces
    return polyhedron

tetrahedron = tetrahedron()
segment = [(0.25,0.25,-1),(0.25,0.25,2)]
print(Polyhedron.tetrahedron_surface_intersect(tetrahedron.verts, segment))

#tetra1 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
#tetra2 = [(0,-.5,0),(1,-.5,0),(0,.5,0),(0,-.5,1)]

#tetra1 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
#tetra2 = [(0.5,-.5,0.25),(1,-.5,0.25),(0.5,.5,0.25),(0.5,-.5,1)]

#tetra1 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
#tetra2 = [(0.25,-.5,0.25),(1,-.5,0.25),(0.25,.5,0.25),(0.5,-.5,1)]

tetra1 = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
#tetra2 = [(0,0,0+1.5),(1,0,0+1.5),(0,1,0+1.5),(0,0,-1+1.5)]
tetra2 = []
for j in range(2):
    point = tuple(random.random() for i in range(3))
    while not Polyhedron.inside_tetrahedron(tetra1, point):
        point = tuple(random.random() for i in range(3))
    tetra2.append(point)
for j in range(2):
    point = tuple(random.random() for i in range(3))
    while Polyhedron.inside_tetrahedron(tetra1, point):
        point = tuple(random.random() for i in range(3))
    tetra2.append(point)

Polyhedron.tetrahedron_subtract(tetra1,tetra2)
