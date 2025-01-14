from sympy import nsolve, solve, symbols, Eq, diff
import sympy
from math import sin, cos, pi, sqrt, copysign
import numpy as np
from decimal import Decimal

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
            if not np.allclose(np.dot(a.astype('float'), x.astype('float')), b.astype('float'), rtol=0.1):
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
        output = []
        p1, p2 = edge1
        p3, p4 = edge2
        alpha, beta = symbols("alpha beta")
        eqs = [Eq(alpha*p1[i]+(1-alpha)*p2[i], beta*p3[i]+(1-beta)*p4[i]) for i in range(3)]
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
                    output.append(tuple(alpha*p1[i]+(1-alpha)*p2[i] for i in range(3)))
            except TypeError:
                    alpha = 0
                    output.append(tuple(alpha*p1[i]+(1-alpha)*p2[i] for i in range(3)))
        return output
    def intersect(self, edge):
        p1, p2 = edge
        circuits = []
        circuits.append([(self.x_min,self.y_min,self.z_min),(self.x_min,self.y_max,self.z_min),(self.x_max,self.y_max,self.z_min),(self.x_max,self.y_min,self.z_min)])
        circuits.append([(self.x_min,self.y_min,self.z_max),(self.x_min,self.y_max,self.z_max),(self.x_max,self.y_max,self.z_max),(self.x_max,self.y_min,self.z_max)])
        circuits.append([(self.x_min,self.y_min,self.z_min),(self.x_min,self.y_min,self.z_max),(self.x_max,self.y_min,self.z_max),(self.x_max,self.y_min,self.z_min)])
        circuits.append([(self.x_min,self.y_min,self.z_max),(self.x_min,self.y_min,self.z_max),(self.x_max,self.y_min,self.z_max),(self.x_max,self.y_min,self.z_min)])
        circuits.append([(self.x_min,self.y_min,self.z_min),(self.x_min,self.y_min,self.z_max),(self.x_min,self.y_max,self.z_max),(self.x_min,self.y_max,self.z_min)])
        circuits.append([(self.x_max,self.y_min,self.z_min),(self.x_max,self.y_min,self.z_max),(self.x_max,self.y_max,self.z_max),(self.x_max,self.y_max,self.z_min)])
        return [x[1] for x in sorted(Polyhedron.circuit_intersect(edge, circuits))]
        
    def delete(self, poly):
        face_lookup = dict()
        for face_index, face in enumerate(poly.faces):
            print(face)
            circuit = Polyhedron.circuit_cut(poly.circuits(face_index))
            print('circuit', circuit)
            new_edges = set()
            box = []
            if circuit[0][0] == circuit[1][0] and circuit[0][0] == circuit[2][0] and circuit[0][0] >= self.x_min and circuit[0][0] <= self.x_max:
                box = [(circuit[0][0], self.y_min, self.z_min),(circuit[0][0], self.y_max, self.z_min),(circuit[0][0], self.y_max, self.z_max),(circuit[0][0], self.y_min, self.z_max)]
            elif circuit[0][1] == circuit[1][1] and circuit[0][1] == circuit[2][1] and circuit[0][1] >= self.y_min and circuit[0][1] <= self.y_max:
                box = [(self.x_min, circuit[0][1], self.z_min),(self.x_max, circuit[0][1], self.z_min),(self.x_max, circuit[0][1], self.z_max),(self.x_min, circuit[0][1], self.z_max)]
            elif circuit[0][2] == circuit[1][2] and circuit[0][2] == circuit[2][2] and circuit[0][2] >= self.z_min and circuit[0][2] <= self.z_max:
                box = [(self.x_min, self.y_min, circuit[0][2]),(self.x_max, self.y_min, circuit[0][2]),(self.x_max, self.y_max, circuit[0][2]),(self.x_min, self.y_max, circuit[0][2])]
            intersections = []
            for i,x in enumerate(circuit):
                p1 = circuit[i-1]
                p2 = x
                last_intersections = []
                for j, box_point in enumerate(box):
                    if not Box.colinear(set([p1,p2,box[j-1],box_point])):
                        for intersection in Box.intersect_segments(frozenset([p1,p2]),frozenset([box[j-1],box_point])):
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
            print('WHO',new_edges)
            if len(intersections):
                print(box)
                for intersection in intersections:
                    for j in range(len(box)):
                        if Box.point_on_segment(frozenset([box[j-1],box[j]]), intersection[0]) and box[j-1] != intersection[0] and box[j] != intersection[0]:
                            box.insert(j,intersection[0])
                            break
                print('intersections', intersections)
                for i, intersection in enumerate(intersections):
                    index = box.index(intersection[0])
                    index += 1
                    index %= len(box)
                    path1 = [intersection[0]]
                    while True:
                        path1.append(box[index])
                        if (box[index],True) in intersections or (box[index],False) in intersections:
                            break
                        index += 1
                        index %= len(box)
                    index = box.index(intersection[0])
                    index -= 1
                    index %= len(box)
                    path2 = [intersection[0]]
                    while True:
                        path2.append(box[index])
                        if (box[index],True) in intersections or (box[index],False) in intersections:
                            break
                        index -= 1
                        index %= len(box)
                    print("path",path1)
                    print("path",path2)
                    #if True:
                    if not intersection[1]:
                        for k,z in enumerate(path1[:-1]):
                            point = tuple((z[j]+path1[k+1][j])/2 for j in range(3))
                            if not any(Polyhedron.inside_triangle(y,point) for y in Polyhedron.triangulate(circuit)):
                                break
                        else:
                            for j,x in enumerate(path1[:-1]):
                                new_edges.add(frozenset([x,path1[j+1]]))
                        for k,z in enumerate(path2[:-1]):
                            point = tuple((z[j]+path2[k+1][j])/2 for j in range(3))
                            if not any(Polyhedron.inside_triangle(y,point) for y in Polyhedron.triangulate(circuit)):
                                break
                        else:
                            for j,x in enumerate(path2[:-1]):
                                new_edges.add(frozenset([x,path2[j+1]]))
            else:
                print('box',box)
                if len(box) and any([Polyhedron.inside_triangle(x,box[0]) for x in Polyhedron.triangulate(circuit)]):
                    for j, box_point in enumerate(box):
                        new_edges.add(frozenset([box[j-1],box_point]))
                for i,x in enumerate(circuit):
                    if circuit[i-1] in poly.verts and x in poly.verts and frozenset([poly.verts.index(circuit[i-1]),poly.verts.index(x)]) in poly.edges:
                        new_edges.add(frozenset([circuit[i-1],x]))
                print('BOO')
            while True:
                combine = False
                for edge1 in new_edges:
                    for edge2 in new_edges:
                        if len(edge1&edge2) == 1 and Box.colinear(edge1|edge2):
                            new_edges.add(edge1^edge2)
                            new_edges.remove(edge1)
                            new_edges.remove(edge2)
                            combine = True
                            break
                    if combine:
                        break
                else:
                    break

            face_lookup[face] = new_edges
            print('face lookup', face, new_edges)
        faces = []
        negative = dict()
        for face in face_lookup:
            if len(face_lookup[face]):
                print('face lookup', face_lookup[face])
                print(set(frozenset(poly.verts[x] for x in poly.edges[y]) for y in face))
                negative[face] = set(frozenset(poly.verts[x] for x in poly.edges[y]) for y in face)
                new_edge = set()
                while True:
                    for edge1 in face_lookup[face]:
                        for edge2 in negative[face]:
                            if edge2 not in new_edge:
                                for point in edge1:
                                    if not Box.point_on_segment(edge2,point):
                                        break
                                else:
                                    negative[face].remove(edge2)
                                    p1, p2 = edge2
                                    edge3 = [p1, sorted((distance(p1,x),x) for x in edge1)[0][1]]
                                    if distance(edge3[0],edge3[1]) > 0:
                                        negative[face].add(frozenset(edge3))
                                        new_edge.add(frozenset(edge3))
                                    edge3 = [p2, sorted((distance(p2,x),x) for x in edge1)[0][1]]
                                    if distance(edge3[0],edge3[1]) > 0:
                                        negative[face].add(frozenset(edge3))
                                        new_edge.add(frozenset(edge3))
                                    break
                                for point in edge2:
                                    if not Box.point_on_segment(edge1,point):
                                        break
                                else:
                                    negative[face].remove(edge2)
                                    p1, p2 = edge1
                                    edge3 = [p1, sorted((distance(p1,x),x) for x in edge2)[0][1]]
                                    if distance(edge3[0],edge3[1]) > 0:
                                        negative[face].add(frozenset(edge3))
                                        new_edge.add(frozenset(edge3))
                                    edge3 = [p2, sorted((distance(p2,x),x) for x in edge2)[0][1]]
                                    if distance(edge3[0],edge3[1]) > 0:
                                        negative[face].add(frozenset(edge3))
                                        new_edge.add(frozenset(edge3))
                                    break
                        else:
                            negative[face].add(edge1)
                            new_edge.add(edge1)
                    else:
                        break
                print('negative', negative[face])
                temp = face_lookup[face]
                print('positive', temp)
                new_edge = set()
                for edge1 in negative[face]:
                    print('edge1',edge1)
                    for face2 in face_lookup:
                        double_break = False
                        for edge2 in face_lookup[face2]:
                            print('edge2',edge2)
                            for point in edge1:
                                if not Box.point_on_segment(edge2,point):
                                    break
                            else:
                                print(edge1,edge2, face2, face)
                                if face2 != face:
                                    temp.add(edge1)
                                    new_edge.add(edge1)
                                    double_break = True
                                    break
                                else:
                                    if edge2 not in new_edge:
                                        temp.remove(edge2)
                                        p1, p2 = edge2
                                        edge3 = [p1, sorted((distance(p1,x),x) for x in edge1)[0][1]]
                                        if distance(edge3[0],edge3[1]) > 0:
                                            temp.add(frozenset(edge3))
                                        edge3 = [p2, sorted((distance(p2,x),x) for x in edge1)[0][1]]
                                        if distance(edge3[0],edge3[1]) > 0:
                                            temp.add(frozenset(edge3))
                                        double_break = True
                                        break
                        if double_break:
                            break
                while True:
                    combine = False
                    for edge1 in temp:
                        for edge2 in temp:
                            if len(edge1&edge2) == 1 and Box.colinear(edge1|edge2):
                                combine = True
                                break
                        if combine:
                            break
                    if combine:
                        temp.add(edge1^edge2)
                        temp.remove(edge1)
                        temp.remove(edge2)
                    else:
                        break
                print(temp)
                faces.append(temp)
            else:
                faces.append(set(frozenset(poly.verts[x] for x in poly.edges[y]) for y in face))
        verts = list(set(point for face in faces for edge in face for point in edge))
        edges = list(set(edge for face in faces for edge in face))
        print('BLAH',edges)
        faces = [frozenset(edges.index(edge) for edge in face) for face in faces]
        edges = [frozenset(verts.index(point) for point in edge) for edge in edges]
        new_poly = Polyhedron()
        new_poly.verts = verts
        new_poly.edges = edges
        new_poly.faces = faces
        for face_index, face in enumerate(new_poly.faces):
            print('face', [frozenset(verts[y] for y in edges[x]) for x in face])
            new_poly.circuits(face_index)
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
                    edge3 = [p1, sorted([(distance(p1,x),x) for x in edge1])[0][1]]
                    print('edge3', edge3)
                    if distance(edge3[0],edge3[1]) > 0:
                        new_new_edges.add(frozenset(edge3))
                        for face in faces:
                            if edge1 in face or edge2 in face:
                                face.add(frozenset(edge3))
                    edge3 = [p2, sorted([(distance(p2,x),x) for x in edge1])[0][1]]
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
                    edge3 = [p1, sorted([(distance(p1,x),x) for x in edge2])[0][1]]
                    print('edge3', edge3)
                    if distance(edge3[0],edge3[1]) > 0:
                        new_new_edges.add(frozenset(edge3))
                        for face in faces:
                            if edge1 in face or edge2 in face:
                                face.add(frozenset(edge3))
                    edge3 = [p2, sorted([(distance(p2,x),x) for x in edge2])[0][1]]
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
                    if len(edge1&edge2)==1 and Box.colinear(edge1|edge2):
                        print(edge1, edge2)
                        combine = True
                        break
                if combine:
                    edges.remove(edge2)
                    for face in faces:
                        if edge2 in face:
                            face.remove(edge2)
                            face.add(edge1^edge2)
                    new_edges.remove(edge1)
                    print("combine",edge1^edge2)
                    new_edges.add(edge1^edge2)
                    break
            else:
                break
        for edge1 in new_edges:
            for face in faces:
                if Box.coplanar(set(point for edge in face for point in edge)|set(edge1)):
                    face.add(edge1)
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
                new_faces[new_face_bool.index(True)].update(face)
                del_faces.add(face)

        print('new_edges', new_edges)
        for edge in new_edges|edges:
            edge = list(edge)
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

        '''
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
        '''
        for face in faces:
            print('HEY', face)
        '''
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
        '''
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
            print(face_index, face)
            print(new_poly.circuits(face_index))
        return new_poly
            
        

        
class Polyhedron:
    def __init__(self):
        self.verts = []
        self.edges = []
        self.faces = []
    def circuits(self, face_index):
        face = self.faces[face_index]
        print([self.edges[x] for x in face])
        if not len(face):
            return []
        edge_lookup = dict()
        for edge_index in face:
            edge = self.edges[edge_index]
            for index in self.edges[edge_index]:
                if index not in edge_lookup:
                    edge_lookup[index] = set()
                edge_lookup[index].add(edge)
        circuits = set()
        edges = {self.edges[x] for x in self.faces[face_index]}
        #print(self.faces[face_index])
        while len(edges):
            start = list(edges)[0]
            edges.remove(start)
            previous = start
            point = list(start)[0]
            circuit = []
            circuit.append(self.verts[point])
            #print([[self.verts[index] for index in self.edges[edge_index]] for edge_index in face])
            print(edge_lookup[point], set([start]))
            current = list(edge_lookup[point] - set([start]))[0]
            #print([[self.verts[index] for index in edge] for edge in self.edges])
            while current != start:
                edges.remove(current)
                point = list(current - previous)[0]
                circuit.append(self.verts[point])
                previous = current
                #print(point, [[self.verts[index] for index in edge] for edge in edge_lookup[point]])
                #print('circuit', circuit)
                for p in edge_lookup:
                    if current in edge_lookup[p]:
                        edge_lookup[p].remove(current)
                current = list(edge_lookup[point])[0]
            circuits.add(tuple(circuit))
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
        print(circuit, is_convex)
        for i in range(len(circuit)):
            if is_convex[i]:
                for j,y in enumerate(circuit):
                    if j != i and j != (i+1)%len(circuit) and j != (i-1)%len(circuit):
                        if Polyhedron.inside_triangle([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]],y):
                            print([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]],y)
                            break
                else:
                    return tuple([circuit[i-1],circuit[i],circuit[(i+1)%len(circuit)]]), tuple([x for k,x in enumerate(circuit) if k != i])
    def triangulate(circuit):
        print(circuit)
        output = []
        remainder = circuit
        while len(remainder) > 3:
            ear, remainder = Polyhedron.clip_ear(remainder)
            output.append(ear)
        output.append(remainder)
        return output
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
                        alpha = 0
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
    def draw(self, pygame, screen):
        screen_width, screen_height = screen.get_size()
        for edge in self.block.edges:
            p1, p2 = tuple(self.block.verts[index] for index in edge)
            #print(p1,p2)
            p1, p2 = camera.project(p1), camera.project(p2)
            p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
            p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
            pygame.draw.line(screen, "white", p1, p2)
        print('faces', self.poly.faces)
        for face_index,face in enumerate(self.poly.faces):
            points = Polyhedron.circuit_cut(self.poly.circuits(face_index))
            centroid = tuple(sum(points[i][j] for i,x in enumerate(points))/len(points) for j in range(3))
            distance = sum((centroid[i]-camera.focal[i])**2 for i in range(3))
            points = [camera.project(x) for x in points]
            points = [(x[0]*1+screen_width/2,x[1]*-1+screen_height/2) for x in points]
            color = "white"
            pygame.draw.polygon(screen, color, points)
        for edge in self.poly.edges:
            p1, p2 = tuple(self.poly.verts[index] for index in edge)
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
    def subdivide(self, mult):
        block = Block(self.size[0], self.size[1], self.size[2], self.unit*mult)
        block.select = self.select
        block.poly = self.poly
        return block

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
    delete = False
    dir_mult1 = 1
    dir_mult2 = 1
    dir_mult3 = 1
    z_forward = True
    block = Block(3,3,3)
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
                if event.key == pygame.K_2:
                    block = block.subdivide(0.5)
                if event.key == pygame.K_3:
                    block = block.subdivide(1/3)
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
                            #while not block.select_by_void():
                            #    block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                        elif block.select[3] == 1:
                            block.select[1] += direction*block.unit
                            #while not block.select_by_void():
                            #    block.select[1] += direction*block.unit
                        elif block.select[3] == 2:
                            block.select[1] += direction*block.unit
                            #while not block.select_by_void():
                            #    block.select[1] += direction*block.unit
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
                            #while not block.select_by_void():
                            #    block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(2*int(z_forward)-1)*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                        elif block.select[3] == 1:
                            block.select[0] += direction*dir_mult1*block.unit
                            #while not block.select_by_void():
                            #    block.select[0] += direction*dir_mult1*block.unit
                        elif block.select[3] == 2:
                            block.select[2] -= direction*dir_mult2*block.unit
                            #while not block.select_by_void():
                            #    block.select[2] -= direction*dir_mult2*block.unit
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
                    box = Box((i,i+size[0]),(j,j+size[1]),(k,k+size[2]))
                    print('box',(i,i+size[0]),(j,j+size[1]),(k,k+size[2]))
                    block.poly, anti_poly = box.delete(block.poly)
                    print("anti_poly")
                    print(anti_poly.verts)
                    if not delete:
                        block.poly = box.add(block.poly)
                    block.select = [block.select[0]+block.select_size[0]-block.unit,block.select[1]+block.select_size[1]-block.unit,block.select[2]+block.select_size[2]-block.unit,block.select[3]]
                    block.select_size = [block.unit,block.unit,block.unit]
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
