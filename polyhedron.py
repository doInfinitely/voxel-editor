import utility
import copy

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
    
    def clockwise_angles(planar):
        output = []
        #temp = []
        for i in range(len(planar)):
            planar = [tuple(y[j]-planar[i][j] for j in range(3)) for y in planar]
            angle = (0,0,-math.atan2(planar[i][1]-planar[i-1][1],planar[i][0]-planar[i-1][0]))
            planar = rotate(planar, angle)
            theta = math.atan2(planar[(i+1)%len(planar)][1],planar[(i+1)%len(planar)][0])
            output.append(theta < 0)
            print(planar, theta/math.pi*180)
            #temp.append(math.acos(dot(planar[i-1],planar[(i+1)%len(planar)])/distance((0,0,0), planar[i-1])/distance((0,0,0),planar[(i+1)%len(planar)]))/math.pi*180)
        #print(temp)
        return output
    def make_planar(circuit):
        p1 = (distance(circuit[0],circuit[1]),0,0)
        angle = [0,0, math.asin((circuit[1][1]-circuit[0][1])/p1[0])]
        p1 = rotate([p1], angle)[0]
        #print(p1,tuple(circuit[1][i]-circuit[0][i] for i in range(3)))
        angle[1] = math.asin((circuit[1][2]-circuit[0][2])/p1[0])
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
        print(circuit, is_convex)
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
                        print(triangulation)
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
                solutions = solve(eqs, dict=True)
                if len(solutions):
                    alpha, beta, gamma = solutions[0][alpha], solutions[0][beta], solutions[0][gamma]
                #(beta*(alpha*triangle[0][i]+(1-alpha)*triangle[1][i])+(1-beta)*triangle[2][i]-position[i]+direction[i]*gamma)**2
                    if alpha >= 0 and alpha <= 1 and beta >= 0 and beta <= 1 and gamma > 0:
                        p = tuple(alpha*triangle[0][i]+(1-alpha)*triangle[1][i] for i in range(3))
                        p = tuple(beta*x+(1-beta)*triangle[2][i] for i,x in enumerate(p))
                        if Polyhedron.inside_triangle(triangle,p):
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
    def union(self, other):
        poly1 = self
        poly
