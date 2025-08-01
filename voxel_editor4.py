from utility import *
from polyhedron import Polyhedron, get_cube
from camera import Camera

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
        print("verts",self.poly.verts)
        for face_index,face in enumerate(self.poly.faces):
            print("face",{frozenset(self.poly.verts[index] for index in self.poly.edges[edge_index]) for edge_index in face})
            points = Polyhedron.circuit_cut(self.poly.circuits(face_index))
            print("face",points)
            print('face', face)
            print()
            if points is None:
                continue
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
            #print('points',points)
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
        min_select = tuple(min(self.select[i],self.select[i]+self.select_size[i]) for i in range(3))
        max_select = tuple(max(self.select[i],self.select[i]+self.select_size[i]) for i in range(3))
        meter = 1
        for x in self.meters:
            meter *= x
        for d1 in range(3):
            for mult in [1,-1]:
                d2, d3 = [i for i in range(3) if i != d1]
                points = [[[0 for i in range(3)] for j in range(abs(int(self.select_size[d3]/meter)))] for k in range(abs(int(self.select_size[d2]/meter)))]
                if mult == 1:
                    select1 = max_select
                    select2 = min_select
                else:
                    select1 = min_select
                    select2 = max_select
                for i in range(len(points)):
                    for j in range(len(points[i])):
                        points[i][j][d2] = select2[d2] + mult*meter*(i+.5)
                        points[i][j][d3] = select2[d3] + mult*meter*(j+.5)
                        points[i][j][d1] = mult*float('inf')
                for face_index,face in enumerate(self.poly.faces):
                    circuit = Polyhedron.circuit_cut(Polyhedron.make_clockwise(self.poly.circuits(face_index)))
                    if circuit is None:
                        continue
                    if circuit[0][d1] == circuit[1][d1] and circuit[1][d1] == circuit[2][d1]:
                        for point_row in points:
                            for point in point_row:
                                projection = list(point)
                                projection[d1] = circuit[0][d1]
                                projection = tuple(projection)
                                if mult*circuit[0][d1] < mult*point[d1] and mult*circuit[0][d1] >= mult*select1[d1] and any(Polyhedron.inside_triangle(x,projection) for x in Polyhedron.triangulate(circuit)):
                                    #print(Polyhedron.triangulate(circuit))
                                    point[d1] = circuit[0][d1]
                #print(d1, points)
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
                                    point[d2] = select2[d2] + mult*i*meter
                                    point[d3] = select2[d3] + mult*j*meter
                                    point = tuple(point)
                                    point_set.add(point)
                                if (i-1,j) not in component or (i,j+1) not in component or (i-1,j+1) not in component:
                                    point = list(points[i][j])
                                    point[d2] = select2[d2] + mult*i*meter
                                    point[d3] = select2[d3] + mult*(j+1)*meter
                                    point = tuple(point)
                                    point_set.add(point)
                                if (i+1,j) not in component or (i,j-1) not in component or (i+1,j-1) not in component:
                                    point = list(points[i][j])
                                    point[d2] = select2[d2] + mult*(i+1)*meter
                                    point[d3] = select2[d3] + mult*j*meter
                                    point = tuple(point)
                                    point_set.add(point)
                                if (i+1,j) not in component or (i,j+1) not in component or (i+1,j+1) not in component:
                                    point = list(points[i][j])
                                    point[d2] = select2[d2] + mult*(i+1)*meter
                                    point[d3] = select2[d3] + mult*(j+1)*meter
                                    point = tuple(point)
                                    point_set.add(point)
                    point = tuple(float('inf') for i in range(3))
                    for x in point_set:
                        if (x[d2],x[d3]) < (point[d2],point[d3]):
                            point = x
                    path = []
                    rotated_point_mapping = {x:x for x in point_set}
                    reverse_point_mapping = {rotated_point_mapping[key]:key for key in rotated_point_mapping}
                    angle = [0,0,0]
                    while not len(path) or point != path[0]:
                        path.append(point)
                        keys = [key for key in rotated_point_mapping]
                        if len(path) > 1:
                            angle[d1] += math.atan2(reverse_point_mapping[path[-1]][d3],reverse_point_mapping[path[-1]][d2])
                            rotated_point_mapping = {x:rotated_point_mapping[keys[i]] for i,x in enumerate(rotate([tuple(y[j]-reverse_point_mapping[path[-1]][j] for j in range(3)) for y in keys], angle))}
                        else:
                            rotated_point_mapping = {tuple(x[j]-reverse_point_mapping[path[-1]][j] for j in range(3)):rotated_point_mapping[x] for x in keys}
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
                        if Polyhedron.colinear({path[i-1],path[i],path[(i+1)%len(path)]}):
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
            if not self.poly.is_inside(point):
                return True
        for vert in self.poly.verts:
            if self.select[0] <= vert[0] and vert[0] <= self.select[0]+self.unit and self.select[1] <= vert[1] and vert[1] <= self.select[1]+self.unit and self.select[2] <= vert[2] and vert[2] <= self.select[2]+self.unit:
                return True
        edges = {frozenset([x,y]) for x in points for y in points if sum(z!=y[k] for k,z in enumerate(x))==1}
        for edge1 in edges:
            for edge2 in [frozenset(self.poly.verts[index] for index in edge) for edge in self.poly.edges]:
                if all(Polyhedron.point_on_segment(edge1, point) for point in edge2) or all(Polyhedron.point_on_segment(edge2, point) for point in edge1):
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
                if Polyhedron.coplanar(set(face1)|{self.poly.verts[index] for edge_index in face2 for index in self.poly.edges[edge_index]}) and Polyhedron.circuit_overlap((face1,),self.poly.circuits(face_index)):
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
                            #while not block.select_by_void():
                            #   block.select[index] += direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
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
                            if block.select_size[index] == 0:
                                block.select_size[index] = direction*dir_mult1*dir_mult3*block.unit*(1-2*int(z_forward and dir_mult1 != dir_mult2))
                                if block.select_size[index] == block.unit:
                                    block.select[index] -= block.unit
                                    block.select_size[index] += block.unit
                                else:
                                    block.select[index] += block.unit
                                    block.select_size[index] -= block.unit
                        elif block.select[3] == 1 or block.select[3] == 2:
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
                                if block.select_size[index] == block.unit:
                                    block.select[index] -= block.unit
                                    block.select_size[index] += block.unit
                                else:
                                    block.select[index] += block.unit
                                    block.select_size[index] -= block.unit
                        elif block.select[3] == 1:
                            block.select_size[0] += direction*dir_mult1*block.unit
                            if block.select_size[0] == 0:
                                block.select_size[0] = direction*dir_mult1*block.unit
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
                    size = (abs(block.select_size[0]),abs(block.select_size[1]),abs(block.select_size[2]))
                    box = get_cube((i+size[0]/2,j+size[1]/2,k+size[2]/2), size)
                    #print('box',(i,i+size[0]),(j,j+size[1]),(k,k+size[2])) 
                    if delete:
                        block.poly = block.poly.subtract(box)
                    else:
                        block.poly = block.poly.add(box)
                    block.poly.round_verts(lambda x: round_point_meter(x, meters))
                    block.flip()
                    #print(meters)
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
