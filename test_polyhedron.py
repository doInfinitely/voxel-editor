import random
from polyhedron import Polyhedron
from utility import *
from camera import Camera
def get_possible_voxels(voxel_poly):
    displacements = []
    for i in range(3):
        for x in [-1,1]:
            d = [0,0,0]
            d[i] = x
            displacements.append(tuple(d))
    occupied = set()
    for i,x in enumerate(voxel_poly):
        for j,y in enumerate(x):
            for k,z in enumerate(y):
                if z == 1:
                    occupied.add((i,j,k))
    possible = {(x[0]+d[0],x[1]+d[1],x[2]+d[2]) for x in occupied for d in displacements}-occupied
    return {x for x in possible if x[0] >= 0 and x[0] < len(voxel_poly) and x[1] >= 0 and x[1] < len(voxel_poly[0]) and x[2] >= 0 and x[2] < len(voxel_poly[0][0])}
'''
def voxel_to_polyhedron(voxel_poly):
    shape = (len(voxel_poly), len(voxel_poly[0]), len(voxel_poly[0][0]))
    displacements = []
    for x in [0,1]:
        for y in [0,1]:
            for z in [0,1]:
                displacements.append((x,y,z))
    points = [[[0 for k in range(shape[2]+1)] for j in range(shape[1]+1)] for i in range(shape[0]+1)]
    for i,x in enumerate(voxel_poly):
        for j,y in enumerate(x):
            for k,z in enumerate(y):
                if z == 1:
                    for d in displacements:
                        coord = (i+d[0],j+d[1],k+d[2])
                        points[coord[0]][coord[1]][coord[2]] += 1
    verts = []
    for i,x in enumerate(points):
        for j,y in enumerate(x):
            for k,z in enumerate(y):
                if z % 2 == 1:
                    verts.append((i,j,k))
    edges = []
    for index1, vert1 in enumerate(verts):
        for index2, vert2 in enumerate(verts):
            if index1 < index2 and sum(int(v1==v2) for v1,v2 in zip(vert1,vert2)) == 2:
                print(vert1,vert2)
                index = [v1==v2 for v1,v2 in zip(vert1,vert2)].index(False)
                for w in range(min(vert1[index],vert2[index])+1,max(vert1[index],vert2[index])):
                    point = list(vert1)
                    point[index] = w
                    print(point, points[point[0]][point[1]][point[2]])
                    if points[point[0]][point[1]][point[2]] % 4 != 2:
                        break
                else:
                    if points[vert1[0]][vert1[1]][vert1[2]] == 1 or points[vert2[0]][vert2[1]][vert2[2]] == 1:
                        edges.append(frozenset([verts.index(vert1),verts.index(vert2)]))
    chains = dict()
    for index1, edge1 in enumerate(edges):
        for index2, edge2 in enumerate(edges):
            if index1 < index2 and len(edge1&edge2) == 1 and Polyhedron.colinear({verts[index] for index in edge1}|{verts[index] for index in edge2}):
                if edge1 not in chains:
                    chains[edge1] = []
                chains[edge1].append(edge2)
                if edge2 not in chains:
                    chains[edge2] = []
                chains[edge2].append(edge1)
    components = set()
    for key in chains:
        q = [key]
        component = set()
        while len(q):
            item = q.pop()
            if item not in component:
                component.add(item)
                q.extend(chains[item])
        components.add(frozenset(component))
    for component in components:
        component = sorted(component, key=lambda x: min(frozenset(verts[index] for index in x)))
        start = min(frozenset(verts[index] for index in component[0]))
        if points[start[0]][start[1]][start[2]] != 1:
            component = list(reversed(component))
        for i,x in enumerate(component):
            if i % 2 == 1:
                edges.remove(x)
    index1 = 0
    while index1 < len(edges):
        index2 = 0
        while index2 < len(edges):
            intersect = Polyhedron.intersect_segments(frozenset(verts[index] for index in edges[index1]),frozenset(verts[index] for index in edges[index2]))
            if index1 < index2 and intersect is not None and intersect not in frozenset(verts[index] for index in edges[index1])|frozenset(verts[index] for index in edges[index2]):
                p1, p2 = frozenset(verts[index] for index in edges[index1]) 
                p3, p4 = frozenset(verts[index] for index in edges[index2])
                verts.append(intersect)
                edges[index1] = frozenset([verts.index(p1), len(verts)-1])
                edges[index2] = frozenset([verts.index(p3), len(verts)-1])
                edges.append(frozenset([verts.index(p2), len(verts)-1]))
                edges.append(frozenset([verts.index(p4), len(verts)-1]))
            index2 += 1
        index1 += 1
    print(edges)
    coplanar_edges = [{edge1, edge2} for edge1 in edges for edge2 in edges if len(edge1&edge2) == 1]
    updated = True
    while updated:
        updated = False
        for index1, edge in enumerate(edges):
            for index2, x in enumerate(coplanar_edges):
                if edge not in x and Polyhedron.coplanar({verts[index] for index in edge}|{verts[index] for y in x for index in y}):
                    coplanar_edges[index2].add(edge)
                    updated = True
    coplanar_edges = {frozenset(x) for x in coplanar_edges}
    face_circuits = []
    for x in coplanar_edges:
        print({frozenset(verts[index] for index in edge) for edge in x})
        circuits = Polyhedron.circuit_helper({frozenset(verts[index] for index in edge) for edge in x})
        print(circuits)
        circuit_sets = [{circuit} for circuit in circuits]
        updated = True
        while updated:
            updated = False
            for index1, circuits1 in enumerate(circuit_sets):
                for index2, circuits2 in enumerate(circuit_sets):
                    if index1 < index2 and Polyhedron.find_exterior_circuit(circuits1|circuits2) is not None:
                        circuit_sets[index1] |= circuits2
                        del circuit_sets[index2]
                        updated = True
                        break
                if updated:
                    break
        face_circuits.extend(circuit_sets)
    faces = []
    for circuits in face_circuits:
        face = set()
        for circuit in circuits:
            for i, p2 in enumerate(circuit):
                p1 = circuit[i-1]
                face.add(edges.index(frozenset([verts.index(p1),verts.index(p2)])))
        faces.append(frozenset(face))
    verts = [tuple(float(x) for x in vert) for vert in verts]
    poly = Polyhedron()
    poly.verts = verts
    poly.edges = edges
    poly.faces = faces
    return poly
'''
def voxel_to_polyhedron(voxel_poly):
    shape = (len(voxel_poly), len(voxel_poly[0]), len(voxel_poly[0][0]))
    faces = []
    for i1 in range(3):
        i2 = (i1 + 1)%3
        i3 = (i1 + 2)%3
        indices = [0,0,0]
        while indices[i1] < shape[i1]+1:
            components = dict()
            indices[i2] = 0
            while indices[i2] < shape[i2]:
                indices[i3] = 0
                while indices[i3] < shape[i3]:
                    indices_temp = list(indices)
                    indices_temp[i1] -= 1
                    if (indices[i1] == shape[i1] and voxel_poly[indices_temp[0]][indices_temp[1]][indices_temp[2]] == 1) or (indices_temp[i1] == -1 and voxel_poly[indices[0]][indices[1]][indices[2]] == 1) or (indices[i1] > 0 and indices[i1] < shape[i1] and voxel_poly[indices[0]][indices[1]][indices[2]]+voxel_poly[indices_temp[0]][indices_temp[1]][indices_temp[2]] == 1):
                        square = (indices[i2],indices[i3])
                        components[square] = frozenset([square])
                    indices[i3] += 1
                indices[i2] += 1
            for key1 in components:
                displacements = [(0,-1),(0,1),(-1,0),(1,0)]
                for d in displacements:
                    neighbor = (key1[0]+d[0], key1[1]+d[1])
                    if neighbor in components:
                        components[neighbor] |= components[key1]
                        for key2 in components[neighbor]:
                            components[key2] = components[neighbor]
            for component in set(components.values()):
                segments = set()
                for square in component:
                    displacements = [(0,0,1,0),(0,0,0,1),(1,1,-1,0),(1,1,0,-1)]
                    for d in displacements:
                        segment = [[0,0,0],[0,0,0]]
                        segment[0][i1] = indices[i1]
                        segment[0][i2] = square[0]+d[0]
                        segment[0][i3] = square[1]+d[1]
                        segment[1][i1] = segment[0][i1]
                        segment[1][i2] = segment[0][i2]+d[2]
                        segment[1][i3] = segment[0][i3]+d[3]
                        segment = frozenset(tuple(x) for x in segment)
                        if segment in segments:
                            segments.remove(segment)
                        else:
                            segments.add(segment)
                segments = list(segments)
                index1 = 0
                while index1 < len(segments):
                    index2 = index1+1
                    did_break = False
                    while index2 < len(segments):
                        if len(segments[index1]&segments[index2]) == 1 and Polyhedron.colinear(segments[index1]|segments[index2]):
                            segments.append(segments[index1]^segments[index2])
                            del segments[index2]
                            del segments[index1]
                            did_break = True
                            break
                        index2 += 1
                    if not did_break:
                        index1 += 1
                if len(segments):
                    faces.append(frozenset(segments))
            indices[i1] += 1
    verts = list(set(point for face in faces for edge in face for point in edge))
    edges = [frozenset(verts.index(point) for point in edge) for face in faces for edge in face]
    faces = [frozenset(edges.index(frozenset(verts.index(point) for point in edge)) for edge in face) for face in faces]
    verts = [tuple(float(x) for x in vert) for vert in verts]
    poly = Polyhedron()
    poly.verts = verts
    poly.edges = edges
    poly.faces = faces
    return poly

def assert_polys_equal(poly1, poly2):
    print(poly1.verts, poly2.verts)
    print(poly1.edges, poly2.edges)
    print(poly1.faces, poly2.faces)
    faces1 = {frozenset(frozenset(tuple(round(x) for x in poly1.verts[index]) for index in poly1.edges[edge_index]) for edge_index in face) for face in poly1.faces}
    faces2 = {frozenset(frozenset(tuple(round(x) for x in poly2.verts[index]) for index in poly2.edges[edge_index]) for edge_index in face) for face in poly2.faces}
    for face in faces1:
        print(face)
    print()
    for face in faces2:
        print(face)
    print()
    for face in faces1-faces2:
        print(face)
    print()
    for face in faces2-faces1:
        print(face)
    print()
    return list(faces1-faces2), list(faces2-faces1)
    assert faces1 == faces2
def test_polyhedron(count1=1000, count2=None, seed=0, shape=(100,100,100)):
    if count2 is None:
        count2 = count1
    assert count1 <= shape[0]*shape[1]*shape[2]
    assert count2 <= shape[0]*shape[1]*shape[2]
    counts = [count1, count2]
    random.seed(seed)
    voxel_poly1 = [[[0 for k in range(shape[2])] for j in range(shape[1])] for i in range(shape[0])]
    voxel_poly2 = [[[0 for k in range(shape[2])] for j in range(shape[1])] for i in range(shape[0])]
    start = tuple(random.randrange(x) for x in shape)
    for voxel_poly in [voxel_poly1, voxel_poly2]:
        voxel_poly[start[0]][start[1]][start[2]] = 1
    voxel_polys = [voxel_poly1, voxel_poly2]
    for index in range(2):
        for _ in range(counts[index]):
            v = random.choice(list(get_possible_voxels(voxel_polys[index])))
            voxel_polys[index][v[0]][v[1]][v[2]] = 1
    voxel_intersection = [[[int(voxel_poly1[i][j][k] and voxel_poly2[i][j][k]) for k in range(shape[2])] for j in range(shape[1])] for i in range(shape[0])]
    voxel_sum = [[[int(voxel_poly1[i][j][k]+voxel_poly2[i][j][k]>0) for k in range(shape[2])] for j in range(shape[1])] for i in range(shape[0])]
    voxel_difference1 = [[[int(voxel_poly1[i][j][k]-voxel_poly2[i][j][k]>0) for k in range(shape[2])] for j in range(shape[1])] for i in range(shape[0])]
    voxel_difference2 = [[[int(voxel_poly2[i][j][k]-voxel_poly1[i][j][k]>0) for k in range(shape[2])] for j in range(shape[1])] for i in range(shape[0])]
    poly1 = voxel_to_polyhedron(voxel_poly1)
    poly2 = voxel_to_polyhedron(voxel_poly2)
    faces1 = {frozenset(frozenset(poly1.verts[index] for index in poly1.edges[edge_index]) for edge_index in face) for face in poly1.faces}
    faces2 = {frozenset(frozenset(poly2.verts[index] for index in poly2.edges[edge_index]) for edge_index in face) for face in poly2.faces}
    for face in faces1:
        print(face)
    print()
    for face in faces2:
        print(face)
    poly3 = voxel_to_polyhedron(voxel_sum) 
    #return poly1, poly2, poly3, poly2.add(poly1)
    assert_polys_equal(voxel_to_polyhedron(voxel_intersection), poly1.intersect(poly2))
    #assert_polys_equal(poly2.add(poly1), poly1.add(poly2))
    assert_polys_equal(voxel_to_polyhedron(voxel_sum), poly2.add(poly1))
    assert_polys_equal(voxel_to_polyhedron(voxel_difference1), poly1.subtract(poly2))
    assert_polys_equal(voxel_to_polyhedron(voxel_difference2), poly2.subtract(poly1))

if __name__ == "__main__":
    test_polyhedron(14,14,0,(3,3,3))
    for i in range(18, 1000):
        #test_polyhedron(seed = i)
        test_polyhedron(14,14,i,(3,3,3))
        print(f'test {i} passed')
if False:
    poly1, poly2, poly3, poly4 = test_polyhedron(14,14,0,(3,3,3))
    #for i in range(1000):
    #    test_polyhedron(seed = i)
    #    print(f'test {i} passed')
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
    print('edges',poly3.edges)
    #colors = {face:tuple(random.randrange(256) for i in range(3)) for face in poly3.faces}
    count = 0
    temp = 0
    temp1 = 0
    temp2 = 0
    faces1, faces2 = assert_polys_equal(poly3, poly4)
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
        '''
        for edge in poly3.edges:
            p1, p2 = tuple(poly3.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "white", p1, p2)
        '''
        for edge in poly4.edges:
            p1, p2 = tuple(poly4.verts[index] for index in edge)
            p1, p2 = camera.project(p1), camera.project(p2)
            if p1 and p2:
                p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                #print(p1,p2)
                pygame.draw.line(screen, "yellow", p1, p2)
        '''
        if count % 1000 == 0:
            temp += 1
            temp %= len(poly4.faces)
        for face_index, face in enumerate(poly4.faces):
            if face_index != temp:
                continue
            for edge_index in face:
                edge = poly4.edges[edge_index]
                p1, p2 = tuple(poly4.verts[index] for index in edge)
                p1, p2 = camera.project(p1), camera.project(p2)
                if p1 and p2:
                    p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                    p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                    #print(p1,p2)
                    pygame.draw.line(screen, "green", p1, p2)
        '''
        '''
        if count % 1000 == 0:
            temp1 += 1
            temp1 %= len(faces1)
        for face_index, face in enumerate(faces1):
            if face_index != temp1:
                continue
            #color = colors[face]
            for edge in face:
                p1, p2 = edge
                p1, p2 = camera.project(p1), camera.project(p2)
                if p1 and p2:
                    p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                    p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                    #print(p1,p2)
                    pygame.draw.line(screen, "green", p1, p2)
        '''
        if count % 1000 == 0:
            temp2 += 1
            temp2 %= len(faces2)
        for face_index, face in enumerate(faces2):
            if face_index != temp2:
                continue
            #color = colors[face]
            for edge in face:
                p1, p2 = edge
                p1, p2 = camera.project(p1), camera.project(p2)
                if p1 and p2:
                    p1 = (p1[0]*1+screen_width/2,p1[1]*-1+screen_height/2)
                    p2 = (p2[0]*1+screen_width/2,p2[1]*-1+screen_height/2)
                    #print(p1,p2)
                    pygame.draw.line(screen, "yellow", p1, p2)
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
