from utility import *
from camera import Camera
from polyhedron import Polyhedron, get_cube

poly1 = get_cube(factors=2, angles=(pi/4,pi/4,pi/4))
poly2 = get_cube((1,1,1), factors=2)

edges = set()
with open("intersects.txt") as f:
    for line in f:
        if not(len(line.strip())):
            continue
        p1_text, p2_text = line.strip().split(" ")
        p1 = tuple(float(x) for x in p1_text.strip('[').strip(']').split(','))
        p2 = tuple(float(x) for x in p2_text.strip('[').strip(']').split(','))
        edges.add(frozenset([p1,p2]))

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
    pygame.display.flip()
