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
