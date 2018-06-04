from math import acos, degrees
from functools import lru_cache
import matplotlib.path as mpltPath

def sarea(A,B,C):
    """Calcula el area signada del triangulo formado por los puntos ordenados A,B,C"""
    A = tuple(A)
    B = tuple(B)
    C = tuple(C)
    return aux_sarea(A,B,C)

@lru_cache(maxsize=100)
def aux_sarea(A,B,C):
    return ((B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1]))/2


def segment_crossing(a,b):
    A,B,C,D = a[0],a[1],b[0],b[1]
    return (sarea(A,B,C)*sarea(A,B,D))<0 and (sarea(C,D,A)*sarea(C,D,B)<0)

def pointInPolygon(Q,p):
    """Devuelve True si el punto Q esta dentro del poligono p"""
    path = mpltPath.Path(p)
    inside = path.contains_points([Q])
    return inside

def get_angles(a,b,c):
    get_angle = lambda a, b, c: degrees(acos((b*b+c*c-a*a)/(float(2*b*c))))
    if a + b <= c or b + c <= a or c + a <= b:
        return [0, 0, 0]
    return [get_angle(a, b, c), get_angle(b, c, a), get_angle(c, a, b)]