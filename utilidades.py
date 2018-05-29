from math import acos, degrees
from functools import lru_cache

def segment_intersection_test(a,b):
	"""Devuelve True si el segmento a se corta con el segmento b"""
	A,B,C,D = a[0],a[1],b[0],b[1]
	if (in_segment(A,b) or in_segment(B,b) or in_segment(C,a) or in_segment(D,a)):
		return True
	return (sarea(A,B,C)*sarea(A,B,D))<0 and (sarea(C,D,A)*sarea(C,D,B)<0)

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

def in_segment(P,s):
	"""Devuelve True si el punto P esta en el segmento s"""
	A,B = s[0],s[1]
	if sarea(P,A,B) != 0: return False
	if A[0] == B[0]: return min(A[1],B[1]) <= P[1] <= max(A[1],B[1])
	return min(A[0],B[0]) <= P[0] <= max(A[0],B[0])

def pointInPolygon(Q,p):
	"""Devuelve True si el punto Q esta dentro del poligono p"""
	Q = list(Q)
	count = 0
	disparo = [Q,[max(p)[0]+1,Q[1]]]
	n = len(p)
	for i in range(len(p)):
		if in_segment(Q,[p[i],p[(i+1)%n]]):
			return True
		if segment_intersection_test(disparo,[p[i],p[(i+1)%n]]):
			if in_segment(p[i],disparo):
				if sarea(Q,p[i],p[i-1])*sarea(Q,p[i],p[(i+1)%n]) < 0:
					count += 1
			else:
				count += 1
	return count%2==1

def get_angles(a,b,c):
    get_angle = lambda a, b, c: degrees(acos((b*b+c*c-a*a)/(float(2*b*c))))
    if a + b <= c or b + c <= a or c + a <= b:
        return [0, 0, 0]
    return [get_angle(a, b, c), get_angle(b, c, a), get_angle(c, a, b)]