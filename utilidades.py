import numpy as np
from functools import cmp_to_key
from math import acos, degrees
import matplotlib.pyplot as plt
from matplotlib import collections  as mc

def dist2(A,B):
	"""Calcula el cuadrado de la distancia entre los puntos A y B"""
	return (A[0]-B[0])**2+(A[1]-B[1])**2

def orientation(A,B,C):
	"""Calcula la orientacion del triangulo formado por los puntos A,B,C"""
	return np.sign(((B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1])))

def angular_sort(p,c):
	"""Devuelve la lista p ordenada angularmente con respecto al punto c"""
	def compare(p,q):
		if p==q:
			return int(0)
		if p[1]>c[1]>=q[1]:
			return int(1)
		if q[1]>c[1]>=p[1]:
			return int(-1)
		orient = orientation(c,p,q)
		if orient==1:
			return int(1)
		if orient==-1:
			return int(-1)
		if dist2(c,p)<dist2(c,q):
			return int(1)
		return int(-1)
	q = p[:]
	q = sorted(q,key=cmp_to_key(compare))
	q.reverse()
	return q

def segment_intersection_test(a,b):
	"""Devuelve True si el segmento a se corta con el segmento b"""
	A,B,C,D = a[0],a[1],b[0],b[1]
	if (in_segment(A,b) or in_segment(B,b) or in_segment(C,a) or in_segment(D,a)):
		return True
	return (sarea(A,B,C)*sarea(A,B,D))<0 and (sarea(C,D,A)*sarea(C,D,B)<0)

def sarea(A,B,C):
	"""Calcula el area signada del triangulo formado por los puntos ordenados A,B,C"""
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

def in_triangle(P,t):
    """Devuelve True si el punto P se encuentra dentro del triangulo t"""
    A,B,C=t[0],t[1],t[2]
    L=[orientation(A,B,P),orientation(B,C,P),orientation(C,A,P)]
    return not(1 in L and -1 in L)

def get_angles(a,b,c):
    get_angle = lambda a, b, c: degrees(acos((b*b+c*c-a*a)/(float(2*b*c))))
    if a + b <= c or b + c <= a or c + a <= b:
        return [0, 0, 0]
    return [get_angle(a, b, c), get_angle(b, c, a), get_angle(c, a, b)]

def random_poly_points(nVertex, nInteriorPoints):
	"""Devuelve [Aristas, puntosInteriores]"""
	p = [[np.random.rand(),np.random.rand()] for i in range(nVertex)]
	Poly = polygonization(p)
	p = []
	while len(p) < nInteriorPoints:
		point = [np.random.rand(),np.random.rand()]
		if pointInPolygon(point,Poly):
			p.append(point)
	return [[[Poly[i],Poly[(i+1)%len(Poly)]] for i in range(len(Poly))],Poly+p]

def in_circle(a,b,c,d):
	"""Devuelve true si el punto d se encuentra dentro del circulo que pasa por a,b,c"""
	sa=sarea(a,b,c)
	if sa==0:
		return
	return -np.sign(sa*svolume(a,b,c,d))

def svolume(a,b,c,d):
	"""Calcula el volumen signado del tetraedro formado por los puntos a,b,c,d"""
	arr = np.array([[1,1,1,1],
					[a[0],b[0],c[0],d[0]],
					[a[1],b[1],c[1],d[1]],
					[a[0]**2+a[1]**2,b[0]**2+b[1]**2, c[0]**2+c[1]**2,d[0]**2+d[1]**2]])
	return np.linalg.det(arr)/6

def polygonization(p):
	"""Devuelve la poligonizacion X-monotona de los puntos p"""
	U = []
	L = []
	maxim = max(p)
	minim = min(p)
	for i in range(len(p)):
		if sarea(maxim,p[i],minim)>=0:
			U.append(p[i])
		else:
			L.append(p[i])
	U = sorted(U)
	L = sorted(L)
	U.reverse()

	return L+U

def plot_poly_points(PolyPoints):
	lines = PolyPoints[0]
	lc = mc.LineCollection(lines, linewidths=2)
	fig, ax = plt.subplots()
	ax.add_collection(lc)
	x = [i[0] for i in PolyPoints[1]]
	y = [i[1] for i in PolyPoints[1]]
	plt.plot(x,y,'ro')
	x = [i[0][0] for i in lines]
	y = [i[0][1] for i in lines]
	plt.plot(x,y,'bo')
	plt.show()