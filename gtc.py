import numpy as np
from functools import cmp_to_key
import matplotlib.pyplot as plt
from matplotlib import collections  as mc

def midPoint(A,B):
	"""Calcula el punto medio entre los puntos A y B"""
	return [(A[0]+B[0])/2,(A[1]+B[1])/2]

def svolume(a,b,c,d):
	"""Calcula el volumen signado del tetraedro formado por los puntos a,b,c,d"""
	arr = np.array([[1,1,1,1],
					[a[0],b[0],c[0],d[0]],
					[a[1],b[1],c[1],d[1]],
					[a[0]**2+a[1]**2,b[0]**2+b[1]**2, c[0]**2+c[1]**2,d[0]**2+d[1]**2]])
	return np.linalg.det(arr)/6

def segmentIntersectionTest(a,b):
	"""Devuelve True si el segmento a se corta con el segmento b"""
	A,B,C,D = a[0],a[1],b[0],b[1]
	if (inSegment(A,b) or inSegment(B,b) or inSegment(C,a) or inSegment(D,a)):
		return True
	return (sarea(A,B,C)*sarea(A,B,D))<0 and (sarea(C,D,A)*sarea(C,D,B)<0)

def segmentCrossing(a,b):
	A,B,C,D = a[0],a[1],b[0],b[1]
	return (sarea(A,B,C)*sarea(A,B,D))<0 and (sarea(C,D,A)*sarea(C,D,B)<0)

def inSegment(P,s):
	"""Devuelve True si el punto P esta en el segmento s"""
	A,B = s[0],s[1]
	if sarea(P,A,B) != 0: return False
	if A[0] == B[0]: return min(A[1],B[1]) <= P[1] <= max(A[1],B[1])
	return min(A[0],B[0]) <= P[0] <= max(A[0],B[0])

def inCircle(a,b,c,d):
	"""Devuelve true si el punto d se encuentra dentro del circulo que pasa por a,b,c"""
	sa=sarea(a,b,c)
	if sa==0:
		return
	return -np.sign(sa*svolume(a,b,c,d))

def sarea(A,B,C):
	"""Calcula el area signada del triangulo formado por los puntos ordenados A,B,C"""
	return ((B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1]))/2

def dist2(A,B):
	"""Calcula el cuadrado de la distancia entre los puntos A y B"""
	return (A[0]-B[0])**2+(A[1]-B[1])**2

def orientation(A,B,C):
	"""Calcula la orientacion del triangulo formado por los puntos A,B,C"""
	return np.sign(((B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1])))

def next_(e,D):
	return D[1][e][3]

def angularSort(p,c):
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

def dcelFromPolygon(P):
	"""Crea un DCEL(Doubly connected edge list) a partir de un conjunto de
	puntos ordenados.

	ESTRUCTURA:
		Vertices = [Coords, arista]
		Aristas  = [origen, gemela, anterior, posterior, cara]
		Caras	= [arista]

		[Lvertices, Laristas, Lcaras]
	"""
	n=len(P)

	#VERTICES [coordenadas, arista_sale]
	V=[[P[i],i] for i in range(len(P))]

	#ARISTAS [origen, gemela, anterior, posterior, cara]
	e=[[i,n+i,(i-1)%n,(i+1)%n,1]for i in range(n)]+[[(i+1)%n,i,n+(i+1)%n,n+(i-1)%n,0]for i in range(n)]

	#CARAS [arista_cara0, arista_cara1]
	f=[n,0]

	return [V,e,f]

def origin(e,D):
	return D[1][e][0]

def originCoords(e,D):
	return D[0][D[1][e][0]][0]

def twin(e,D):
	return D[1][e][1]

def face(e,D):
	return D[1][e][4]

def faceEdges(c,D):
	f=D[2][c]   #la arista que le identifica
	C=[f]	   #lista de aristas
	f=next_(f,D) #meto las siguientes
	while f != C[0]:
		C.append(f)
		f=next_(f,D)
	return C

def faceVertices(c,D):
	aristas = faceEdges(c,D)
	return [origin(i,D) for i in aristas]

def prev(e,D):
	return D[1][e][2]

def isEdgeInD(edge,D):
	for i in range(len(D[1])):
		if originCoords(i,D) in edge and originCoords(next_(i,D),D) in edge:
			return True
	else:
		return False


def splitFace(e1,e2,D):

	# si no son aristas de la misma cara o si son adyacentes sus origenes no definen una diagonal
	if face(e1,D) != face(e2,D) or origin(e2,D) == origin(twin(e1,D),D) or origin(e1,D) == origin(twin(e2,D),D):
		return "no diagonal"

	ne, nf = len(D[1]), len(D[2])
	preve1 = prev(e1,D)
	preve2 = prev(e2,D)
	k=face(e1,D)

	# añadimos las aristas nuevas [origen, gemela, anterior, post, cara]
	D[1].append([origin(e1,D),ne+1,preve1,e2,k])
	D[1].append([origin(e2,D),ne,preve2,e1,nf])

	# modificamos aristas afectadas
	D[1][preve1][3]=ne
	D[1][e1][2]=ne+1
	D[1][preve2][3]=ne+1
	D[1][e2][2]=ne
	i=e1
	while i != ne+1:
		D[1][i][4]=nf
		i=D[1][i][3]

	#modificamos la cara afectada
	D[2][k]=ne

	# añadimos la nueva cara
	D[2].append(ne+1)

def triangulacionDCEL(D):
	"""Triangula la cara interna del DCEL D"""
	ini = D[2][1]
	j = ini
	mincoordX = originCoords(j,D)[0]
	k = j
	while origin(next_(j,D),D) != origin(ini,D):
		for i in range(len(D[0])):
			x = originCoords(i,D)[0]
			if x <= mincoordX:
				if (x < mincoordX):
					mincoordX = x
					ini = origin(i,D)
					k = i
				elif originCoords(i,D)[1] < originCoords(k,D)[1]:
					ini = origin(i,D)
					k = i
		j = D[1][j][3]
	noHeTerminado = True
	while (origin(next_(k,D),D) != ini) or noHeTerminado:
		if sarea(D[0][D[1][k][0]][0],originCoords(next_(k,D),D),originCoords(next_(next_(k,D),D),D))>0:
			if splitFace(k,next_(next_(k,D),D),D)=="no diagonal":
				break
			k = D[2][1]
			if origin(D[1][k][3],D) == ini:
				noHeTerminado = True
		else:
			k = D[1][k][3]
			noHeTerminado = False

def convexHullDCEL(D):
	"""Crea el cierre convexo del DCEL D aplicando Graham en la cara externa"""
	i = D[2][0]
	aristaDelHull = 0
	j = i
	minCoordX = originCoords(i,D)[0]
	while origin(next_(j,D),D) != origin(i,D):
		x = originCoords(j,D)[0]
		if x < minCoordX:
			minCoordX = x
			aristaDelHull = j
		j = next_(j,D)
	if (aristaDelHull != 0):
		j = aristaDelHull
	ini = origin(j,D)
	noHeTerminado = True
	while (origin(next_(j,D),D) != ini) or noHeTerminado:
		if sarea(originCoords(j,D),originCoords(next_(j,D),D),originCoords(next_(next_(j,D),D),D))>0:
			splitFace(j,next_(next_(j,D),D),D)
			j = prev(D[2][0],D)
			if origin(next_(j,D),D) == ini:
				noHeTerminado = True
		else:
			j = next_(j,D)
			noHeTerminado = False
	while sarea(originCoords(j,D),originCoords(next_(j,D),D),originCoords(next_(next_(j,D),D),D))>0:
			splitFace(j,next_(next_(j,D),D),D)
			j = prev(D[2][0],D)

def triangulation(p):
	"""Crea un DCEL triangulado a partir de un conjunto de puntos p"""
	P=angularSort(p,min(p))
	D=dcelFromPolygon(P)
	triangulacionDCEL(D)
	convexHullDCEL(D)
	return D

def flip(a,D):
	"""Cambia la diagonal a de D por su contraria"""
	#ARISTAS [origen, gemela, anterior, posterior, cara]
	ga=D[1][a][1] #gemela
	oga=D[1][ga][0] #origen gemela
	oa=D[1][a][0] #origen a
	ca=D[1][a][4] #cara
	aa=D[1][a][2] #anterior
	pa=D[1][a][3] #posterior
	cb=D[1][ga][4] #cara de la gemela
	ab=D[1][ga][2] #anterior gemela
	pb=D[1][ga][3] #posterior gemela
	D[1][a]=[D[1][aa][0],ga,pa,ab,ca]
	D[1][ga]=[D[1][ab][0],a,pb,aa,cb]#bien
	D[1][pa][2]=ab
	D[1][pa][3]=a
	D[1][aa][2]=ga
	D[1][aa][3]=pb
	D[1][aa][4]=cb
	D[1][pb][2]=aa
	D[1][pb][3]=ga
	D[1][ab][2]=a
	D[1][ab][3]=pa
	D[1][ab][4]=ca
	D[2][ca]=a #faltaba
	D[2][cb]=ga #faltaba
	D[0][oga][1]=pa #faltaba
	D[0][oa][1]=pb #faltaba
	return

def flipable(a,D):
	"""Devuelve True si la diagonal a de D tiene contraria"""
	if face(a,D) == 0 or face(twin(a,D),D) == 0:
		return False
	A = originCoords(a,D)
	B = originCoords(prev(twin(a,D),D),D)
	C = originCoords(next_(a,D),D)
	D = originCoords(next_(next_(a,D),D),D)
	return -1 not in [orientation(A,B,C),orientation(A,C,D),orientation(B,C,D),orientation(B,D,A)]

def legal(a,D):
	"""Devuelve True si la arista a cumple la condicion de Delaunay"""
	if face(a,D) == 0 or face(twin(a,D),D) == 0:
		return True
	A = originCoords(a,D)
	B = originCoords(prev(twin(a,D),D),D)
	C = originCoords(next_(a,D),D)
	D = originCoords(next_(next_(a,D),D),D)
	if -1 in [orientation(A,B,C),orientation(A,C,D),orientation(B,C,D),orientation(B,D,A)]:
		return True
	return inCircle(A,C,D,B)==-1

def legalize(T):
	"""Cambia todas las aristas de T por aristas legales"""
	heFlipao = True
	lastChanged = []
	while heFlipao:
		changed = []
		heFlipao = False
		for i in range(len(T[1])):
			if not legal(i,T):
				changed.append(i)
				flip(i,T)
				heFlipao = True
		if changed in lastChanged:
			return
		else:
			lastChanged.append(changed)
	return

def delone(p):
	"""Crea la triangulacion de Delaunay de un conjunto de puntos"""
	T = triangulation(p)
	legalize(T)
	return T

def plotDCEL(D,edges=None):
	nAristas = len(D[1])
	lines =  [[originCoords(i,D), originCoords(twin(i,D),D)] for i in range(nAristas)]
	lc = mc.LineCollection(lines, linewidths=2)
	fig, ax = plt.subplots()
	ax.add_collection(lc)
	x = [i[0][0] for i in lines]
	y = [i[0][1] for i in lines]
	if edges:
		ax.add_collection(edges)
	plt.plot(x,y,'ro')

def optimizer(p):
	delone(p)

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

def plotPolyPoints(PolyPoints):
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

def pointInPolygon(Q,p):
	"""Devuelve True si el punto Q esta dentro del poligono p"""
	count = 0
	disparo = [Q,[max(p)[0]+1,Q[1]]]
	n = len(p)
	for i in range(len(p)):
		if inSegment(Q,[p[i],p[(i+1)%n]]):
			return True
		if segmentIntersectionTest(disparo,[p[i],p[(i+1)%n]]):
			if inSegment(p[i],disparo):
				if sarea(Q,p[i],p[i-1])*sarea(Q,p[i],p[(i+1)%n]) < 0:
					count += 1
			else:
				count += 1
	return count%2==1

def randomPolyPoints(nVertex, nInteriorPoints):
	"""Devuelve [Aristas, puntosInteriores]"""
	p = [[np.random.rand(),np.random.rand()] for i in range(nVertex)]
	Poly = polygonization(p)
	p = []
	while len(p) < nInteriorPoints:
		point = [np.random.rand(),np.random.rand()]
		if pointInPolygon(point,Poly):
			p.append(point)
	return [[[Poly[i],Poly[(i+1)%len(Poly)]] for i in range(len(Poly))],Poly+p]

def constrainedDelaunay(puntos,aristas):
	D = delone(puntos)
	for (Vi,Vj) in aristas:
		if isEdgeInD([Vi,Vj],D):
			continue
		newEdges = []
		crossingEdges = []
		for i in range(len(D[1])):
			if twin(i,D) in crossingEdges:
				continue
			if segmentCrossing([Vi,Vj],[originCoords(i,D),originCoords(next_(i,D),D)]):
				crossingEdges.append(i)
		while len(crossingEdges) > 0:
			e = crossingEdges.pop()
			if not flipable(e,D):
				crossingEdges.insert(0,e)
			else:
				flip(e,D)
				if segmentCrossing([Vi,Vj],[originCoords(e,D),originCoords(next_(e,D),D)]):
					crossingEdges.insert(0,e)
				else:
					newEdges.append(e)
		swap = True
		while swap:
			swap = False
			for e in newEdges:
				if [originCoords(e,D),originCoords(next_(e,D),D)] == [Vi,Vj] or [originCoords(e,D),originCoords(next_(e,D),D)] == [Vj,Vi]:
					continue
				if not legal(e,D):
					flip(e,D)
					swap = True
	return D

def dcelFromDelaunay(tri):
	"""
	Crea un DCEL(Doubly connected edge list) a partir de un Delaunay.

	ESTRUCTURA:
	Vertices = [Coords, arista]
	Aristas  = [origen, gemela, anterior, posterior, cara]
	Caras	= [arista]
	"""
	##########rehacer con dictionary
	V = tri.points.copy()
	e = []
	f = []
	for i,(a,b,c) in enumerate(tri.simplices):
		n = len(e)
		A,B,C = 0,1,2
		ab,ba,bc,cb,ca,ac = n,n+1,n+2,n+3,n+4,n+5
		nuevaCara = i+1
		vecinos = list(map(lambda x: x+1,tri.neighbors[i]))
		#a-b
		e.append([a,ba,ca,bc,nuevaCara])
		#b-a
		e.append([b,ab,cb,ac,vecinos[C]])
		#b-c
		e.append([b,cb,ab,ca,nuevaCara])
		#c-b
		e.append([c,bc,ac,ba,vecinos[A]])
		#c-a
		e.append([c,ac,bc,ab,nuevaCara])
		#a-c
		e.append([a,ca,ba,cb,vecinos[B]])
		f.append(ab)
	for index,edge in enumerate(e):
		if edge[-1]==0:
			f.insert(0,index)
			break
	return [V,e,f]

def triangulatePolyPoints(polyP):
	D = constrainedDelaunay(polyP[1],polyP[0])
	polygon = [i[0] for i in polyP[0]]
	triangulos = [faceVertices(i,D) for i in range(1,len(D[2]))]
	borrar = []
	for i,(a,b,c) in enumerate(triangulos):
		 x = (D[0][a][0][0]+D[0][b][0][0]+D[0][c][0][0])/3
		 y = (D[0][a][0][1]+D[0][b][0][1]+D[0][c][0][1])/3
		 if not pointInPolygon([x,y],polygon):
			 borrar.append(i)
	triangulos = [t for i,t in enumerate(triangulos) if i not in borrar]
	return [np.array([D[0][i][0] for i in range(len(D[0]))]),triangulos]
