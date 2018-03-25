import matplotlib.pyplot as plt
import numpy as np
import gtc
from matplotlib import collections  as mc

class Vertex:
	
	def __init__(self,coords,edge=None):
		self.coords = coords
		self.edge = edge
	
class Edge:
	def __init__(self,origin,twin=None,prev=None,
			  next_=None,face=None):
		
		self.origin = origin
		self.twin = twin
		self.prev = prev
		self.next = next_
		self.face = face
	
	def flip(self):
		
		#ARISTAS [origen, gemela, anterior, posterior, cara]
		ga = self.twin #gemela
		oga = ga.origin #origen gemela
		oa = self.origin #origen a
		ca = self.face #cara
		aa = self.prev #anterior
		pa = self.next #posterior
		cb = ga.face #cara de la gemela
		ab = ga.prev #anterior gemela
		pb = ga.next #posterior gemela
		
		self.origin = aa.origin
		self.prev = pa
		self.next = ab
		
		ga.origin = ab.origin
		ga.prev = pb
		ga.next = aa
		
		pa.prev = ab
		pa.next = self
		
		aa.prev = ga
		aa.next = pb
		aa.face = cb
		
		pb.prev = aa
		pb.next = ga
		
		ab.prev = self
		ab.next = pa
		ab.face = ca
		
		ca.edge = self
		cb.edge = ga
		oga.edge = pa
		oa.edge = pb
		return
	
	def isLegal(self,face0):
		if self.face == face0 or self.twin.face == face0:
			return True
		A = self.origin.coords
		B = self.twin.prev.origin.coords
		C = self.next.origin.coords
		D = self.next.next.origin.coords
		if -1 in [gtc.orientation(A,B,C),gtc.orientation(A,C,D),gtc.orientation(B,C,D),gtc.orientation(B,D,A)]:
			return True
		return gtc.inCircle(A,C,D,B)==-1
	
class Face:
	
	def __init__(self,edge):
		self.edge = edge
		
class Dcel:
	
	def __init__(self,orderedPoints):
		n = len(orderedPoints)
		self.vertices = [Vertex(orderedPoints[i]) for i in range(n)]
		self.edges = [Edge(self.vertices[i]) for i in range(n)]
		self.edges += [Edge(self.vertices[(i+1)%n]) for i in range(n)]
		
		for i in range(n):
			edge = self.edges[i]
			self.vertices[i].edge = edge
			edge.twin = self.edges[n+i]
			edge.prev = self.edges[(i-1)%n]
			edge.next = self.edges[(i+1)%n]
			
			edge = self.edges[n+i]
			edge.twin = self.edges[i]
			edge.prev = self.edges[n+(i+1)%n]
			edge.next = self.edges[n+(i-1)%n]
			""" END  Vertices """
		
		self.faces = [Face(self.edges[n]),Face(self.edges[0])]
		""" END  Faces """
		
		for i in range(n):
			self.edges[i].face = self.faces[1]
			self.edges[n+i].face = self.faces[0]
		""" END edges """
	
	@classmethod
	def deloneFromPoints(cls,points):
		P = gtc.angularSort(points,min(points))
		D = cls(P)
		D.triangulateInterior()
		D.triangulateExterior()
		D.plot()
		D.legalize()
		return D
	
	def plot(self, edges=None):
		lines =  [[edge.origin.coords, edge.twin.origin.coords] for edge in self.edges]
		lc = mc.LineCollection(lines, linewidths=2)
		fig, ax = plt.subplots()
		ax.add_collection(lc)
		x = [i[0][0] for i in lines]
		y = [i[0][1] for i in lines]
		if edges:
			ax.add_collection(edges)
		plt.plot(x,y,'ro')
	
	def plotWithVertexNumber(self):
		lines =  [[edge.origin.coords, edge.twin.origin.coords] for edge in self.edges]
		lc = mc.LineCollection(lines, linewidths=2)
		fig, ax = plt.subplots()
		ax.add_collection(lc)
		x = []
		y = []
		for i,vertex in enumerate(self.vertices):
			x.append(vertex.coords[0])
			y.append(vertex.coords[1])
			plt.text(x[-1],y[-1]+0.02,i)
		plt.plot(x,y,'ro')
	
	def splitFace(self,e1,e2):
		if e1.face != e2.face or e2.origin == e1.twin.origin or e1.origin == e2.twin.origin:
			return "no diagonal"
		
		newEdge = Edge(e1.origin, None, e1.prev, e2, e1.face)
		twinNewEdge = Edge(e2.origin, None, e2.prev, e1, None)
		
		newEdge.twin = twinNewEdge
		twinNewEdge.twin = newEdge
		
		newFace = Face(twinNewEdge)
		twinNewEdge.face = newFace
		
		k = e1.face
		
		self.edges.append(newEdge)
		self.edges.append(twinNewEdge)
		
		e1.prev.next = newEdge
		e1.prev = twinNewEdge
		e2.prev.next = twinNewEdge
		e2.prev = newEdge
		
		i = e1
		while i != twinNewEdge:
			i.face = newFace
			i = i.next
	
		k.edge = newEdge
	
		self.faces.append(newFace)
	
	def triangulateInterior(self):
		"""Triangula la cara interna del DCEL D"""
		ini = self.faces[1].edge
		mostLeftEdge = ini
		iterator = ini
		while True:
			if iterator.origin.coords[0] < mostLeftEdge.origin.coords[0]:
				mostLeftEdge = iterator
			iterator = iterator.next
			if iterator == ini:
				break
		iterator = mostLeftEdge
		while True:
			if iterator.face == self.faces[0]:
				break
			if gtc.orientation(iterator.origin.coords,iterator.next.origin.coords,iterator.next.next.origin.coords) == 1:
				self.splitFace(iterator,iterator.next.next)
				iterator = iterator.prev.twin
			else:
				iterator = iterator.next
	
	def triangulateExterior(self):
		ini = self.faces[0].edge
		mostLeftEdge = ini
		iterator = ini
		while True:
			if iterator.origin.coords[0] < mostLeftEdge.origin.coords[0]:
				mostLeftEdge = iterator
			iterator = iterator.next
			if iterator == ini:
				break
		iterator = mostLeftEdge
		while True:
			if gtc.orientation(iterator.origin.coords,iterator.next.origin.coords,iterator.next.next.origin.coords) == 1:
				self.splitFace(iterator,iterator.next.next)
				iterator = iterator.prev.twin.prev
				continue
			else:
				iterator = iterator.next
			if iterator.next.origin == mostLeftEdge.origin:
				break
			
	def legalize(self):
		flipped = True
		lastChanged = []
		while flipped:
			changed = []
			flipped = False
			for edge in self.edges:
				if not edge.isLegal(self.faces[0]):
					changed.append(edge)
					edge.flip()
					flipped = True
			if changed in lastChanged:
				return
			else:
				lastChanged.append(changed)
		return
	
points = [[np.random.rand(),np.random.rand()] for i in range(100)]
D = Dcel.deloneFromPoints(points)
D.plot()