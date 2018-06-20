import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
import numpy as np
from matplotlib import animation
import collections
from scipy.spatial import Delaunay
from math import acos, degrees

#FUNCIONES AUXILIARES
def sarea(A,B,C):
    """Calcula el area signada del triangulo formado por los puntos ordenados A,B,C"""
    return ((B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1]))/2

def segment_crossing(a,b):
    """Devuelve True si el segmento a interseca con el segmento b"""
    A,B,C,D = a[0],a[1],b[0],b[1]
    return (sarea(A,B,C)*sarea(A,B,D))<0 and (sarea(C,D,A)*sarea(C,D,B)<0)

def pointInPolygon(Q,p):
    """Devuelve True si el punto Q esta dentro del poligono p"""
    path = mpltPath.Path(p)
    inside = path.contains_points([Q])
    return inside

def get_angles(a,b,c):
    """Devuelve una lista de los ángulos del triángulo a,b,c"""
    get_angle = lambda a, b, c: degrees(acos((b*b+c*c-a*a)/(float(2*b*c))))
    if a + b <= c or b + c <= a or c + a <= b:
        return [0, 0, 0]
    return [get_angle(a, b, c), get_angle(b, c, a), get_angle(c, a, b)]

def orientation(A,B,C):
	"""Calcula la orientacion del triangulo formado por los puntos A,B,C"""
	return np.sign(((B[0]-A[0])*(C[1]-A[1])-(C[0]-A[0])*(B[1]-A[1])))

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

class Vertex:
    """ 2-D Vertex with coordinates and an Edge """
    
    def __init__(self,coords,edge=None):
        self.coords = coords
        self.edge = edge
    
    def get_edges(self):
        e0 = self.edge
        edges = [e0]
        e = e0.prev.twin
        while e != e0:
            edges.append(e)
            e = e.prev.twin
        return edges
    
    def get_neighbours(self):
        edges = self.get_edges()
        return [edge.next.origin for edge in edges]
    
    def get_faces(self):
        edges = self.get_edges()
        return [edge.face for edge in edges]
    
    def _get_force_vector(self):
        vectors = []
        for vertex in self.get_neighbours():
            x,y = vertex.coords
            my_x, my_y = self.coords
            vectors.append([x-my_x, y-my_y])
        vector = [sum(v[0] for v in vectors),sum(v[1] for v in vectors)]
        return vector
    
    def add_force_vector(self, point_projection=None, coef=0.05):
        force = self._get_force_vector()
        if point_projection:
            vector_project = [point_projection[0]-self.coords[0],
                          point_projection[1]-self.coords[1]]
            aux = (vector_project[0]*force[0] + vector_project[1]*force[1]) \
                  /(vector_project[0]**2 + vector_project[1]**2)
            force = [aux*vector_project[0],aux*vector_project[1]]
        self.coords[0] = self.coords[0] + force[0]*coef
        self.coords[1] = self.coords[1] + force[1]*coef
    
class Edge:
    """ 2-D Edge with an Origin Vertex, twin edge, previous edge and next edge """
    
    def __init__(self, origin, twin=None, prev=None, next_=None, face=None):
        self.origin = origin
        self.twin = twin
        self.prev = prev
        self.next = next_
        self.face = face
    
    def get_length(self):
        vector = [self.origin.coords[0]-self.next.origin.coords[0],
                  self.origin.coords[1]-self.next.origin.coords[1]]
        return np.linalg.norm(vector)
    
    def destino(self):
        return self.next.origin
    
    def flip(self):
        
        twin = self.twin
        origin = self.origin
        face = self.face
        prev = self.prev
        next_ = self.next
        
        origin_twin = twin.origin
        face_twin = twin.face
        prev_twin = twin.prev
        next_twin = twin.next
        
        self.origin = prev.origin
        self.prev = next_
        self.next = prev_twin
        
        twin.origin = prev_twin.origin
        twin.prev = next_twin
        twin.next = prev
        
        next_.prev = prev_twin
        next_.next = self
        
        prev.prev = twin
        prev.next = next_twin
        prev.face = face_twin
        
        next_twin.prev = prev
        next_twin.next = twin
        
        prev_twin.prev = self
        prev_twin.next = next_
        prev_twin.face = face
        
        face.edge = self
        face_twin.edge = twin
        origin_twin.edge = next_
        origin.edge = next_twin
        return
    
    def is_flippable(self,face0):
        if self.face == face0 or self.twin.face == face0:
            return False
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        return -1 not in np.sign([sarea(A,B,C),sarea(A,C,D),
                                  sarea(B,C,D),sarea(B,D,A)])
    
    def remove(self):
        twin = self.twin
        prev = self.prev
        prev.next = twin.next
        prev.face = twin.face
        
        next_ = self.next
        next_.prev = twin.prev
        next_.face = twin.face
        
        twin.next.prev = self.prev
        twin.prev.next = self.next
        twin.face.edge = twin.prev
        
        self.origin.edge = twin.next
        twin.origin.edge = self.next
    
    def mid_point(self):
        coords_1 = self.origin.coords
        coords_2 = self.next.origin.coords
        x = (coords_1[0]+coords_2[0]) / 2
        y = (coords_1[1]+coords_2[1]) / 2
        return [x,y]
    
    def is_legal(self,face0):
        if self.face == face0 or self.twin.face == face0:
            return True
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        if -1 in [orientation(A,B,C),orientation(A,C,D),
                  orientation(B,C,D),orientation(B,D,A)]:
            return True
        else:
            return in_circle(A,C,D,B)==-1
    
class Face:
    
    def __init__(self,edge=None):
        self.edge = edge
    
    def get_edges(self):
        edge = self.edge
        edges= [edge]
        edge = edge.next
        while edge != edges[0]:
            edges.append(edge)
            edge = edge.next
        return edges
    
    def getVertices(self):
        edges = self.get_edges()
        return [edge.origin for edge in edges]
        
class Dcel:
    
    def __init__(self,points):
        tesselation = Delaunay(points)
        self.alpha = None
        self.min_x = None
        self.max_x = None
        self.min_y = None
        self.max_y = None
        edges = collections.defaultdict()
        self.vertices = []
        self.faces = [Face()]
        self.edges = []
        self.splitted = []
        for point in points:
            self.vertices.append(Vertex(point))
            if not self.min_x or point[0] < self.min_x: self.min_x = point[0]
            if not self.max_x or point[0] > self.max_x: self.max_x = point[0]
            if not self.min_y or point[1] < self.min_y: self.min_y = point[1]
            if not self.max_y or point[1] > self.max_y: self.max_y = point[1]
        for a,b,c in tesselation.simplices:
            edges[(a,b)] = Edge(self.vertices[a])
            edges[(b,c)] = Edge(self.vertices[b])
            edges[(c,a)] = Edge(self.vertices[c])
            self.edges.append(edges[(a,b)])
            self.edges.append(edges[(b,c)])
            self.edges.append(edges[(c,a)])
            
            self.vertices[a].edge = edges[(a,b)]
            self.vertices[b].edge = edges[(b,c)]
            self.vertices[c].edge = edges[(c,a)]
            
            face = Face(edges[(a,b)])
            self.faces.append(face)
            
            edges[(a,b)].face = face
            edges[(b,c)].face = face
            edges[(c,a)].face = face
            
            edges[(a,b)].prev = edges[(c,a)]
            edges[(a,b)].next = edges[(b,c)]
            
            edges[(b,c)].prev = edges[(a,b)]
            edges[(b,c)].next = edges[(c,a)]
            
            edges[(c,a)].prev = edges[(b,c)]
            edges[(c,a)].next = edges[(a,b)]
            
            if (b,a) in edges:
                edges[(a,b)].twin = edges[(b,a)]
                edges[(b,a)].twin = edges[(a,b)]
            if (c,b) in edges:
                edges[(b,c)].twin = edges[(c,b)]
                edges[(c,b)].twin = edges[(b,c)]
            if (a,c) in edges:
                edges[(a,c)].twin = edges[(c,a)]
                edges[(c,a)].twin = edges[(a,c)]
        hull = []
        for a,b in tesselation.convex_hull:
            if (a,b) not in edges:
                edges[(a,b)] = Edge(self.vertices[a],twin=edges[b,a],face=self.faces[0])
                edges[(b,a)].twin = edges[(a,b)]
                hull.append((a,b))
            elif (b,a) not in edges:
                edges[(b,a)] = Edge(self.vertices[a],twin=edges[a,b],face=self.faces[0])
                edges[(a,b)].twin = edges[(b,a)]
                hull.append((b,a))
        for a,b in hull:
            for c,d in hull:
                if c == b:
                    edges[a,b].next = edges[c,d]
                    edges[c,d].prev = edges[a,b]
                    break
                
                
    
    @classmethod
    def deloneFromPolygonFile(cls,fileName):
        f = open(fileName,"r")
        points = []
        number_of_points = int(f.readline())
        for i in range(number_of_points):
            line = f.readline()
            x,y = line.split(" ")
            points.append([float(x),float(y)])
        D = cls(points)
        D.polygon = [[points[i],points[(i+1)%len(points)]] for i in range(len(points))]
        D.enforce_edges()
        f.close()
        return D
    
    def plotPolygon(self):
        if self.polygon:
            plt.axes().set_aspect('equal')
            points = np.array([vertex.coords for vertex in self.vertices])
            simplices = []
            for face in self.get_interior_triangles():
                a,b,c = face.getVertices()
                simplices.append([self.vertices.index(a),self.vertices.index(b),self.vertices.index(c)])
            plt.triplot(points[:,0], points[:,1], simplices, 'o-')
            plt.show()
    
    def plot(self):
        plotted = []
        plt.axes().set_aspect('equal')
        for edge in self.edges:
            if edge.twin not in plotted:
                plotted.append(edge)
                plt.plot([edge.origin.coords[0],edge.next.origin.coords[0]],
                         [edge.origin.coords[1],edge.next.origin.coords[1]],'bo-')
        plt.show()
        
    def splitFace(self,e1,e2):
        if (e1.face != e2.face or 
            e2.origin == e1.twin.origin or 
            e1.origin == e2.twin.origin):
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
    
    def contains_edge(self,searched_edge):
        for edge in self.edges:
            if (edge.origin.coords in searched_edge and 
                edge.next.origin.coords in searched_edge):
                return True
        return False
    
    def enforce_edges(self):
        for (Vi,Vj) in self.polygon:
            if self.contains_edge([Vi,Vj]):
                continue
            newEdges = []
            crossingEdges = []
            for edge in self.edges:
                if edge.twin in crossingEdges:
                    continue
                if segment_crossing([Vi,Vj],
                                    [edge.origin.coords,edge.next.origin.coords]):
                    crossingEdges.append(edge)
            while len(crossingEdges) > 0:
                e = crossingEdges.pop()
                if not e.is_flippable(self.faces[0]):
                    crossingEdges.insert(0,e)
                else:
                    e.flip()
                    if segment_crossing([Vi,Vj],
                                           [e.origin.coords,e.next.origin.coords]):
                        crossingEdges.insert(0,e)
                    else:
                        newEdges.append(e)
            swap = True
            while swap:
                swap = False
                for e in newEdges:
                    if (e.origin.coords in [Vi,Vj] and 
                        e.next.origin.coords in [Vi,Vj]):
                        continue
                    if not e.is_legal:
                        e.flip()
                        swap = True
    
    def iterate_forces(self):
        polygon_vertices = [arista[0] for arista in self.polygon]
        for i,vertex in enumerate(self.vertices):
            if i in self.splitted:
                for a,b in self.polygon:
                    if a == vertex.coords:
                        vertex.add_force_vector(b)
                        break
                    if b == vertex.coords:
                        vertex.add_force_vector(a)
                        break
            elif vertex.coords not in polygon_vertices:
                vertex.add_force_vector()
        D = Dcel([vertex.coords for vertex in self.vertices])
        self.vertices = D.vertices
        self.edges = D.edges
        self.faces = D.faces
        self.enforce_edges()
        
    
    def animate_main(self):
        fig = plt.figure()
        ax = plt.axes(xlim=(self.min_x-1,self.max_x+1), ylim=(self.min_y-1, self.max_y+1))
        angle_text = plt.text(self.max_x/2, self.max_y, '', fontsize=10)
        iteration = plt.text(self.max_x/2, self.max_y-1, '', fontsize=10)
        
        lines = [plt.plot([], [],'bo-')[0] for _ in range(len(self.edges))]
        def init():
            for line in lines:
                line.set_data([], [])
            return lines
        def animate(frame):
            if frame%5 == 0:
                self.add_point()
            else:
                self.iterate_forces()
            angle = self.get_minimun_angle()
            if angle >= self.alpha:
                ani.event_source.stop()
            angle_text.set_text("min_angle: "+str(int(angle)))
            iteration.set_text("iter: "+str(frame))
            edges = []
            for face in self.get_interior_triangles():
                for edge in face.get_edges():
                    if [edge.next.origin.coords,edge.origin.coords] not in edges:
                        edges.append([edge.origin.coords, edge.next.origin.coords])
            for i,edge in enumerate(edges):
                if (len(lines) > i):
                    lines[i].set_data([edge[0][0],edge[1][0]],[edge[0][1],edge[1][1]])
                else:
                    lines.append(plt.plot([edge[0][0],edge[1][0]],[edge[0][1],edge[1][1]],'bo-')[0])
            return lines+[angle_text,iteration]
        ani = animation.FuncAnimation(fig, animate, init_func=init,interval=10, blit=True)
        # ani.save('main.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
        plt.show()
    
    def get_interior_triangles(self):
        """ Returns [points, interiorSimplices] """
        poly = [i[0] for i in self.polygon]
        triangles = [face.getVertices() for face in self.faces[1:]]
        faces = []
        for i,(a,b,c) in enumerate(triangles):
             x = (a.coords[0]+b.coords[0]+c.coords[0])/3
             y = (a.coords[1]+b.coords[1]+c.coords[1])/3
             if pointInPolygon([x,y],poly):
                 faces.append(self.faces[i+1])
        return faces
    
    def isConstrained(self,edge):
        if [edge.origin.coords,edge.destino().coords] in self.polygon:
            return True
        if [edge.destino().coords,edge.origin.coords] in self.polygon:
            return True
        return False
    
    def add_point(self):
        new_point = None
        face = self.get_face_with_min_angle()
        a,b,c = face.getVertices()
        a1 = np.linalg.norm([a.coords[0]-b.coords[0],
          a.coords[1]-b.coords[1]])
        a2 = np.linalg.norm([a.coords[0]-c.coords[0],
          a.coords[1]-c.coords[1]])
        a3 = np.linalg.norm([c.coords[0]-b.coords[0],
          c.coords[1]-b.coords[1]])
        for angle in get_angles(a1,a2,a3):
            if angle < self.alpha:
                for edge in face.get_edges():
                    if self.isConstrained(edge):
                        if edge.get_length() > 2*edge.next.get_length() \
                            or edge.get_length() > 2*edge.prev.get_length() \
                            or (edge.get_length() > edge.next.get_length() and \
                                edge.get_length() > edge.prev.get_length()):
                            self.splitEdge(edge)
                            return
                else:
                    x = (a.coords[0]+b.coords[0]+c.coords[0])/3
                    y = (a.coords[1]+b.coords[1]+c.coords[1])/3
                    new_point = [x,y]
        if new_point == None:
            return
        puntos = [vertex.coords for vertex in self.vertices]+[new_point]
        D = Dcel(puntos)
        D.polygon = self.polygon
        self.vertices = D.vertices
        self.edges = D.edges
        self.faces = D.faces
            
    def get_minimun_angle(self):
        angles = []
        if self.polygon:
            for face in self.get_interior_triangles():
                a,b,c = face.getVertices()
                a1 = np.linalg.norm([a.coords[0]-b.coords[0],
                  a.coords[1]-b.coords[1]])
                a2 = np.linalg.norm([a.coords[0]-c.coords[0],
                  a.coords[1]-c.coords[1]])
                a3 = np.linalg.norm([c.coords[0]-b.coords[0],
                  c.coords[1]-b.coords[1]])
                angles += get_angles(a1,a2,a3)
        else:
            for face in self.faces[1:]:
                a,b,c = (edge.get_length() for edge in face.get_edges())
                angles += get_angles(a,b,c)
        return min(angles)
    
    def splitEdge(self,split):
        new_point = split.mid_point()
        a,b = split.origin.coords, split.destino().coords
        for i, edge in enumerate(self.polygon):
            if (edge[0] == a and edge[1] == b):
                edge[1] = new_point
                self.polygon.insert(i+1,[new_point,b])
                break
            if (edge[0] == b and edge[1] == a):
                edge[1] = new_point
                self.polygon.insert(i+1,[new_point,a])
                
        puntos = [vertex.coords for vertex in self.vertices]+[new_point]
        self.splitted.append(len(puntos)-1)
        D = Dcel(puntos)
        D.polygon = self.polygon
        self.vertices = D.vertices
        self.edges = D.edges
        self.faces = D.faces
    
    def get_face_with_min_angle(self):
        min_face = None
        min_angle = 1000
        for face in self.get_interior_triangles():
            a,b,c = face.getVertices()
            a1 = np.linalg.norm([a.coords[0]-b.coords[0],
              a.coords[1]-b.coords[1]])
            a2 = np.linalg.norm([a.coords[0]-c.coords[0],
              a.coords[1]-c.coords[1]])
            a3 = np.linalg.norm([c.coords[0]-b.coords[0],
              c.coords[1]-b.coords[1]])
            minimo = min(get_angles(a1,a2,a3))
            if minimo < min_angle:
                min_face = face
                min_angle = minimo
        return min_face
    
    def generate_mesh(self, max_iterations=500):
        iteration = 0
        while self.get_minimun_angle() < self.alpha and iteration < max_iterations:
            if iteration%10 == 0:
                self.add_point()
            else:
                self.iterate_forces()
            iteration += 1