import matplotlib.pyplot as plt
import numpy as np
import gtc
from matplotlib import collections  as mc
from matplotlib import animation


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
    
    def add_force_vector(self,face0, coef=0.01):
        if face0 not in self.get_faces():
            vector = self._get_force_vector()
            self.coords[0] = self.coords[0] + vector[0]*coef
            self.coords[1] = self.coords[1] + vector[1]*coef
    
class Edge:
    """ 2-D Edge with an Origin Vertex, twin edge, previous edge and next edge """
    
    def __init__(self, origin, twin=None, prev=None, next_=None, face=None):
        self.origin = origin
        self.twin = twin
        self.prev = prev
        self.next = next_
        self.face = face
    
    def get_length(self):
        vector = [self.origin.coords[0]-self.twin.origin.coords[0],
                  self.origin.coords[1]-self.twin.origin.coords[1]]
        return np.linalg.norm(vector)
    
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
    
    def is_legal(self,face0):
        if self.face == face0 or self.twin.face == face0:
            return True
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        if -1 in [gtc.orientation(A,B,C),gtc.orientation(A,C,D),
                  gtc.orientation(B,C,D),gtc.orientation(B,D,A)]:
            return True
        else:
            return gtc.inCircle(A,C,D,B)==-1
    
    def isFlippable(self,face0):
        if self.face == face0 or self.twin.face == face0:
            return False
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        return -1 not in [gtc.orientation(A,B,C),gtc.orientation(A,C,D),
                          gtc.orientation(B,C,D),gtc.orientation(B,D,A)]
    
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
    
class Face:
    
    def __init__(self,edge):
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
        D.triangulate_interior()
        D.triangulate_exterior()
        D.legalize()
        return D
    
    @classmethod
    def triangulatePolygonWithPoints(cls,points,polygon):
        D = Dcel.deloneFromPoints(points)
        D.enforce_edges(polygon)
        D.polygon = polygon
        return D
    
    def plotPolygon(self):
        if self.polygon:
            points, simplices = D.get_interior_triangles(self.polygon)
            plt.triplot(points[:,0], points[:,1], simplices)
            plt.plot(points[:,0], points[:,1], 'bo')
            plt.show()
    
    def plot(self):
        for edge in self.edges:
            plt.plot([edge.origin.coords[0],edge.twin.origin.coords[0]],
                     [edge.origin.coords[1],edge.twin.origin.coords[1]],'bo-')
        plt.show()
    
    def plotWithEdges(self, edges):
        lines =  [[edge.origin.coords, edge.twin.origin.coords] for edge in self.edges]
        lc = mc.LineCollection(lines, linewidths=2)
        fig, ax = plt.subplots()
        ax.add_collection(lc)
        x = [i[0][0] for i in lines]
        y = [i[0][1] for i in lines]
        lc = mc.LineCollection(edges, linewidths=2,color='k')
        ax.add_collection(lc)
        plt.plot(x,y,'ro')
        plt.show()
    
    def plot_with_vertex_number(self):
        for edge in self.edges:
            plt.plot([edge.origin.coords[0],edge.twin.origin.coords[0]],
                     [edge.origin.coords[1],edge.twin.origin.coords[1]],'b')
        x,y = [],[]
        for i,vertex in enumerate(self.vertices):
            x.append(vertex.coords[0])
            y.append(vertex.coords[1])
            plt.text(x[-1],y[-1]+0.02,i)
        plt.plot(x,y,'ro')
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
    
    def triangulate_interior(self):
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
            if gtc.orientation(iterator.origin.coords,
                               iterator.next.origin.coords,
                               iterator.next.next.origin.coords) == 1:
                self.splitFace(iterator,iterator.next.next)
                iterator = iterator.prev.twin
            else:
                iterator = iterator.next
    
    def triangulate_exterior(self):
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
            if gtc.orientation(iterator.origin.coords,
                               iterator.next.origin.coords,
                               iterator.next.next.origin.coords) == 1:
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
                if not edge.is_legal(self.faces[0]):
                    changed.append(edge)
                    edge.flip()
                    flipped = True
            if changed in lastChanged:
                return
            else:
                lastChanged.append(changed)
        return
    
    def contains_edge(self,searched_edge):
        for edge in self.edges:
            if (edge.origin.coords in searched_edge and 
                edge.next.origin.coords in searched_edge):
                return True
        else:
            return False
    
    def enforce_edges(self,edges):
        for (Vi,Vj) in edges:
            if self.contains_edge([Vi,Vj]):
                continue
            newEdges = []
            crossingEdges = []
            for edge in self.edges:
                if edge.twin in crossingEdges:
                    continue
                if gtc.segmentCrossing([Vi,Vj],
                                       [edge.origin.coords,edge.next.origin.coords]):
                    crossingEdges.append(edge)
            while len(crossingEdges) > 0:
                e = crossingEdges.pop()
                if not e.isFlippable(self.faces[0]):
                    crossingEdges.insert(0,e)
                else:
                    e.flip()
                    if gtc.segmentCrossing([Vi,Vj],
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
    
    def animate_forces(self):
        fig = plt.figure()
        ax = plt.axes(xlim=(0, 1.3), ylim=(0, 1.3))
        angle_text = plt.text(1, 0.9, '', fontsize=10)
        iteration = plt.text(1, 0.8, '', fontsize=10)
        
        lines = [plt.plot([], [],'bo-')[0] for _ in range(len(self.edges))]
        def init():
            for line in lines:
                line.set_data([], [])
            return lines
        def animate(frame):
            for vertex in self.vertices:
                vertex.add_force_vector(self.faces[0])
            self.legalize()
            angle_text.set_text("min_angle: "+str(self.get_minimun_angle())+"º")
            iteration.set_text("iter: "+str(frame))
            edges =  [[edge.origin.coords, edge.twin.origin.coords] for edge in self.edges]
            for i,edge in enumerate(edges):
                lines[i].set_data([edge[0][0],edge[1][0]],[edge[0][1],edge[1][1]])
            return lines+[angle_text,iteration]
        ani = animation.FuncAnimation(fig, animate, init_func=init,interval=3, blit=True)
        plt.show()
#       ani.save('mover_legalizar.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    
    def iterate_forces(self):
        for vertex in self.vertices:
            vertex.add_force_vector(self.faces[0])
        self.legalize()
    
    def get_interior_triangles(self, polygon):
        """ Returns [points, interiorSimplices] """
        poly = [i[0] for i in polygon]
        triangles = [face.getVertices() for face in self.faces if face != self.faces[0]]
        exterior = []
        for i,(a,b,c) in enumerate(triangles):
             x = (a.coords[0]+b.coords[0]+c.coords[0])/3
             y = (a.coords[1]+b.coords[1]+c.coords[1])/3
             if not gtc.pointInPolygon([x,y],poly):
                 exterior.append(i)
        triangles = [t for i,t in enumerate(triangles) if i not in exterior]
        return [np.array([vertex.coords for vertex in self.vertices]),
          [[self.vertices.index(a),self.vertices.index(b),self.vertices.index(c)] for (a,b,c) in triangles]]
    
    def remove_edge(self, edge):
        if type(edge) is Edge:
            if edge.face == self.faces[0]:
                edge = edge.twin
            self.edges.remove(edge)
            self.edges.remove(edge.twin)
            if edge.twin.face != edge.face:
                self.faces.remove(edge.face)
            edge.remove()
            
        elif type(edge) is int:
            edge = self.edges[edge]
            if edge.face == self.faces[0]:
                edge = edge.twin
            self.edges.remove(edge)
            self.edges.remove(edge.twin)
            if edge.twin.face != edge.face:
                self.faces.remove(edge.face)
            edge.remove()
        
    def remove_vertex(self, vertex):
        face = None
        if type(vertex) == Vertex:
            self.vertices.remove(vertex)
            for edge in vertex.get_edges():
                face = edge.next.face
                self.remove_edge(edge)
        if type(vertex) == int:
            vertex = self.vertices[vertex]
            self.vertices.remove(vertex)
            for edge in vertex.get_edges():
                face = edge.next.face
                self.remove_edge(edge)
        return face
    
    def triangulate_face(self, face):
        vertices = [vertex.coords for vertex in face.getVertices()]
        if len(vertices) == 3: return
        iterator = face.edge
        while True:
            if iterator.next.next.next == iterator:
                break
            A = iterator.origin.coords
            B = iterator.next.origin.coords
            C = iterator.next.next.origin.coords
            if gtc.sarea(A,B,C) < 0:
                iterator = iterator.next
                continue
            for v in vertices:
                if v in [A,B,C]: continue
                if gtc.inTriangle(v,[A,B,C]): 
                    iterator = iterator.next
                    break
            else:
                self.splitFace(iterator,iterator.next.next)
                iterator = iterator.prev.twin
                continue
    
    def get_minimun_angle(self):
        angles = []
        for face in self.faces[1:]:
            a,b,c = (edge.get_length() for edge in face.get_edges())
            angles += gtc.get_angles(a,b,c)
        return min(angles)
            
""" NORMAL DELONE """
points = [list(np.random.uniform(0,1,2)) for i in range(7)]
D = Dcel.deloneFromPoints(points)
D.plot_with_vertex_number()

""" ANIMATION DELONE """
#points = [list(np.random.uniform(0,1,2)) for i in range(30)]
#D = Dcel.deloneFromPoints(points)
#D.animate_forces()

""" Polygon """
#polyP = gtc.randomPolyPoints(20,5)
#gtc.plotPolyPoints(polyP)
#D = Dcel.triangulatePolygonWithPoints(polyP[1],polyP[0])
#D.plotPolygon()