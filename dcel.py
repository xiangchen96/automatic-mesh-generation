import matplotlib.pyplot as plt
import numpy as np
import utils
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
    
    def add_force_vector(self, coef=0.05):
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
    
    def is_legal(self,face0):
        if self.face == face0 or self.twin.face == face0:
            return True
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        if -1 in [utils.orientation(A,B,C),utils.orientation(A,C,D),
                  utils.orientation(B,C,D),utils.orientation(B,D,A)]:
            return True
        else:
            return utils.in_circle(A,C,D,B)==-1
    
    def is_flippable(self,face0):
        if self.face == face0 or self.twin.face == face0:
            return False
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        return -1 not in [utils.orientation(A,B,C),utils.orientation(A,C,D),
                          utils.orientation(B,C,D),utils.orientation(B,D,A)]
    
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
        coords_2 = self.twin.origin.coords
        x = (coords_1[0]+coords_2[0]) / 2
        y = (coords_1[1]+coords_2[1]) / 2
        return [x,y]
    
    
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
        self.min_x = None
        self.max_x = None
        self.min_y = None
        self.max_y = None
        self.vertices = []
        for point in orderedPoints:
            self.vertices.append(Vertex(point))
            if not self.min_x or point[0] < self.min_x:
                self.min_x = point[0]
            elif not self.max_x or point[0] > self.max_x:
                self.max_x = point[0]
            if not self.min_y or point[1] < self.min_y:
                self.min_y = point[1]
            elif not self.max_y or point[1] > self.max_y:
                self.max_y = point[1]
            
        self.edges = [Edge(self.vertices[i]) for i in range(n)]
        self.edges += [Edge(self.vertices[(i+1)%n]) for i in range(n)]
        self.polygon = None
        self.alpha = None
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
    def deloneFromPolygonFile(cls,fileName):
        f = open(fileName,"r")
        points = []
        number_of_points = int(f.readline())
        for i in range(number_of_points):
            line = f.readline()
            x,y = line.split(" ")
            points.append([float(x),float(y)])
        D = Dcel.deloneFromPoints(points)
        D.polygon = [[points[i],points[(i+1)%len(points)]] for i in range(len(points))]
        D.enforce_edges(D.polygon)
        f.close()
        return D
    
    @classmethod
    def deloneFromPoints(cls,points):
        P = utils.angular_sort(points,min(points))
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
            points = np.array([vertex.coords for vertex in self.vertices])
            simplices = []
            for face in self.get_interior_triangles(self.polygon):
                a,b,c = face.getVertices()
                simplices.append([self.vertices.index(a),self.vertices.index(b),self.vertices.index(c)])
            plt.triplot(points[:,0], points[:,1], simplices)
            plt.plot(points[:,0], points[:,1], 'bo')
            plt.show()
    
    def plot(self):
        plotted = []
        for edge in self.edges:
            if edge.twin not in plotted:
                plotted.append(edge)
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
            if utils.orientation(iterator.origin.coords,
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
            if utils.orientation(iterator.origin.coords,
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
                if utils.segment_crossing([Vi,Vj],
                                       [edge.origin.coords,edge.next.origin.coords]):
                    crossingEdges.append(edge)
            while len(crossingEdges) > 0:
                e = crossingEdges.pop()
                if not e.is_flippable(self.faces[0]):
                    crossingEdges.insert(0,e)
                else:
                    e.flip()
                    if utils.segment_crossing([Vi,Vj],
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
        for vertex in self.vertices:
            if vertex.coords not in polygon_vertices:
                vertex.add_force_vector()
        self.legalize()
        self.enforce_edges(self.polygon)
        
    
    def animate_main(self):
        fig = plt.figure()
        ax = plt.axes(xlim=(self.min_x-1,self.max_x+1), ylim=(self.min_y-1, self.max_y+1))
        angle_text = plt.text(self.max_x-2, self.max_y, '', fontsize=10)
        iteration = plt.text(self.max_x-2, self.max_y-1, '', fontsize=10)
        
        lines = [plt.plot([], [],'bo-')[0] for _ in range(len(self.vertices)**2)]
        def init():
            for line in lines:
                line.set_data([], [])
            return lines
        def animate(frame):
            if frame%10 == 0:
                self.add_point()
            else:
                self.iterate_forces()
            angle_text.set_text("min_angle: "+str(self.get_minimun_angle())+"º")
            iteration.set_text("iter: "+str(frame))
            edges = []
            for face in self.get_interior_triangles(self.polygon):
                for edge in face.get_edges():
                    if [edge.twin.origin.coords,edge.origin.coords] not in edges:
                        edges.append([edge.origin.coords, edge.twin.origin.coords])
            for i,edge in enumerate(edges):
                lines[i].set_data([edge[0][0],edge[1][0]],[edge[0][1],edge[1][1]])
            return lines+[angle_text,iteration]
        ani = animation.FuncAnimation(fig, animate, init_func=init,interval=0, blit=True)
#        ani.save('mover_legalizar.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
        plt.show()
        
    def get_interior_triangles(self, polygon):
        """ Returns [points, interiorSimplices] """
        poly = [i[0] for i in polygon]
        triangles = [face.getVertices() for face in self.faces[1:]]
        faces = []
        for i,(a,b,c) in enumerate(triangles):
             x = (a.coords[0]+b.coords[0]+c.coords[0])/3
             y = (a.coords[1]+b.coords[1]+c.coords[1])/3
             if utils.pointInPolygon([x,y],poly):
                 faces.append(self.faces[i+1])
        return faces
    
    def add_point(self):
        if self.polygon:
            new_point = None
            for face in self.get_interior_triangles(self.polygon):
                a,b,c = face.getVertices()
                a1 = np.linalg.norm([a.coords[0]-b.coords[0],
                  a.coords[1]-b.coords[1]])
                a2 = np.linalg.norm([a.coords[0]-c.coords[0],
                  a.coords[1]-c.coords[1]])
                a3 = np.linalg.norm([c.coords[0]-b.coords[0],
                  c.coords[1]-b.coords[1]])
                for angle in utils.get_angles(a1,a2,a3):
                    if angle < self.alpha:
#                        for edge in face.get_edges():
#                            if edge.twin.face == self.faces[0]:
#                                self.splitEdge(edge)
#                                break
#                        else:
                        x = (a.coords[0]+b.coords[0]+c.coords[0])/3
                        y = (a.coords[1]+b.coords[1]+c.coords[1])/3
                        new_point = [x,y]
                        break
                if new_point:
                    break
            if new_point == None:
                return
            puntos = [vertex.coords for vertex in self.vertices]+[new_point]
            D = Dcel.deloneFromPoints(puntos)
            D.polygon = self.polygon
            D.enforce_edges(D.polygon)
            self.vertices = D.vertices
            self.edges = D.edges
            self.faces = D.faces
            self.nuevos_vertices = []
    
    def get_minimun_angle(self):
        angles = []
        if self.polygon:
            for face in self.get_interior_triangles(self.polygon):
                a,b,c = face.getVertices()
                a1 = np.linalg.norm([a.coords[0]-b.coords[0],
                  a.coords[1]-b.coords[1]])
                a2 = np.linalg.norm([a.coords[0]-c.coords[0],
                  a.coords[1]-c.coords[1]])
                a3 = np.linalg.norm([c.coords[0]-b.coords[0],
                  c.coords[1]-b.coords[1]])
                angles += utils.get_angles(a1,a2,a3)
        else:
            for face in self.faces[1:]:
                a,b,c = (edge.get_length() for edge in face.get_edges())
                angles += utils.get_angles(a,b,c)
        return min(angles)
    
    def splitEdge(self,edge):
        new_point = edge.mid_point()
        a,b = edge.origin.coords, edge.destino().coords
        for i, edge in enumerate(self.polygon):
            if (edge[0] == a and edge[1] == b):
                aux = edge[1]
                edge[1] = new_point
                self.polygon.append([new_point,b])
                break
            if (edge[0] == b and edge[1] == a):
                aux = edge[1]
                edge[1] = new_point
                self.polygon.append([new_point,a])
                
        puntos = [vertex.coords for vertex in self.vertices]+[new_point]
        D = Dcel.deloneFromPoints(puntos)
        D.polygon = self.polygon
        D.enforce_edges(D.polygon)
        self.vertices = D.vertices
        self.edges = D.edges
        self.faces = D.faces
        self.nuevos_vertices = []
    
    def generate_mesh(self):
        iteration = 0
        while self.get_minimun_angle() < self.alpha and iteration < 1000:
            if iteration%10 == 0:
                self.add_point()
            else:
                self.iterate_forces()
            iteration += 1
        print("angulo ",self.get_minimun_angle())
        print("vertices ", len(self.vertices))
    
""" Meshing """
#D = Dcel.deloneFromPolygonFile("puntos")
#D.plotPolygon()
#D.alpha = 30
#D.generate_mesh()
#D.plotPolygon()

""" Demo Animation """
D = Dcel.deloneFromPolygonFile("puntos")
D.alpha = 28
D.animate_main()