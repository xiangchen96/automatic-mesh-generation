import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
import collections
from scipy.spatial import Delaunay
from utilities import *


class Vertex:
    """Nodo vertice de un DCEL"""

    def __init__(self, coords, edge=None):
        self.coords = coords
        self.edge = edge

    @property
    def edges(self):
        yield self.edge
        e = self.edge.prev.twin
        while e != self.edge:
            yield e
            e = e.prev.twin

    @property
    def neighbours(self):
        for edge in self.edges:
            yield edge.next.origin

    @property
    def force_vector(self):
        """Vector suma de vertices adyacentes"""
        vector = [0, 0]
        for vertex in self.neighbours:
            x, y = vertex.coords
            my_x, my_y = self.coords
            vector[0] += x-my_x
            vector[1] += y-my_y
        return vector

    def add_force_vector(self, point=None, coef=0.05):
        """Mueve el punto en la direccion del vector suma (proyectado)"""
        force = self.force_vector

        if point:
            vector_project = [point[0]-self.coords[0], point[1]-self.coords[1]]
            force = project_vector(force, vector_project)

        self.coords[0] = self.coords[0] + force[0]*coef
        self.coords[1] = self.coords[1] + force[1]*coef


class Edge:
    """Semiarista de un DCEL"""

    def __init__(self, origin, twin=None, prev=None, next_=None, face=None):
        self.origin = origin
        self.twin = twin
        self.prev = prev
        self.next = next_
        self.face = face

    @property
    def length(self):
        vector = [self.origin.coords[0]-self.next.origin.coords[0],
                  self.origin.coords[1]-self.next.origin.coords[1]]
        return np.linalg.norm(vector)

    @property
    def destination(self):
        return self.next.origin

    def flip(self):
        """
        Gira la arista, la cambia por la otra diagonal del cuadrilatero que
        forman sus triangulos
        """
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

    def is_flippable(self, face0):
        if self.face == face0 or self.twin.face == face0:
            return False
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        return -1 not in np.sign([sarea(A, B, C), sarea(A, C, D),
                                  sarea(B, C, D), sarea(B, D, A)])

    def mid_point(self):
        coords_1 = self.origin.coords
        coords_2 = self.next.origin.coords
        x = (coords_1[0]+coords_2[0]) / 2
        y = (coords_1[1]+coords_2[1]) / 2
        return [x, y]

    def is_legal(self, face0):
        if self.face == face0 or self.twin.face == face0:
            return True
        A = self.origin.coords
        B = self.twin.prev.origin.coords
        C = self.next.origin.coords
        D = self.next.next.origin.coords
        if -1 in np.sign([sarea(A, B, C), sarea(A, C, D),
                          sarea(B, C, D), sarea(B, D, A)]):
            return True
        else:
            return in_circle(A, C, D, B) == -1


class Face:

    def __init__(self, edge=None):
        self.edge = edge

    @property
    def edges(self):
        yield self.edge
        edge = self.edge.next
        while edge != self.edge:
            yield edge
            edge = edge.next

    @property
    def vertices(self):
        for edge in self.edges:
            yield edge.origin


class Dcel:

    def __init__(self, points):
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
            if not self.min_x or point[0] < self.min_x:
                self.min_x = point[0]
            if not self.max_x or point[0] > self.max_x:
                self.max_x = point[0]
            if not self.min_y or point[1] < self.min_y:
                self.min_y = point[1]
            if not self.max_y or point[1] > self.max_y:
                self.max_y = point[1]
        for a, b, c in tesselation.simplices:
            edges[(a, b)] = Edge(self.vertices[a])
            edges[(b, c)] = Edge(self.vertices[b])
            edges[(c, a)] = Edge(self.vertices[c])
            self.edges.append(edges[(a, b)])
            self.edges.append(edges[(b, c)])
            self.edges.append(edges[(c, a)])

            self.vertices[a].edge = edges[(a, b)]
            self.vertices[b].edge = edges[(b, c)]
            self.vertices[c].edge = edges[(c, a)]

            face = Face(edges[(a, b)])
            self.faces.append(face)

            edges[(a, b)].face = face
            edges[(b, c)].face = face
            edges[(c, a)].face = face

            edges[(a, b)].prev = edges[(c, a)]
            edges[(a, b)].next = edges[(b, c)]

            edges[(b, c)].prev = edges[(a, b)]
            edges[(b, c)].next = edges[(c, a)]

            edges[(c, a)].prev = edges[(b, c)]
            edges[(c, a)].next = edges[(a, b)]

            if (b, a) in edges:
                edges[(a, b)].twin = edges[(b, a)]
                edges[(b, a)].twin = edges[(a, b)]
            if (c, b) in edges:
                edges[(b, c)].twin = edges[(c, b)]
                edges[(c, b)].twin = edges[(b, c)]
            if (a, c) in edges:
                edges[(a, c)].twin = edges[(c, a)]
                edges[(c, a)].twin = edges[(a, c)]
        hull = []
        for a, b in tesselation.convex_hull:
            if (a, b) not in edges:
                edges[(a, b)] = Edge(self.vertices[a],
                                     twin=edges[b, a],
                                     face=self.faces[0])
                edges[(b, a)].twin = edges[(a, b)]
                hull.append((a, b))
            elif (b, a) not in edges:
                edges[(b, a)] = Edge(self.vertices[a],
                                     twin=edges[a, b],
                                     face=self.faces[0])
                edges[(a, b)].twin = edges[(b, a)]
                hull.append((b, a))
        for a, b in hull:
            for c, d in hull:
                if c == b:
                    edges[a, b].next = edges[c, d]
                    edges[c, d].prev = edges[a, b]
                    break

    @classmethod
    def deloneFromPolygonFile(cls, fileName):
        f = open(fileName, 'r')
        points = []
        number_of_points = int(f.readline())
        for i in range(number_of_points):
            line = f.readline()
            x, y = line.split(' ')
            points.append([float(x), float(y)])
        D = cls(points)
        D.polygon = [[points[i], points[(i+1) % len(points)]] for i in range(len(points))]
        D.enforce_edges()
        f.close()
        return D

    @classmethod
    def deloneFromPolygon(cls, points):
        D = cls(points)
        D.polygon = [[points[i], points[(i+1) % len(points)]] for i in range(len(points))]
        D.enforce_edges()
        return D

    def plotPolygon(self):
        plt.axes().set_aspect('equal')
        for face in self.get_interior_triangles():
            a, b, c = face.vertices
            plt.triplot([a.coords[0], b.coords[0], c.coords[0]],
                        [a.coords[1], b.coords[1], c.coords[1]], 'bo-')
        plt.show()

    def plot(self):
        plt.axes().set_aspect('equal')
        for face in self.faces[1:]:
            a, b, c = face.vertices
            plt.triplot([a.coords[0], b.coords[0], c.coords[0]],
                        [a.coords[1], b.coords[1], c.coords[1]], 'bo-')
        plt.show()

    def contains_edge(self, searched_edge):
        for edge in self.edges:
            if [edge.origin.coords, edge.next.origin.coords] == searched_edge:
                return True
            if [edge.next.origin.coords, edge.origin.coords] == searched_edge:
                return True
        return False

    def enforce_edges(self):
        for (Vi, Vj) in self.polygon:
            if self.contains_edge([Vi, Vj]):
                continue
            newEdges = []
            crossingEdges = []
            for edge in self.edges:
                if edge.twin in crossingEdges:
                    continue
                if segment_crossing([Vi, Vj],
                                    [edge.origin.coords, edge.next.origin.coords]):
                    crossingEdges.append(edge)
            while len(crossingEdges) > 0:
                e = crossingEdges.pop()
                if not e.is_flippable(self.faces[0]):
                    crossingEdges.insert(0, e)
                else:
                    e.flip()
                    if segment_crossing([Vi, Vj],
                                        [e.origin.coords, e.next.origin.coords]):
                        crossingEdges.insert(0, e)
                    else:
                        newEdges.append(e)
            swap = True
            while swap:
                swap = False
                for e in newEdges:
                    if (e.origin.coords in [Vi, Vj] and
                       e.next.origin.coords in [Vi, Vj]):
                        continue
                    if not e.is_legal:
                        e.flip()
                        swap = True

    def iterate_forces(self):
        polygon_vertices = [arista[0] for arista in self.polygon]
        for i, vertex in enumerate(self.vertices):
            if i in self.splitted:
                for a, b in self.polygon:
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
        plt.axes(xlim=(self.min_x-1, self.max_x+1), ylim=(self.min_y-1, self.max_y+1))
        angle_text = plt.text((self.max_x+self.min_x)/2, self.max_y, '', fontsize=10)

        lines = [plt.plot([], [], 'bo-')[0] for _ in range(len(self.edges))]

        def init():
            for line in lines:
                line.set_data([], [])
            return lines

        def animate(frame):
            if frame % 5 == 0:
                self.add_point()
            else:
                self.iterate_forces()
            angle = self.minimum_angle
            if angle >= self.alpha:
                ani.event_source.stop()
            angle_text.set_text(f'min_angle: {angle:.2f} iter: {frame}')
            edges = []
            for face in self.get_interior_triangles():
                for edge in face.edges:
                    if [edge.next.origin.coords, edge.origin.coords] not in edges:
                        edges.append([edge.origin.coords, edge.next.origin.coords])
            for i, edge in enumerate(edges):
                if (len(lines) > i):
                    lines[i].set_data([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]])
                else:
                    lines.append(plt.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], 'bo-')[0])
            return lines+[angle_text]
        ani = animation.FuncAnimation(fig, animate, init_func=init, interval=10, blit=True)
        # ani.save('main.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
        plt.show()

    def get_interior_triangles(self):
        poly = [i[0] for i in self.polygon]
        triangles = [face.vertices for face in self.faces[1:]]
        for i, (a, b, c) in enumerate(triangles):
            x = (a.coords[0]+b.coords[0]+c.coords[0])/3
            y = (a.coords[1]+b.coords[1]+c.coords[1])/3
            if pointInPolygon([x, y], poly):
                yield self.faces[i+1]

    def isConstrained(self, edge):
        a = edge.origin.coords
        b = edge.destination.coords
        if [a, b] in self.polygon or [b, a] in self.polygon:
            return True
        else:
            return False

    def add_point(self):
        """TODO"""
        new_point = None
        face, _ = self.get_face_with_min_angle()
        e1, e2, e3 = (edge.length for edge in face.edges)
        a, b, c = face.vertices
        for angle in get_angles(e1, e2, e3):
            if angle < self.alpha:
                for edge in face.edges:
                    if self.isConstrained(edge):
                        if (edge.length > 2*edge.next.length
                           or edge.length > 2*edge.prev.length
                           or (edge.length > edge.next.length and edge.length > edge.prev.length)):
                            self.splitEdge(edge)
                            return
                else:
                    x = (a.coords[0]+b.coords[0]+c.coords[0])/3
                    y = (a.coords[1]+b.coords[1]+c.coords[1])/3
                    new_point = [x, y]
        if new_point is None:
            return
        puntos = [vertex.coords for vertex in self.vertices]+[new_point]
        D = Dcel(puntos)
        D.polygon = self.polygon
        self.vertices = D.vertices
        self.edges = D.edges
        self.faces = D.faces
        self.enforce_edges()

    @property
    def minimum_angle(self):
        min_angle = 1000
        for face in self.get_interior_triangles():
            e1, e2, e3 = (edge.length for edge in face.edges)
            minimo = min(get_angles(e1, e2, e3))
            if minimo < min_angle:
                min_angle = minimo
        return min_angle

    def splitEdge(self, split):
        new_point = split.mid_point()
        a, b = split.origin.coords, split.destination.coords
        for i, edge in enumerate(self.polygon):
            if (edge[0] == a and edge[1] == b):
                edge[1] = new_point
                self.polygon.insert(i+1, [new_point, b])
                break
            if (edge[0] == b and edge[1] == a):
                edge[1] = new_point
                self.polygon.insert(i+1, [new_point, a])
                break
        puntos = [vertex.coords for vertex in self.vertices]+[new_point]
        self.splitted.append(len(puntos)-1)
        D = Dcel(puntos)
        D.polygon = self.polygon
        self.vertices = D.vertices
        self.edges = D.edges
        self.faces = D.faces
        self.enforce_edges()

    def get_face_with_min_angle(self):
        min_face = None
        min_angle = 1000
        for face in self.get_interior_triangles():
            e1, e2, e3 = (edge.length for edge in face.edges)
            minimo = min(get_angles(e1, e2, e3))
            if minimo < min_angle:
                min_face = face
                min_angle = minimo
        return min_face, min_angle

    def generate_mesh(self, alpha=20, max_iterations=500):
        self.alpha = alpha
        iteration = 0
        while self.minimum_angle < self.alpha and iteration < max_iterations:
            if iteration % 10 == 0:
                self.add_point()
            else:
                self.iterate_forces()
            iteration += 1
