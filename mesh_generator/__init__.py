import matplotlib.pyplot as plt
import numpy as np
import collections

from matplotlib import animation
from scipy.spatial import Delaunay

from .utilities import (
    sarea,
    segment_crossing,
    point_in_polygon,
    get_angles,
    in_circle,
    project_vector,
)


class Vertex:
    """DCEL Vertex"""

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
        return (edge.next.origin for edge in self.edges)

    @property
    def force_vector(self):
        """Vector suma de vertices adyacentes"""
        vector = [0, 0]
        for vertex in self.neighbours:
            x, y = vertex.coords
            my_x, my_y = self.coords
            vector[0] += x - my_x
            vector[1] += y - my_y
        return vector

    def add_force_vector(self, point=None, step=0.05):
        """Mueve el punto en la direccion del vector suma (proyectado)"""
        force = self.force_vector

        if point:
            vector_project = [point[0] - self.coords[0], point[1] - self.coords[1]]
            force = project_vector(force, vector_project)

        self.coords[0] = self.coords[0] + force[0] * step
        self.coords[1] = self.coords[1] + force[1] * step


class Edge:
    """Half Edge"""

    def __init__(self, origin, twin=None, prev=None, next_=None, face=None):
        self.origin = origin
        self.twin = twin
        self.prev = prev
        self.next = next_
        self.face = face

    @property
    def length(self):
        vector = [
            self.origin.coords[0] - self.next.origin.coords[0],
            self.origin.coords[1] - self.next.origin.coords[1],
        ]
        return np.linalg.norm(vector)

    @property
    def destination(self):
        return self.next.origin

    @property
    def mid_point(self):
        origin = self.origin.coords
        destination = self.next.origin.coords
        x = (origin[0] + destination[0]) / 2
        y = (origin[1] + destination[1]) / 2
        return [x, y]

    def flip(self):
        """
        Flip an edge into the oposite diagonal of the quadilateral
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
        a = self.origin.coords
        b = self.twin.prev.origin.coords
        c = self.next.origin.coords
        d = self.next.next.origin.coords
        return -1 not in np.sign(
            [sarea(a, b, c), sarea(a, c, d), sarea(b, c, d), sarea(b, d, a)]
        )

    def is_legal(self, face0):
        if self.face == face0 or self.twin.face == face0:
            return True
        a = self.origin.coords
        b = self.twin.prev.origin.coords
        c = self.next.origin.coords
        d = self.next.next.origin.coords
        if -1 in np.sign(
            [sarea(a, b, c), sarea(a, c, d), sarea(b, c, d), sarea(b, d, a)]
        ):
            return True
        else:
            return in_circle(a, c, d, b) == -1


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
        return (edge.origin for edge in self.edges)


class Dcel:
    def __init__(self, points):
        tesselation = Delaunay(points)
        self.polygon = None
        self.alpha = None
        self.min_x = points[0][0]
        self.max_x = points[0][0]
        self.min_y = points[0][1]
        self.max_y = points[0][1]
        edges = collections.defaultdict()
        self.vertices = []
        self.faces = [Face()]
        self.edges = []
        self.splitted = []
        for point in points:
            self.vertices.append(Vertex(point))
            if point[0] < self.min_x:
                self.min_x = point[0]
            if point[0] > self.max_x:
                self.max_x = point[0]
            if point[1] < self.min_y:
                self.min_y = point[1]
            if point[1] > self.max_y:
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
                edges[(a, b)] = Edge(
                    self.vertices[a], twin=edges[b, a], face=self.faces[0]
                )
                edges[(b, a)].twin = edges[(a, b)]
                hull.append((a, b))
            elif (b, a) not in edges:
                edges[(b, a)] = Edge(
                    self.vertices[a], twin=edges[a, b], face=self.faces[0]
                )
                edges[(a, b)].twin = edges[(b, a)]
                hull.append((b, a))
        for a, b in hull:
            for c, d in hull:
                if c == b:
                    edges[a, b].next = edges[c, d]
                    edges[c, d].prev = edges[a, b]
                    break

    @property
    def minimum_angle(self):
        result = 1000
        for face in self.interior_triangles:
            e1, e2, e3 = (edge.length for edge in face.edges)
            minimum = min(get_angles(e1, e2, e3))
            if minimum < result:
                result = minimum
        return result

    @property
    def interior_triangles(self):
        poly = [i[0] for i in self.polygon]
        triangles = [face.vertices for face in self.faces[1:]]
        for i, (a, b, c) in enumerate(triangles):
            x = (a.coords[0] + b.coords[0] + c.coords[0]) / 3
            y = (a.coords[1] + b.coords[1] + c.coords[1]) / 3
            if point_in_polygon([x, y], poly):
                yield self.faces[i + 1]

    @property
    def face_with_min_angle(self):
        min_face = None
        min_angle = 1000
        for face in self.interior_triangles:
            e1, e2, e3 = (edge.length for edge in face.edges)
            minimo = min(get_angles(e1, e2, e3))
            if minimo < min_angle:
                min_face = face
                min_angle = minimo
        return min_face, min_angle

    @classmethod
    def delone_from_file(cls, filename):
        f = open(filename, "r")
        points = []
        number_of_points = int(f.readline())
        for i in range(number_of_points):
            line = f.readline()
            x, y = line.split(" ")
            points.append([float(x), float(y)])
        d = cls(points)
        d.polygon = [
            [points[i], points[(i + 1) % len(points)]] for i in range(len(points))
        ]
        d.enforce_edges()
        f.close()
        return d

    @classmethod
    def delone_from_polygon(cls, points):
        d = cls(points)
        d.polygon = [
            [points[i], points[(i + 1) % len(points)]] for i in range(len(points))
        ]
        d.enforce_edges()
        return d

    def plot_polygon(self):
        plt.axes().set_aspect("equal")
        for face in self.interior_triangles:
            a, b, c = face.vertices
            plt.triplot(
                [a.coords[0], b.coords[0], c.coords[0]],
                [a.coords[1], b.coords[1], c.coords[1]],
                "bo-",
            )
        plt.show()

    def plot(self):
        plt.axes().set_aspect("equal")
        for face in self.faces[1:]:
            a, b, c = face.vertices
            plt.triplot(
                [a.coords[0], b.coords[0], c.coords[0]],
                [a.coords[1], b.coords[1], c.coords[1]],
                "bo-",
            )
        plt.show()

    def contains_edge(self, searched_edge):
        for edge in self.edges:
            o = edge.origin.coords
            d = edge.next.origin.coords
            if [o, d] == searched_edge or [d, o] == searched_edge:
                return True
        return False

    def enforce_edges(self):
        for (Vi, Vj) in self.polygon:
            if self.contains_edge([Vi, Vj]):
                continue
            new_edges = []
            crossing_edges = []
            for edge in self.edges:
                if edge.twin in crossing_edges:
                    continue
                if segment_crossing(
                    [Vi, Vj], [edge.origin.coords, edge.next.origin.coords]
                ):
                    crossing_edges.append(edge)
            while len(crossing_edges) > 0:
                e = crossing_edges.pop()
                if not e.is_flippable(self.faces[0]):
                    crossing_edges.insert(0, e)
                else:
                    e.flip()
                    if segment_crossing(
                        [Vi, Vj], [e.origin.coords, e.next.origin.coords]
                    ):
                        crossing_edges.insert(0, e)
                    else:
                        new_edges.append(e)
            swap = True
            while swap:
                swap = False
                for e in new_edges:
                    if e.origin.coords in [Vi, Vj] and e.next.origin.coords in [Vi, Vj]:
                        continue
                    if not e.is_legal:
                        e.flip()
                        swap = True

    def iterate_forces(self):
        polygon_vertices = [edge[0] for edge in self.polygon]
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
        d = Dcel([vertex.coords for vertex in self.vertices])
        self.vertices = d.vertices
        self.edges = d.edges
        self.faces = d.faces
        self.enforce_edges()

    def animate_main(self, output=None):
        fig = plt.figure()
        plt.axes(
            xlim=(self.min_x - 1, self.max_x + 1), ylim=(self.min_y - 1, self.max_y + 1)
        )
        angle_text = plt.text(
            (self.max_x + self.min_x) / 2, self.max_y, "", fontsize=10
        )

        lines = [plt.plot([], [], "bo-")[0] for _ in range(len(self.edges))]

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
            angle_text.set_text(f"min_angle: {angle:.2f} iter: {frame}")
            edges = []
            for face in self.interior_triangles:
                for edge in face.edges:
                    if [edge.next.origin.coords, edge.origin.coords] not in edges:
                        edges.append([edge.origin.coords, edge.next.origin.coords])
            for i, edge in enumerate(edges):
                if len(lines) > i:
                    lines[i].set_data(
                        [edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]]
                    )
                else:
                    lines.append(
                        plt.plot(
                            [edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], "bo-"
                        )[0]
                    )
            return lines + [angle_text]

        ani = animation.FuncAnimation(
            fig, animate, init_func=init, interval=10, blit=True
        )
        if output:
            if not output.endswith(".mp4"):
                output += ".mp4"
            ani.save(output, fps=30, extra_args=["-vcodec", "libx264"])
        else:
            plt.show()

    def is_constrained(self, edge):
        o = edge.origin.coords
        d = edge.destination.coords
        for (a, b) in self.polygon:
            if (o, d) == (a, b) or (o, d) == (b, a):
                return True
        return False

    def add_point(self):
        """TODO"""
        new_point = None
        face, _ = self.face_with_min_angle
        e1, e2, e3 = (edge.length for edge in face.edges)
        a, b, c = face.vertices
        for angle in get_angles(e1, e2, e3):
            if angle < self.alpha:
                for edge in face.edges:
                    if self.is_constrained(edge):
                        if (
                            edge.length > 2 * edge.next.length
                            or edge.length > 2 * edge.prev.length
                            or (
                                edge.length > edge.next.length
                                and edge.length > edge.prev.length
                            )
                        ):
                            self.split_edge(edge)
                            return
                else:
                    x = (a.coords[0] + b.coords[0] + c.coords[0]) / 3
                    y = (a.coords[1] + b.coords[1] + c.coords[1]) / 3
                    new_point = [x, y]
        if new_point is None:
            return
        puntos = [vertex.coords for vertex in self.vertices] + [new_point]
        d = Dcel(puntos)
        d.polygon = self.polygon
        self.vertices = d.vertices
        self.edges = d.edges
        self.faces = d.faces
        self.enforce_edges()

    def split_edge(self, split):
        new_point = split.mid_point
        a, b = split.origin.coords, split.destination.coords
        for i, edge in enumerate(self.polygon):
            if edge[0] == a and edge[1] == b:
                edge[1] = new_point
                self.polygon.insert(i + 1, [new_point, b])
                break
            if edge[0] == b and edge[1] == a:
                edge[1] = new_point
                self.polygon.insert(i + 1, [new_point, a])
                break
        puntos = [vertex.coords for vertex in self.vertices] + [new_point]
        self.splitted.append(len(puntos) - 1)
        d = Dcel(puntos)
        d.polygon = self.polygon
        self.vertices = d.vertices
        self.edges = d.edges
        self.faces = d.faces
        self.enforce_edges()

    def generate_mesh(self, alpha=20, max_iterations=500):
        self.alpha = alpha
        iteration = 0
        while self.minimum_angle < self.alpha and iteration < max_iterations:
            if iteration % 10 == 0:
                self.add_point()
            else:
                self.iterate_forces()
            iteration += 1
