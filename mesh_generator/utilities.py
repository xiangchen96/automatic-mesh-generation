import matplotlib.path as mpltPath
from math import acos, degrees


def sarea(a, b, c):
    """Calculate the signed area of 3 points."""
    return 0.5 * ((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1]))


def segment_crossing(s1, s2):
    """Return true if the two segments intersect."""
    a, b, c, d = s1[0], s1[1], s2[0], s2[1]
    a1 = sarea(a, b, c)*sarea(a, b, d) < 0
    a2 = sarea(c, d, a)*sarea(c, d, b) < 0
    return a1 and a2


def pointInPolygon(point, polygon):
    """Return True if the point is interior to the polygon"""
    path = mpltPath.Path(polygon)
    return path.contains_points([point])


def get_angles(a, b, c):
    """Return the angles of a triangle with the provided side lengths."""
    def get_angle(a, b, c): return degrees(acos((b*b+c*c-a*a)/(float(2*b*c))))
    if a + b <= c or b + c <= a or c + a <= b:
        return [0, 0, 0]
    return [get_angle(*_) for _ in ((a, b, c), (b, c, a), (c, a, b))]


def in_circle(a, b, c, point):
    """Return true if d is in the circle a,b,c"""
    sa = sarea(a, b, c)
    if sa == 0:
        return
    return -np.sign(sa * svolume(a, b, c, point))


def svolume(a, b, c, d):
    """Calculate the signed volume"""
    last_line = [a[0]**2+a[1]**2, b[0]**2+b[1]**2,
                 c[0]**2+c[1]**2, d[0]**2+d[1]**2]
    arr = np.array([[1, 1, 1, 1],
                    [a[0], b[0], c[0], d[0]],
                    [a[1], b[1], c[1], d[1]],
                    last_line])
    return np.linalg.det(arr)/6


def project_vector(vector, vector_project):
    """Proyeccion vectorial de vector sobre vector_project"""
    aux = vector_project[0]*vector[0] + vector_project[1]*vector[1]
    aux = aux / (vector_project[0]**2 + vector_project[1]**2)
    return [aux*vector_project[0], aux*vector_project[1]]
