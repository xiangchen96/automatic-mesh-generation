import triangle
import matplotlib.pyplot as plt
import numpy as np
from dcel import Dcel
import utilidades as utils

def filePoints(fileName):
    f = open(fileName,"r")
    points = []
    number_of_points = int(f.readline())
    for i in range(number_of_points):
        line = f.readline()
        x,y = line.split(" ")
        points.append([float(x),float(y)])
    f.close()
    return {'vertices':np.array(points),
            'segments':np.array([[i,(i+1)%len(points)] for i in range(len(points))])} 


fileName = "puntos1"
points_and_segments = filePoints(fileName)
angle = 30

triangulation = triangle.triangulate(points_and_segments, 'q'+str(angle))
plt.axes().set_aspect('equal')
plt.triplot(triangulation['vertices'][:,0], triangulation['vertices'][:,1], triangulation['triangles'],'ro-')
plt.show()
print("vertices ", len(triangulation['vertices']))

D = Dcel.deloneFromPolygonFile(fileName)
D.alpha = angle
D.generate_mesh()
D.plot()

print("angulo ",D.get_minimun_angle())
print("vertices ", len(D.vertices))

#puntos = list(map(list,np.random.rand(10,2)))
#puntos = utils.polygonization(puntos)
#
#D = Dcel.deloneFromPolygon(puntos)
#D.plot()