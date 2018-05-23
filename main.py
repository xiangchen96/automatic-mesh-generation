import dcel
import numpy as np
import utilidades as utils
import triangle
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
fileName = input("Introduzca el nombre del archivo de puntos, se generara un poligono aleatorio si se deja vacio:")
if fileName == '':
    puntos = list(map(list,np.random.rand(30,2)))
    puntos = utils.polygonization(puntos)
    D = dcel.Dcel.deloneFromPolygon(puntos)
else:
    D = dcel.Dcel.deloneFromPolygonFile(fileName)
alpha = int(input("Introduzca el alpha (angulo minimo):"))
D.alpha = alpha
print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()

points_and_segments = {'vertices':np.array(puntos), 
                       'segments':np.array([[i,(i+1)%len(puntos)] for i in range(len(puntos))])} 
triangulation = triangle.triangulate(points_and_segments, 'q'+str(angle))
plt.axes().set_aspect('equal')
lines =  [[puntos[i],puntos[(i+1)%len(puntos)]] for i in range(len(puntos))]
lc = mc.LineCollection(lines, linewidths=2)
plt.triplot(triangulation['vertices'][:,0], triangulation['vertices'][:,1], triangulation['triangles'],'ro-')
plt.axes().add_collection(lc)
plt.show()

D.generate_mesh()
print("\nMalla generada con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()
print("angulo ",D.get_minimun_angle())
print("vertices ", len(D.vertices))