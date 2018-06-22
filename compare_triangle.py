import dcel
import numpy as np
import triangle
import triangle.plot as plot
import matplotlib.pyplot as plt
import time

fileName = input("Introduzca el nombre del archivo de puntos:")

D = dcel.Dcel.deloneFromPolygonFile(fileName)
f = open(fileName,"r")
puntos = []
number_of_points = int(f.readline())
for i in range(number_of_points):
    line = f.readline()
    x,y = line.split(" ")
    puntos.append([float(x),float(y)])
    
alpha = int(input("Introduzca el alpha (angulo minimo):"))

D.alpha = alpha

print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()

points_and_segments = {'vertices':np.array(puntos), 
                      'segments':np.array([[i,(i+1)%len(puntos)] for i in range(len(puntos))])} 

triangulation = triangle.triangulate(points_and_segments, 'pq'+str(alpha))

ax = plt.axes()
plot.plot(ax, **triangulation)
plt.show()
print("vertices nuevos:",len(triangulation['vertices'])-number_of_points)

D.generate_mesh()
D.plotPolygon()

print("angulo ",D.get_minimun_angle())
print("vertices nuevos", len(D.vertices)-number_of_points)