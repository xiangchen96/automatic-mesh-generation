import dcel
import numpy as np
import triangle
import triangle.plot as plot
import matplotlib.pyplot as plt
import time

#fileName = input("Introduzca el nombre del archivo de puntos:")
fileName = "puntos2"

D = dcel.Dcel.deloneFromPolygonFile(fileName)
f = open(fileName,"r")
puntos = []
number_of_points = int(f.readline())
for i in range(number_of_points):
    line = f.readline()
    x,y = line.split(" ")
    puntos.append([float(x),float(y)])
    
#alpha = int(input("Introduzca el alpha (angulo minimo):"))
alpha = 20

D.alpha = alpha

print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()

#points_and_segments = {'vertices':np.array(puntos), 
#                       'segments':np.array([[i,(i+1)%len(puntos)] for i in range(len(puntos))])} 
#
#t0 = time.time()
#triangulation = triangle.triangulate(points_and_segments, 'pq'+str(alpha))
#print("Tiempo triangle: ",time.time()-t0)
#
#ax = plt.axes()
#plot.plot(ax, **triangulation)
#plt.show()
#print("vertices añadidos:",len(triangulation['vertices'])-number_of_points)

D.animate_main()
#t0 = time.time()
# D.generate_mesh()
#print("Tiempo fuerzas: ",time.time()-t0)
# D.plotPolygon()

print("angulo ",D.get_minimun_angle())
print("vertices añadidos", len(D.vertices)-number_of_points)