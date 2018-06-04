import dcel
import numpy as np
import triangle
import triangle.plot as plot
import matplotlib.pyplot as plt

fileName = input("Introduzca el nombre del archivo de puntos:")
#fileName = "puntos2"

D = dcel.Dcel.deloneFromPolygonFile(fileName)
f = open(fileName,"r")
puntos = []
number_of_points = int(f.readline())
for i in range(number_of_points):
    line = f.readline()
    x,y = line.split(" ")
    puntos.append([float(x),float(y)])
    
alpha = int(input("Introduzca el alpha (angulo minimo):"))
#alpha = 20

D.alpha = alpha

print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plot()

points_and_segments = {'vertices':np.array(puntos), 
                       'segments':np.array([[i,(i+1)%len(puntos)] for i in range(len(puntos))])} 
triangulation = triangle.triangulate(points_and_segments, 'q'+str(alpha))

ax = plt.axes()
ax.set_aspect('equal')
plot.plot(ax, **triangulation)
plt.show()
print("vertices:",len(triangulation['vertices']))

#D.animate_main()

D.generate_mesh()
D.plot()

print("angulo ",D.get_minimun_angle())
print("vertices ", len(D.vertices))