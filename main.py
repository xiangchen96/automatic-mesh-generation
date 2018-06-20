import dcel
import numpy as np
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

D.generate_mesh()
D.plotPolygon()

print("angulo:",D.get_minimun_angle())
print("vertices nuevos:", len(D.vertices)-number_of_points)