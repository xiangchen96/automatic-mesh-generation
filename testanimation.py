import dcel
import numpy as np
import matplotlib.pyplot as plt
import time

fileName = "demo"

D = dcel.Dcel.deloneFromPolygonFile(fileName)
f = open(fileName,"r")
puntos = []
number_of_points = int(f.readline())
for i in range(number_of_points):
    line = f.readline()
    x,y = line.split(" ")
    puntos.append([float(x),float(y)])
    
alpha = 20

D.alpha = alpha

D.animate_main()
