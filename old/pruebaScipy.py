# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 19:06:32 2018

@author: Stradi
"""
import numpy as np
import matplotlib.pyplot as plt
import gtc


points = np.array([[np.random.rand(),np.random.rand()] for _ in range(4)])
from scipy.spatial import Delaunay
tri = Delaunay(points,False,True)
#tri.add_points(np.array([[1,0],[1,1]]))

plt.triplot(tri.points[:,0], tri.points[:,1], tri.simplices.copy())
plt.plot(tri.points[:,0], tri.points[:,1], 'ro')
for i,(a,b) in enumerate(tri.points):
    plt.text(tri.points[i,0],tri.points[i,1]+0.01,i)
plt.show()


print(tri.simplices)
print(tri.neighbors)
D = gtc.dcelFromDelaunay(tri)
#print(D[1])