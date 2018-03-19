# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 19:06:32 2018

@author: Stradi
"""
import numpy as np
import matplotlib.pyplot as plt


points = np.array([[np.random.rand(),np.random.rand()] for _ in range(50)])
from scipy.spatial import Delaunay
tri = Delaunay(points,False,True)
tri.add_points(np.array([[1,0],[1,1]]))
points = np.append(points,[[1,0],[1,1]],axis=0)
plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.show()