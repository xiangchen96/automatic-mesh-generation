# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 10:54:50 2018

@author: Stradi
"""
import numpy as np
import matplotlib.pyplot as plt
import gtc
from matplotlib import collections  as mc
from scipy.spatial import Delaunay


#PP = gtc.randomPolyPoints(10,40)
#gtc.plotPolyPoints(PP)
#
#gtc.constrainedDelaunay(PP[1],[])

#points = [[0,0],[0,1],[1,0],[1,1],[0.5,0.4],[0.5,2]]
#gtc.plotDCEL(gtc.triangulation(points))
cont = 0
while True:
    points = [[np.random.rand(),np.random.rand()] for i in range(50)]
    gtc.plotDCEL(gtc.delone(points))
#    print("mio")
#    tri = Delaunay(points)
#    print("suyo")
    cont = cont +1
    print(cont)

#gtc.constrainedDelaunay(points,[])
