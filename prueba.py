# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 10:54:50 2018

@author: Stradi
"""
import random
import numpy as np
import matplotlib.pyplot as plt
import gtc
from matplotlib import collections  as mc
from scipy.spatial import Delaunay

def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return [pool[i] for i in indices]

polyP = gtc.randomPolyPoints(20,30)
"""PolyP, Delone, CDT"""
lines =  polyP[0]
D = gtc.constrainedDelaunay(polyP[1],polyP[0])
lc = mc.LineCollection(lines, linewidths=2,color='k')
gtc.plotPolyPoints(polyP)
gtc.plotDCEL(gtc.delone(polyP[1]))
gtc.plotDCEL(D,lc)

"""Triangulate polygon"""
#gtc.plotPolyPoints(polyP)
#points, simplices = gtc.triangulatePolyPoints(polyP)
#plt.triplot(points[:,0], points[:,1], simplices)
#plt.plot(points[:,0], points[:,1], 'bo')
#plt.show()