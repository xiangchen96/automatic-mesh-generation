import dcel
import numpy as np
import matplotlib.pyplot as plt
from time import time

fileName = "demo"
alpha = 10
tiempos = []
for i in range(1):
	t0 = time()
	D = dcel.Dcel.deloneFromPolygonFile(fileName)
	D.alpha = alpha
	D.generate_mesh()
	tiempos.append(time()-t0)
	
print(np.mean(tiempos))

