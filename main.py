import dcel
import numpy as np
import utilidades as utils


fileName = input("Introduzca el nombre del archivo de puntos, se generara un poligono aleatorio si se deja vacio:")
if fileName == '':
    puntos = list(map(list,np.random.rand(10,2)))
    puntos = utils.polygonization(puntos)
    D = dcel.Dcel.deloneFromPolygon(puntos)
else:
    D = dcel.Dcel.deloneFromPolygonFile(fileName)
alpha = int(input("Introduzca el alpha (angulo minimo):"))
D.alpha = alpha
print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()
D.generate_mesh()
print("\nMalla generada con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()
print("angulo ",D.get_minimun_angle())
print("vertices ", len(D.vertices))