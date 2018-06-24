import dcel

fileName = input("Introduzca el nombre del archivo de puntos:")

D = dcel.Dcel.deloneFromPolygonFile(fileName)
    
alpha = int(input("Introduzca el alpha (angulo minimo):"))
D.alpha = alpha

print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()

D.generate_mesh()
D.plotPolygon()

print("angulo:",D.get_minimun_angle())
print("vertices nuevos:", len(D.vertices)-number_of_points)