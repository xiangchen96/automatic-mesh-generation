import mesh_generator

fileName = input("Introduzca el nombre del archivo de puntos:")

D = mesh_generator.Dcel.deloneFromPolygonFile(fileName)
with open(fileName,'r') as f:
    number_of_points = int(f.readline())

alpha = float(input("Introduzca el alpha (angulo minimo):"))

print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()

D.generate_mesh(alpha,1000)
D.plotPolygon()

print("angulo:",D.get_minimun_angle())
print("vertices nuevos:", len(D.vertices)-number_of_points)