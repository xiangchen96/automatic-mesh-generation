import mesh_generator

fileName = input("Introduzca el nombre del archivo de puntos:")
D = mesh_generator.Dcel.deloneFromPolygonFile(fileName)

alpha = float(input("Introduzca el alpha (angulo minimo): "))
D.alpha = alpha

print("\nPoligono inicial con angulo minimo %d:\n"%D.get_minimun_angle())
D.plotPolygon()

D.animate_main()
