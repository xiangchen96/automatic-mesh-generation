import mesh_generator

fileName = input('Introduzca el nombre del archivo de puntos:')
D = mesh_generator.Dcel.deloneFromPolygonFile(fileName)

alpha = float(input('Introduzca el alpha (angulo minimo): '))
D.alpha = alpha

print(f'\nPoligono inicial con angulo minimo {D.minimum_angle:.2f}:\n')
D.plotPolygon()

D.animate_main()
