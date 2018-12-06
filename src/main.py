import mesh_generator

fileName = input('Introduzca el nombre del archivo de puntos:')

D = mesh_generator.Dcel.deloneFromPolygonFile(fileName)
with open(fileName, 'r') as f:
    number_of_points = int(f.readline())

alpha = float(input('Introduzca el alpha (angulo minimo):'))

print(f'\nPoligono inicial con angulo minimo {D.minimum_angle:.2f}:\n')
D.plotPolygon()

D.generate_mesh(alpha, 1000)
D.plotPolygon()

print(f'angulo: {D.minimum_angle:.2f}')
print(f'vertices nuevos: {len(D.vertices)-number_of_points}')
