# Generador de mallas de regiones poligonales
Trabajo Fin de Grado
- Autor: Xiang Chen Chen
- Tutor: Manuel Abellanas Oar

<img src="https://thumbs.gfycat.com/AppropriateHighlevelBobwhite-size_restricted.gif" width="400" height="300"/>
## Descripción
Dado un polígono simple y un ángulo *α*, genera una malla de triángulos con ángulo mínimo >= *α*. Hace uso de la inserción de puntos de Steiner, partición de aristas restringidas y método de las fuerzas.

## Ejemplo de uso
```python
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
```
## Fichero de entrada
El fichero de entrada debe contener el número de puntos y sus coordenadas x,y.
```
4
-2.42 3.02
-2.82 3.02
-3.22 2.96
-3.62 2.9
```
