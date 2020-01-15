# Generador de mallas de regiones poligonales
<img src="https://thumbs.gfycat.com/AppropriateHighlevelBobwhite-size_restricted.gif" width="400" height="300"/>

Dado un polígono simple y un ángulo `α`,
genera una malla de triángulos con ángulo mínimo `>= α`.
Hace uso de la inserción de puntos de Steiner,
partición de aristas restringidas y método de las fuerzas.

Trabajo Fin de Grado
- Autor: Xiang Chen Chen
- Tutor: Manuel Abellanas Oar

## Ejemplos de uso
Básico:
```
python main.py data/demo
```
Ángulo minimo 8:
```
python main.py data/demo -a 8
```
## Fichero de entrada
El polígono estará representado por un conjunto de puntos ordenados. El fichero de entrada debe contener el número total de puntos y sus coordenadas x,y.
```
4
-2.42 3.02
-2.82 3.02
-3.22 2.96
-3.62 2.9
```

## Librerías externas
```
matplotlib
numpy
scipy
```