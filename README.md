# Generador de mallas de regiones poligonales
<img src="https://thumbs.gfycat.com/AppropriateHighlevelBobwhite-size_restricted.gif" width="400" height="300"/>

Trabajo: https://oa.upm.es/51503/

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

## Referencias
* Chew, L. Paul (1987). Constrained Delaunay Triangulations. Proceedings of the
Third Annual Symposium on Computational Geometry.
* Sloan, S.W. (1993) A fast algorithm for generating constrained Delanay triangu-
lations. Computers & Structures Vol.47, No.3.
* Triangle. A Two-Dimensional Quality Mesh Generator and Delaunay Triangu-
lator. Online. Available: https://www.cs.cmu.edu/~quake/triangle.html
* Ruppert, J. (1995). A Delaunay Refinement Algorithm for Quality 2-Dimensional
Mesh Generation. Journal of Algorithms.
