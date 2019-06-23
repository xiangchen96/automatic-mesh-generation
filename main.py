import argparse

from mesh_generator.classes import Dcel

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file")
parser.add_argument("-A", "--animate", help="animate", action="store_true")
parser.add_argument("-a", "--alpha", type=float, default=15, help="minimum output angle")
parser.add_argument("-o", "--output", type=str, default=None, help="output file")
args = parser.parse_args()


D = Dcel.deloneFromPolygonFile(args.input)
with open(args.input, 'r') as f:
    number_of_points = int(f.readline())

print(f'\nPoligono inicial con angulo minimo {D.minimum_angle:.2f}:\n')
D.plotPolygon()

if args.animate:
    D.alpha = args.alpha
    D.animate_main(args.output)
else:
    D.generate_mesh(args.alpha, 1000)
    D.plotPolygon()

    print(f'angulo: {D.minimum_angle:.2f}')
    print(f'vertices nuevos: {len(D.vertices)-number_of_points}')
