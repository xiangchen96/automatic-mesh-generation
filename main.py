import sys
import argparse

from mesh_generator import Dcel

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input", type=str, help="input file")
    parser.add_argument(
        "-a", "--alpha", type=float, default=15, help="minimum output angle"
    )

    # Help if no args
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    D = Dcel.delone_from_file(args.input)
    D.alpha = args.alpha
    D.animate_main()
