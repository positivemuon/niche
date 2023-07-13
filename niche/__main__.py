import sys
import argparse
import numpy as np
from pymatgen.core import Structure

from . import Niche

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute interstitial positions")
    parser.add_argument(
        "--spacing",
        metavar="S",
        type=float,
        help="Maximum grid point spacing (Ang.)",
        default=1.0,
    )
    parser.add_argument(
        "--distance",
        metavar="D",
        type=float,
        default=1.0,
        help="Distance between interstitial grid points and lattice points (Ang.)",
    )
    parser.add_argument(
        "--atom",
        metavar="A",
        type=str,
        default="H",
        help="Atom to be inserted in structures",
    )
    parser.add_argument(
        "--covalent-atom",
        metavar="C",
        type=str,
        default="",
        help="Atom symbol to be cosidered for covalence radius calc.",
    )
    parser.add_argument(
        "--supercell",
        metavar="SC",
        type=str,
        default="1 0 0  0 1 0  0 0 1",
        help="Supercell to be built",
    )

    parser.add_argument("input_structure")

    args = parser.parse_args()
    spacing = args.spacing
    distance = args.distance
    atom = args.atom
    cov_atom = args.covalent_atom
    sc_matrix = (
        np.array([float(x) for x in args.supercell.split()]).reshape(-1, 3).squeeze()
    )

    # load structure with pymatgen
    st = Structure.from_file(args.input_structure)

    nc = Niche(st, atom)
    st_w_imps = nc.apply(spacing, distance)
    st_w_imps.to(filename="positions.cif")

    i = 1
    for site in st_w_imps:
        if site.species_string != atom:
            continue
        sc = nc.build_supercell(st, atom, site.frac_coords, sc_matrix)
        sc.perturb(0.0, 0.01)
        sc.to(filename="positions_{}.cif".format(i))
        i += 1
