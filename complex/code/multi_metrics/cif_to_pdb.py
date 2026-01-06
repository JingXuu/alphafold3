#!/usr/bin/env python
"""
Convert AF3-generated .cif / .mmcif structure to .pdb format.

Usage:
    python cif_to_pdb.py input.cif [-o output.pdb]

也可以在其他脚本中：
    from cif_to_pdb import cif_to_pdb
    pdb_path = cif_to_pdb("xxx_model.cif")
"""

import os
import argparse
from biotite.structure.io import load_structure
from biotite.structure.io.pdb import PDBFile


def cif_to_pdb(cif_path: str, pdb_path: str | None = None) -> str:
    """
    Convert a .cif / .mmcif file to .pdb using biotite.

    - If pdb_path is None, output is <cif_path_without_ext>.pdb
    - If the pdb_path already exists, it will be reused (not overwritten).

    Returns:
        The path to the PDB file.
    """
    cif_path = os.path.abspath(cif_path)
    if pdb_path is None:
        base, _ = os.path.splitext(cif_path)
        pdb_path = base + ".pdb"
    else:
        pdb_path = os.path.abspath(pdb_path)

    if os.path.isfile(pdb_path):
        # 已经存在就直接复用，避免重复转换
        print(f"[cif_to_pdb] Reuse existing PDB: {pdb_path}")
        return pdb_path

    if not os.path.isfile(cif_path):
        raise FileNotFoundError(f"Input CIF not found: {cif_path}")

    print(f"[cif_to_pdb] Converting {cif_path} -> {pdb_path}")
    atom_array = load_structure(cif_path)
    pdb_file = PDBFile()
    pdb_file.set_structure(atom_array)
    pdb_file.write(pdb_path)

    return pdb_path


def main():
    parser = argparse.ArgumentParser(description="Convert .cif / .mmcif to .pdb")
    parser.add_argument("input", help="input .cif / .mmcif file")
    parser.add_argument("-o", "--output", help="output .pdb file (optional)")
    args = parser.parse_args()

    pdb_path = cif_to_pdb(args.input, args.output)
    print(f"[cif_to_pdb] Done: {pdb_path}")


if __name__ == "__main__":
    main()
