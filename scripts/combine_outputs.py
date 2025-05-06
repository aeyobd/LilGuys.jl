#!/usr/bin/env python

import h5py
from glob import glob
import re
import numpy as np

from sys import argv
import argparse
import os.path


def main():
    args = parse_args()
    ids, filenames = get_filenames(args.input)
    write_combined(args.output, ids, filenames, args.centres)
    

def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine multiple Gadget outputs into a single file."
        )
    parser.add_argument("input", help="Path to input files")
    parser.add_argument("-c",  "--centres", 
        help="Path to file containing the centre of the simulation, optional", 
        default=None)
    parser.add_argument("-o", "--output", 
        help="Output file name", default="combined.hdf5")

    args = parser.parse_args()
    return args


def get_filenames(path):
    filenames = glob(path + "/snapshot*.hdf5")
    ids = [re.search(r'\d+', os.path.basename(f)).group(0) for f in filenames]
    ids = np.array(ids, dtype=int)
    idx = np.argsort(ids)
    ids = ids[idx]
    filenames = np.array(filenames)[idx]
    if len(filenames) < 1:
        raise Exception("no output files found ")

    return ids, filenames

def write_combined(outfile, ids, filenames, centre_path):
    if os.path.exists(outfile):
        os.remove(outfile)
    with h5py.File(outfile, "w") as f:
        for i, filename in zip(ids, filenames):
            f["snap{0}".format(i)] = h5py.ExternalLink(filename, "/")

        if centre_path:
            f["x_cen"] = h5py.ExternalLink(centre_path, "/positions")
            f["v_cen"] = h5py.ExternalLink(centre_path, "/velocities")



if __name__ == "__main__":
    main()
