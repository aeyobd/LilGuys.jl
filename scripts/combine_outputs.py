#!/usr/bin/env python

import h5py
from glob import glob
import re
import numpy as np

path = "."
outfile = "combined.hdf5"

filenames = glob(path + "/snapshot*.hdf5")
ids = [re.search(r'\d+', f).group(0) for f in filenames]
ids = np.array(ids, dtype=int)
idx = np.argsort(ids)
ids = ids[idx]
filenames = np.array(filenames)[idx]

with h5py.File(outfile, "w") as f:
    for i, filename in zip(ids, filenames):
        f["snap{0}".format(i)] = h5py.ExternalLink(filename, "/")
