import agama
import numpy as np



def sample_potential(pot, time, N, r_max):

    def f(pos):
        np.maximum(pot.density(pos, t=time), 0)

    positions, _, _, _ = agama.sampleNdim(f, N, -r_max * np.ones(3), r_max * np.ones(3)) 

    return positions
