#
# Christoph Lehner 2019
#
import numpy as np
from sqc.state import state
from sqc.operator import operator
import sqc.noise

def seed(s):
    np.random.seed(s)

def sample(s, n, mask = None, save_states = False):
    rc={}
    rs={}
    if mask is None:
        mask=range(s.nbits)
    for i in range(n):
        v=0
        s0=s
        for l,j in enumerate(mask):
            s1,v1=s0.measure(j)
            v=v + 2**l * v1
            s0=s1
        if not v in rc:
            rc[v]=1
            if save_states:
                rs[v] = s0
        else:
            rc[v]+=1
    if save_states:
        return rc,rs
    return rc
