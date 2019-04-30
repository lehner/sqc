#
# Christoph Lehner 2019
#
import numpy as np
from sqc.state import state
from sqc.operator import operator

def seed(s):
    np.random.seed(s)

def sample(s, n):
    rc={}
    for i in range(n):
        v=0
        s0=s
        for j in range(s.nbits):
            s1,v1=s0.measure(j)
            v=v + 2**j * v1
            s0=s1
        if not v in rc:
            rc[v]=1
        else:
            rc[v]+=1
    return rc
