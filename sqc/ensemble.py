#
# Christoph Lehner 2019
#
import numpy as np
from sqc.state import fmtstr
from sqc.state import state

class ensemble:
    def __init__(self, v):
        if type(v) == state:
            self.states = [ (1,v) ]
        else:
            self.states = v
            assert( abs(sum( e[0] for e in self.states ) - 1) < 1e-13 )
        self.nbits = self.states[0][1].nbits
        for e in self.states:
            assert(e[1].nbits == self.nbits)

    def __str__(self):
        return str(self.states)

    def __repr__(self):
        return self.__str__()

    def measure(self, b):
        x=np.random.uniform()
        p0=0
        for e in self.states:
            p1=p0+e[0]
            if p0 <= x and x < p1:
                return e[1].measure(b)
            p0=p1
        assert(0)
