#
# Christoph Lehner 2019
#
from sqc.operator import _X
import numpy as np

class generic:
    def __init__(self, noisy_gates, readout_errors):
        self.noisy_gates = noisy_gates
        self.readout_errors = readout_errors

    def sample(self, op, st):
        st0=st.clone()

        for o in op.m:
            if not o[2] in self.noisy_gates:
                print("Noise model incomplete: %s" % o[2])
                assert(0)
            st0=self.noisy_gates[o[2]](o[0],o[1],st0,op)
            
        for i,e in enumerate(self.readout_errors):
            assert(e <= 0.5)
            x=np.random.uniform()
            if x < e:
                st0=_X(i,st0,op)

        return st0

