#
# Christoph Lehner 2019
#
import numpy as np

class state:
    def __init__(self, nbits, v = None, basis = None):
        self.nbits = nbits
        self.N = 2**nbits
        if v == None:
            self.v = np.array(([ 0.0 ] * (self.N - 1)) + [ 1.0 ], dtype=np.cdouble)
        else:
            self.v = np.array(v)
        if basis == None:
            self.basis = [ self.default_basis_string(i) for i in range(self.N) ]
        else:
            self.basis = basis
        assert(self.nbits < 31) # don't run out of bits

    def clone(self):
        return state(self.nbits, v = self.v.copy(), basis = self.basis)

    def default_basis_string(self, i):
        res=""
        for j in range(self.nbits):
            b=2**j
            if i & b == b:
                res+="|1>"
            else:
                res+="|0>"
        return res

    def __str__(self):
        coef=dict([ (i,"%s" % (str(c))) for (i,c) in enumerate(self.v) if abs(c) != 0.0 ])
        if len(coef) == 0:
            return "0"

        slen=max([ len(c) for c in coef.values() ])
        lines=[ "%s%s * %s" % (coef[i]," "*(slen-len(coef[i])),self.basis[i]) for i in coef ]
        return "   " + "\n + ".join(lines)
