#
# Christoph Lehner 2019
#
import numpy as np

def fmtstr(c):
    eps=1e-13
    if type(c) == type(0.0):
        r=c
        i=0
    else:
        r=c.real
        i=c.imag
    if abs(r) > eps:
        if i > eps:
            return "(%g+%gj)" % (r,i)
        elif i < -eps:
            return "(%g%gj)" % (r,i)
        else:
            if r > 0.0:
                return "%g" % r
            else:
                return "(%g)" % r
    else:
        if abs(i) > eps:
            return "%gj" % i
        else:
            return "0"


class state:
    def __init__(self, nbits, v = None, basis = None):
        self.nbits = nbits
        self.N = 2**nbits
        if v is None:
            self.v = np.array([ 1.0 ] + ([ 0.0 ] * (self.N - 1)), dtype=np.cdouble)
        else:
            self.v = np.array(v)
        if basis is None:
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
                res="1"+res
            else:
                res="0"+res
        return "|" + res + ">"

    def _chop(self):
        tol = 1e-14
        self.v[abs(self.v) < tol] = 0.0
        self.v = np.real_if_close(self.v,tol=100)
        self.v /= np.sqrt(np.linalg.norm(self.v))
        return self

    def measure(self, b):
        A=sum([ self.v[i]*self.v[i].conj() for i in range(self.N) if i & 2**b != 0 ]).real
        B=sum([ self.v[i]*self.v[i].conj() for i in range(self.N) if i & 2**b == 0 ]).real
        assert(abs(1.0 - A - B) < 1e-11)
        x=np.random.uniform()
        if x < A:
            # project to |1>
            assert(A > 0.0)
            l = 1.0 / np.sqrt(A)
            v = np.array([ l*self.v[i] if i & 2**b != 0 else 0 for i in range(self.N) ])
            r = 1
        else:
            # project to |0>
            assert(B > 0.0)
            l = 1.0 / np.sqrt(B)
            v = np.array([ l*self.v[i] if i & 2**b == 0 else 0 for i in range(self.N) ])
            r = 0

        return (state(self.nbits, v=v, basis=self.basis),r)

    def __str__(self):
        coef=dict([ (i,"%s" % (fmtstr(c))) for (i,c) in enumerate(self.v) if abs(c) > 1e-13 ])
        if len(coef) == 0:
            return "0"

        slen=max([ len(c) for c in coef.values() ])
        lines=[ "%s%s * %s" % (coef[i]," "*(slen-len(coef[i])),self.basis[i]) for i in coef ]
        return "   " + "\n + ".join(lines)

    def __repr__(self):
        return self.__str__()
