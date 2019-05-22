#
# Christoph Lehner 2019
#
import numpy as np
import copy
import time
from sqc.state import state

def _notbit(l,i):
    return l ^ (2**i)

def getbit(l,i):
    return 1 if (l & 2**i != 0) else 0

def setbit(l,i,v):
    m = 2**i
    l &= ~m
    if v:
        l |= m
    return l

_cache={}
def _cache_init(nbits):
    p=[1.0,-1.0]
    N=2**nbits
    if not nbits in _cache:
        c={}
        c["notmask"]=[ np.array([ _notbit(l,i) for l in range(N) ]) for i in range(nbits) ]
        c["cnotmask"]=[ [ np.array([ _notbit(l,j) if getbit(l,i) else l for l in range(N) ]) for j in range(nbits) ] for i in range(nbits) ]
        c["onemask"]=[ np.array([ 1.0 if l & 2**i != 0 else 0.0 for l in range(N) ]) for i in range(nbits) ]
        c["zeromask"]=[ np.array([ 1.0 if l & 2**i == 0 else 0.0 for l in range(N) ]) for i in range(nbits) ]
        c["sbmask0"]=[ np.array([ setbit(l,i,0) for l in range(N) ]) for i in range(nbits) ]
        c["sbmask1"]=[ np.array([ setbit(l,i,1) for l in range(N) ]) for i in range(nbits) ]
        c["hmask"]=[ np.array([ p[getbit(l,i)] for l in range(N) ]) for i in range(nbits) ]
        _cache[nbits] = c

def _H(i,s,o):
    n=s.nbits
    c=1.0/2.0**0.5
    return state(n,c*(s.v[o.sbmask0[i]] + o.hmask[i]*s.v[o.sbmask1[i]]),basis=s.basis)

def _X(i,s,o):
    n=s.nbits
    return state(n,s.v[o.notmask[i]],basis=s.basis)

def _Rz(i,phi,s,o):
    n=s.nbits
    return state(n,(o.onemask[i]*np.exp(1j*phi)
                    + o.zeromask[i])*s.v,basis=s.basis)

def _CNOT(i,j,s,o): # i is control, j is target
    n=s.nbits
    return state(n,s.v[o.cnotmask[i][j]],basis=s.basis)

class operator:
    def __init__(self, nbits, m = None):
        global _cache
        self.nbits = nbits
        self.N = 2**nbits
        if not m is None:
            self.m = m
        else:
            self.m = []

        _cache_init(self.nbits)
        self.notmask = _cache[self.nbits]["notmask"]
        self.cnotmask = _cache[self.nbits]["cnotmask"]
        self.onemask = _cache[self.nbits]["onemask"]
        self.zeromask = _cache[self.nbits]["zeromask"]
        self.sbmask0 = _cache[self.nbits]["sbmask0"]
        self.sbmask1 = _cache[self.nbits]["sbmask1"]
        self.hmask = _cache[self.nbits]["hmask"]
        self.verbose = False
        self.cache={}

    def clone(self):
        return operator(self.nbits, copy.deepcopy(self.m))

    def __mul__(self, s):
        assert(s.nbits == self.nbits)
        timing={}
        count={}
        for i,ops in enumerate(self.m):
            t0=time.time()
            tag=ops[2]
            s=ops[0](*ops[1],s,self)
            if i % 100 == 0:
                s._chop()
            t1=time.time()
            if not tag in timing:
                timing[tag]=0.0
            if not tag in count:
                count[tag]=0
            timing[tag]+=t1-t0
            count[tag]+=1
        if self.verbose:
            print("-----------------------------------------------------------")
            print(" Ran circuit of length=%d" % len(self.m))
            print("-----------------------------------------------------------")
            for tag in count:
                print(" - %d %s with total time of %gs" % (count[tag],tag,timing[tag]))
            print("-----------------------------------------------------------")
        return s

    def __add__(self, o):
        assert(self.nbits == o.nbits)
        return operator(self.nbits, self.m + o.m )

    def op(self,_f,_a,_s):
        return operator(self.nbits, self.m + [ (_f,_a,_s) ] )

    def H(self, i):
        return self.op(_H,[i],"H")

    def NOT(self, i):
        return self.op(_X,[i],"X")

    def X(self, i):
        return self.NOT(i)

    # i is control bit and j is target bit
    def CNOT(self, i, j):
        assert(i!=j)
        return self.op(_CNOT,[i,j],"CNOT")

    def Rz(self, i, phi):
        return self.op(_Rz,[i,phi],"Rz")
