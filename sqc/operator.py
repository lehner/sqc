#
# Christoph Lehner 2019
#
import numpy as np
from sqc.state import state

stock_H=1.0 / np.sqrt(2.0) * np.array([ [1,1], [1,-1] ])
stock_Z=np.array([ [ 1,0], [0,-1] ])
stock_I=np.array([ [ 1,0], [0,1] ])
stock_X=np.array([ [ 0,1], [1,0] ])
stock_Y=np.array([ [ 0,-1j], [1j,0] ])

# control bit is 1, target bit is 0
stock_CNOT=np.array( [ [1,0,0,0],
                       [0,1,0,0],
                       [0,0,0,1],
                       [0,0,1,0] ] )


class operator:
    def __init__(self, nbits, m = None):
        self.nbits = nbits
        self.N = 2**nbits
        if not m is None:
            self.m = m
        else:
            self.m = np.identity(self.N)

    def clone(self):
        return operator(self.nbits, self.m.copy())

    def __mul__(self, s):
        assert(s.nbits == self.nbits)
        return state(self.nbits,v=np.dot(self.m,s.v),basis=s.basis)

    def nonunitarity(self):
        return np.linalg.norm(np.dot(self.m,np.transpose(np.conjugate(self.m))) - np.identity(self.N))

    def H(self, i):
        return self.gate1(i, stock_H)

    def Z(self, i):
        return self.gate1(i, stock_Z)

    def Y(self, i):
        return self.gate1(i, stock_Y)

    def X(self, i):
        return self.gate1(i, stock_X)

    def NOT(self, i):
        return self.gate1(i, stock_X)

    def I(self, i):
        return self.gate1(i, stock_I)

    # i is control bit and j is target bit
    def CNOT(self, i, j):
        return self.gate2(i, j, stock_CNOT)

    def CX(self, i, j):
        return self.CNOT(i,j)

    def _u3v(self, theta, phi, lam):
        return np.array( [ [ np.cos(theta/2.0), -np.exp(1j*lam)*np.sin(theta/2.0) ],
                           [ np.exp(1j*phi)*np.sin(theta/2.0), np.exp(1j*(lam+phi))*np.cos(theta/2.0) ] ] )

    def u3(self, i, theta, phi, lam):
        u3v=self._u3v(theta,phi,lam)
        return self.gate1(i, u3v)

    def u2(self, i, phi, lam):
        return self.u3(i, np.pi/2.0,phi,lam)

    def u1(self, i, lam):
        return self.u3(i, 0.0, 0.0, lam)

    def S(self, i):
        return self.u1(i, np.pi/2.0)

    def Sdg(self, i):
        return self.u1(i, -np.pi/2.0)

    def T(self, i):
        return self.u1(i, np.pi/4.0)

    def Tdg(self, i):
        return self.u1(i, -np.pi/4.0)

    def Rx(self, i, theta):
        return self.u3(i,theta, -np.pi/2.0, np.pi/2.0) 

    def Ry(self, i, theta):
        return self.u3(i,theta, 0.0, 0.0)

    def Rz(self, i, phi):
        return self.u1(i, phi)

    def _idx2bit(self, i):
        return [ (i & 2**l) != 0 for l in range(self.nbits) ]

    def _bit2idx(self, b):
        return sum([ 2**i if b[i] else 0 for i in range(self.nbits) ])

    def _swapbit(self, l, i, j):
        c=[ a for a in l ]
        t=c[i]
        c[i]=c[j]
        c[j]=t
        return c
        
    def _swapmat(self, i, j):
        order=[ self._bit2idx(self._swapbit(self._idx2bit(b),i,j)) for b in range(self.N) ]
        return np.array( [ [ 1.0 if order[a] == b else 0.0 for a in range(self.N) ] for b in range(self.N) ] )
        
    def _chop(self):
        tol = 1e-14
        self.m[abs(self.m) < tol] = 0.0
        self.m = np.real_if_close(self.m,tol=100)
        return self

    def gate1(self, i, b):
        assert(self.N % 2 == 0)
        P=self._swapmat(0,i)
        o=np.identity(self.N // 2)
        t=np.dot(np.dot(np.transpose(P),np.concatenate(np.concatenate(np.multiply.outer(o,b),-2),-1)),P)
        r=operator(self.nbits, m = np.dot(t,self.m))
        r._chop()
        return r

    def gate2(self, i, j, b):
        assert(i!=j)
        assert(self.N % 4 == 0)

        P=self._swapmat(1,i)
        if j == 1:
            j=i
        P=np.dot(self._swapmat(0,j),P)
        o=np.identity(self.N // 4)
        t=np.dot(np.dot(np.transpose(P),np.concatenate(np.concatenate(np.multiply.outer(o,b),-2),-1)),P)
        r=operator(self.nbits, m = np.dot(t,self.m))
        r._chop()
        return r
