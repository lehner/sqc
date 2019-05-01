#!/usr/bin/env python
import sqc

Nbits=2

def Uf_const(o): # N=1, f(x)=1
    return o.NOT(1)

def Uf_balanced(o): # N=1, f(x)=x
    return o.CNOT(0,1)

for Uf,n in [ (Uf_const,"f(x)=1"), (Uf_balanced,"f(x)=x") ]:
    print("Run Deutsch-Jozsa for %s" % n)
    print("--------------------------------------")

    # Psi0
    psi0=sqc.operator(Nbits).X(1) * sqc.state(Nbits, basis=[ "|%d>|%d>" % (i%2,i//2) for i in range(4) ])
    print("|Psi0> = ")
    print(psi0)

    # Psi1
    psi1=sqc.operator(Nbits).H(0).H(1) * psi0
    print("|Psi1> = ")
    print(psi1)

    # Psi2
    psi2=Uf(sqc.operator(Nbits)) * psi1
    print("|Psi2> = ")
    print(psi2)

    # Psi3
    psi3=sqc.operator(Nbits).H(0) * psi2
    print("|Psi3> = ")
    print(psi3)

    # Measurement
    psi4,r=psi3.measure(0)
    print("Result: %d, %s is %s" % (r,n,["constant","balanced"][r]))

    print("")
