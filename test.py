#!/usr/bin/env python
import sqc, sys

Nbits=3

sqc.seed(13)

s=sqc.operator(Nbits).H(0).CNOT(0,1) * sqc.state(Nbits)

print s

for n in range(20):
    sp0,v0=s.measure(0)
    sp1,v1=sp0.measure(1)
    sp2,v2=sp1.measure(2)
    
    print sp2, v2,v1,v0, "\n\n"


