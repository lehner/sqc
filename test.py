#!/usr/bin/env python
import sqc

Nbits=3

s=sqc.state(Nbits)
s=s.clone()

o=sqc.operator(Nbits)

o=o.H(2).H(1)

print s
s=o * s

print s
