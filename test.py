#!/usr/bin/env python
import sqc

Nbits=3

sqc.seed(13)

s=sqc.state(Nbits)
s=s.clone()

o=sqc.operator(Nbits)

#o=o.H(1)
o=o.u3(1,0.5,0.222,0.315)

print o.nonunitarity()
print s, "\n--"
s=o * s

print s, "\n--"

for n in range(10):
    sp,v1=s.measure(1)
    print sp, v1, "\n\n"


