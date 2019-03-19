#!/usr/bin/env python
import sqc

Nbits=2

# Create Bell state
s=sqc.operator(Nbits).H(0).CNOT(0,1) * sqc.state(Nbits)

# Print state
print "Bell state\n", s

# Perform 10 measurements
for n in range(10):
    s0,v0=s.measure(0)
    s1,v1=s0.measure(1)
    
    # And print results of state and classical bits
    print "Measurement #", n, " ", s1, v1, v0
