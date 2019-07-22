#!/usr/bin/env python
import sqc

Nbits=3

# Create Greenberger-Horne-Zeilinger state
op=sqc.operator(Nbits).H(0).CNOT(0,1).CNOT(0,2)
print(op.toQASM())
s=op * sqc.state(Nbits)

# Print state
print("GHZ state")
print(s)

# Perform 10 measurements
for n in range(10):
    s0,v0=s.measure(0)
    s1,v1=s0.measure(1)
    s2,v2=s1.measure(2)
    
    # And print results of state and classical bits
    print("Measurement # %d: %s | %d%d%d" % (n,s2,v2,v1,v0))
