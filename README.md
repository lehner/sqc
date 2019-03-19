# sqc
This repository contains a simple digital quantum computer
simulator.  Its main purpose is to illustrate concepts
introduced in an accompanying lecture on quantum computing.

As an example, creating an entangled Bell state and
performing ten classical measurements can be done
by the simple code

```
import sqc

Nbits=2

# Create Bell state
s=sqc.operator(Nbits).H(0).CNOT(0,1) * sqc.state(Nbits)

# Print state
print s

# Perform 10 measurements
for n in range(10):
    s0,v0=s.measure(0)
    s1,v1=s0.measure(1)
    
    # And print results of state and classical bits
    print s1, v1, v0

```

