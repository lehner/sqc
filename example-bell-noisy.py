#!/usr/bin/env python
import sqc

Nbits=2

# Create Bell state
op=sqc.operator(Nbits).H(0).CNOT(0,1)

# Create noise model
nm=sqc.noise.model.simple(
    T1 = 10000,
    gate_times = { "X" : 0.2, "CNOT" : 0.5, "Rz" : 0.7, "H": 0.3 }, # times are in units of mus
    qubit_readout_errors = [ 0.1, 0.0 ],
    gate_depolarization_p = { "X" : 0.0, "CNOT" : 0.0, "Rz" : 0.0, "H" : 0.0 }
)

# Initial state
s0=sqc.state(Nbits)

# Sample
res=sqc.noise.sample(nm,op,s0,100)
print(res)
