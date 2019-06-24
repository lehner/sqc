#!/usr/bin/env python
import sqc

Nbits=2

# Create Bell state
op=sqc.operator(Nbits).H(0).CNOT(0,1)

# Create noise model
nm=sqc.noise.model.simple(
    T1 = 50,
    gate_times = { "X" : 0.5, "CNOT" : 0.9, "Rz" : 0.5, "H": 0.5 },
    qubit_readout_errors = [ 0.04, 0.05 ],
    gate_depolarization_p = { "X" : 0.04, "CNOT" : 0.06, "Rz" : 0.0, "H" : 0.03 }
)

# Initial state
s0=sqc.state(Nbits)

# Sample
res=sqc.noise.sample(nm,op,s0,1000)
print(res)
