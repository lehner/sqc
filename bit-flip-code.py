#!/usr/bin/env python
import sqc

Nbits=5

# Circuit
op=sqc.operator(Nbits)

def ERR(op):
    return op.X(2)

op=ERR(op.CNOT(0,1).CNOT(0,2)).CNOT(0,3).CNOT(1,3).CNOT(0,4).CNOT(2,4).M(3,0).M(4,1)

# If
op=op.IF(1).X(1).X(3).ENDIF()

op=op.IF(2).X(2).X(4).ENDIF()

op=op.IF(3).X(0).X(3).X(4).ENDIF()

print(op.toQASM())

# Initial state
s0=sqc.operator(Nbits).H(0).Rz(0,0.5)*sqc.state(Nbits)
print(s0)

s=op*s0

print(s)
