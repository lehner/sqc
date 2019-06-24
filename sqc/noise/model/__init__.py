#
# Christoph Lehner 2019
#
from sqc.noise.model.generic import generic
from sqc.state import state
import numpy as np

def one_qubit_diag(i, s, zerofac, onefac, o):
    n=s.nbits
    return state(n,(o.onemask[i]*onefac + o.zeromask[i]*zerofac)*s.v,basis=s.basis)

def one_qubit_odiag(i, s, zerofac, onefac, o):
    n=s.nbits
    return state(n,(o.onemask[i]*onefac + o.zeromask[i]*zerofac)*s.v[o.notmask[i]],basis=s.basis)
    
def kraus_one_qubit(Psik):
    nk=[ np.dot(np.conj(np.transpose(Psik[i].v)),Psik[i].v) for i in range(len(Psik)) ]
    eps=1 - sum(nk)
    assert(abs(eps) < 1e-13)
    x=np.random.uniform()
    p0=0.0
    for i in range(len(Psik)):
        p1=p0 + nk[i]
        if x < p1:
            break
        p0=p1
    f=1.0 / np.sqrt(nk[i])
    return state(Psik[i].nbits, Psik[i].v * f,basis=Psik[i].basis)

def qubit_relax(i, tOT1, s, o):
    #print(one_qubit_odiag(i, s, 1, 0, o))
    #print(one_qubit_odiag(i, s, 0, 1, o))
    #assert(0)
    gam=1 - np.exp(-tOT1)
    return kraus_one_qubit(
        [ one_qubit_diag(i, s, 1, np.sqrt(1-gam), o),
          one_qubit_odiag(i, s, np.sqrt(gam), 0, o) ])

def qubit_depol(i, p, s, o):
    fac_unit = 1/2*np.sqrt(1 + 3*p)
    fac_z    = 1/2*np.sqrt(1 - p)
    fac_x    = 1/2*np.sqrt(1 - p)
    fac_y    = 1/2*np.sqrt(1 - p)

    return kraus_one_qubit(
        [ one_qubit_diag(i, s, fac_unit, fac_unit, o),
          one_qubit_diag(i, s, fac_z, -fac_z, o),
          one_qubit_odiag(i, s, fac_x, fac_x, o),
          one_qubit_odiag(i, s, -1j*fac_y, 1j*fac_y, o) ] )    

def qubit_depol_two(i, j, p, s, o):
    # for now just depol two qubits individually
    return qubit_depol(i, p, qubit_depol(j, p, s, o), o)

def simple(qubit_readout_errors, T1, gate_times, gate_depolarization_p):
    ng = {
        "H" : lambda f, p, s, o: qubit_depol(p[0], 1-gate_depolarization_p["H"], qubit_relax(p[0],gate_times["H"] / T1, f(*p,s,o), o), o),
        "X" : lambda f, p, s, o: qubit_depol(p[0], 1-gate_depolarization_p["X"], qubit_relax(p[0],gate_times["X"] / T1, f(*p,s,o), o), o),
        "Rz" : lambda f, p, s, o: qubit_depol(p[0], 1-gate_depolarization_p["Rz"], qubit_relax(p[0],gate_times["Rz"] / T1, f(*p,s,o), o), o),
        "CNOT" : lambda f, p, s, o: qubit_depol_two(p[0],p[1], 1-gate_depolarization_p["CNOT"], qubit_relax(p[1],gate_times["CNOT"] / T1, qubit_relax(p[0],gate_times["CNOT"] / T1, f(*p,s,o), o), o), o),
        }
    return generic(noisy_gates = ng, readout_errors = qubit_readout_errors)
