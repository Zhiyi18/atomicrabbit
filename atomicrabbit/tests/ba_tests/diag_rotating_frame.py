"""
Diagnostic: check assign_rotating_frame_energies() for Ba-137.
Run from project root: python -m atomicrabbit.tests.diag_rotating_frame
"""
import numpy as np
from atomicrabbit.physics import states, lasers, transitions, atom

# ── Ba-137 setup (same as ba_137_dynamics_test.py) ───────────────────────────
state1 = states.FineState(configuration='6s2',  L=0, S=0, J=0,
                          energy_cm=0,        parity='',
                          A_hyperfine=-109.5, B_hyperfine=0)
state2 = states.FineState(configuration='6s6p', L=1, S=0, J=1,
                          energy_cm=18060.261, parity='',
                          A_hyperfine=-109.5,  B_hyperfine=0)
state3 = states.FineState(configuration='6s6d', L=2, S=0, J=2,
                          energy_cm=30236.826, parity='',
                          A_hyperfine=-109.5,  B_hyperfine=0)

t12 = transitions.Transition(state1, state2, linewidth=20e6, A=1.19e8)
t23 = transitions.Transition(state2, state3, linewidth=20e6, A=5.00e6)
t13 = transitions.Transition(state1, state3, linewidth=20e6, A=1e5)

laser1 = lasers.Laser(frequency=t12.frequency, intensity=1e5,
                      linewidth=20e6, polarization=1)
laser2 = lasers.Laser(frequency=t23.frequency, intensity=1e5,
                      linewidth=20e6, polarization=1)

Ba137 = atom.Atom(species='Ba I', isotope='137',
                  transition_list=[t12, t23, t13],
                  laser_list=[laser1, laser2],
                  I=3/2, mu_I=0.1, B=1e-4)

import contextlib, io

def silent(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return fn(*args, **kwargs)

# ── 1. atom_graph ─────────────────────────────────────────────────────────────
G = silent(Ba137.atom_graph)
print(f'atom_graph(): {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

# ── 2. assign_rotating_frame_energies ────────────────────────────────────────
energies = silent(Ba137.assign_rotating_frame_energies)
print(f'\nassign_rotating_frame_energies() returned {len(energies)} entries:')

for s, E in sorted(energies.items(), key=lambda x: x[1].real if hasattr(x[1], 'real') else x[1]):
    passes = abs(E) <= 1e10
    label = f'F={s.F} mF={s.mF} config={s.parent.parent.configuration}'
    print(f'  E={E:+.6e} rad/s  |E|<=1e10: {passes}  {label}')

# ── 3. build_hamiltonian diagonal ─────────────────────────────────────────────
H, index = silent(Ba137.build_hamiltonian)

diag = np.diag(H)
n_nonzero = np.sum(diag != 0)
print(f'\nbuild_hamiltonian() diagonal: {n_nonzero}/{len(diag)} nonzero entries')
for s, i in sorted(index.items(), key=lambda x: x[1]):
    d = diag[i]
    if d != 0:
        print(f'  [{i}] {s}  diag={d:+.4e}')
