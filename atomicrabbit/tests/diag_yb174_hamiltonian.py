"""
Test build_hamiltonian() vs manual rotating-frame H for Yb-174 at the 399 nm peak.
Run from project root: python -m atomicrabbit.tests.diag_yb174_hamiltonian
"""
import contextlib, io
import numpy as np
from pathlib import Path
from scipy.constants import c
from qutip import Qobj, steadystate

from atomicrabbit.physics.atom import Atom
from atomicrabbit.physics.lasers import Laser

DB_PATH = str(Path(__file__).parent.parent / "db" / "spectroscopy.db")
INTENSITY = 10.0        # W/m^2
LINEWIDTH = 30e6        # Hz
B_FIELD   = 1e-4        # T  (1 Gauss)

def silent(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return fn(*args, **kwargs)

# ── find resonance frequency from bare atom ───────────────────────────────────
bare = Atom.make_atom('Yb', 174, DB_PATH, laser_list=[], B=B_FIELD)
t399 = next(t for t in bare.transition_list if abs(t.upper.energy_cm - 25068) < 100)
nu0  = (t399.upper.energy_cm - t399.lower.energy_cm) * c * 100
print(f'Yb-174 1S0→1P1 frequency: {nu0:.6e} Hz')

# ── atom with pi-polarization laser exactly on resonance ─────────────────────
laser = Laser(frequency=nu0, intensity=INTENSITY, linewidth=LINEWIDTH, polarization=0)
atom  = Atom.make_atom('Yb', 174, DB_PATH, laser_list=[laser], B=B_FIELD)

zeeman_states, index = silent(atom.basis)
n = len(zeeman_states)
print(f'Basis: {n} states')

is_1P1 = np.array([20000 < s.energy_cm < 30000 for s in zeeman_states])
is_gnd = np.array([s.energy_cm < 1000           for s in zeeman_states])
excited_idx = [i for i, s in enumerate(zeeman_states) if 20000 < s.energy_cm < 30000]
print(f'Ground states: {is_gnd.sum()},  1P1 states: {is_1P1.sum()}')

c_ops = silent(atom.build_collapse_operators)

# ── Method 1: build_hamiltonian() ─────────────────────────────────────────────
H_arr, _ = silent(atom.build_hamiltonian)
rho_ss    = steadystate(Qobj(H_arr), c_ops)
Pe_build  = sum(rho_ss.full()[i, i].real for i in excited_idx)
print(f'\nMethod 1  build_hamiltonian():  P_excited = {Pe_build:.6f}')

# ── Method 2: manual rotating-frame H (polarization-scan approach) ────────────
driven = silent(atom.driven_transitions)
couplings = [(index[d.transition.lower], index[d.transition.upper], d.rabi_frequency / 2)
             for d in driven]

state_E_omega = np.array([s.energy_omega for s in zeeman_states])
omega_L = 2 * np.pi * nu0

H_offdiag = np.zeros((n, n), dtype=complex)
for gi, ei, omega_half in couplings:
    H_offdiag[gi, ei] += omega_half
    H_offdiag[ei, gi] += np.conj(omega_half)

diag = np.where(is_1P1, state_E_omega - omega_L,
       np.where(is_gnd, state_E_omega, 0.0))
H_manual = H_offdiag.copy()
np.fill_diagonal(H_manual, diag)

rho_ss2  = steadystate(Qobj(H_manual), c_ops)
Pe_manual = sum(rho_ss2.full()[i, i].real for i in excited_idx)
print(f'Method 2  manual rotating frame: P_excited = {Pe_manual:.6f}')

print(f'\nDifference: {abs(Pe_build - Pe_manual):.2e}')
