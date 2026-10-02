# -*- coding: utf-8 -*-
"""
Yb_spectrum_test.py
Tests Atom.make_atom for all stable Yb isotopes, builds Hamiltonians,
and drives the 399 nm (1S0 -> 1P1) transition with a resonant laser.
"""

import io
import contextlib
from pathlib import Path
import numpy as np
from atomicrabbit.physics.atom import Atom
from atomicrabbit.physics.lasers import Laser

DB_PATH = str(Path(__file__).parent.parent / "db" / "spectroscopy.db")

YB_ISOTOPES = [168, 170, 171, 172, 173, 174, 176]

# 399 nm laser parameters
INTENSITY_W_m2 = 10.0       # W/m^2 (weak probe)
LINEWIDTH_HZ   = 30e6       # ~natural linewidth of 1P1
POLARIZATION   = 0          # pi (drives Delta_mF = 0)
B_FIELD        = 1e-4       # Tesla (1 Gauss -- resolves Zeeman structure)


def find_399_transition(atom):
    """Return the FineTransition whose upper level is the 1P1 state (~25068 cm^-1)."""
    return next(
        (t for t in atom.transition_list
         if abs(t.upper.energy_cm - 25068) < 100),
        None,
    )


# ── Section 1: level structure for all isotopes ───────────────────────────────

print("=" * 60)
print("  SECTION 1 -- Level structure (no laser)")
print("=" * 60)

summary = {}

for A in YB_ISOTOPES:
    atom = Atom.make_atom('Yb', A, DB_PATH, laser_list=[], B=0)

    print(f"\nYb-{A}   I={atom.I}   mu_I={atom.mu_I}")
    for t in atom.transition_list:
        dE = t.upper.energy_cm - t.lower.energy_cm
        wl = 1e7 / dE
        print(f"  {t.lower.configuration} J={t.lower.J:.0f}"
              f" -> {t.upper.configuration} J={t.upper.J:.0f}"
              f"   {wl:.1f} nm   A={t.A:.3e} s^-1")
        if t.upper.A_hyperfine is not None:
            print(f"    upper A_hf={t.upper.A_hyperfine:.3f} MHz"
                  f"  B_hf={t.upper.B_hyperfine} MHz")

    zeeman_states, _ = atom.basis()
    summary[A] = {'atom': atom, 'n_basis': len(zeeman_states)}
    print(f"  Basis: {len(zeeman_states)} Zeeman states")


# ── Section 2: 399 nm laser test ──────────────────────────────────────────────

print()
print("=" * 60)
print("  SECTION 2 -- 399 nm driven transitions")
print(f"  I={INTENSITY_W_m2} W/m^2   q={POLARIZATION}   B={B_FIELD*1e4:.0f} G")
print("=" * 60)

import dataclasses

for A in YB_ISOTOPES:
    # Load atom without laser to retrieve the isotope-specific 399 nm frequency
    atom0 = Atom.make_atom('Yb', A, DB_PATH, laser_list=[], B=B_FIELD)
    t399  = find_399_transition(atom0)
    if t399 is None:
        print(f"\nYb-{A}: no 399 nm transition in DB -- skipping")
        continue

    from scipy.constants import c
    nu399 = (t399.upper.energy_cm - t399.lower.energy_cm) * c * 100

    laser = Laser(
        frequency   = nu399,
        intensity   = INTENSITY_W_m2,
        linewidth   = LINEWIDTH_HZ,
        polarization= POLARIZATION,
    )

    # Rebuild with laser (frozen dataclass -- use replace)
    atom = dataclasses.replace(atom0, laser_list=[laser])

    print(f"\nYb-{A}   nu={nu399/1e12:.4f} THz")

    # Capture verbose debug prints from driven_transitions / build_hamiltonian
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        driven = atom.driven_transitions()
        H, index = atom.build_hamiltonian()

    print(f"  Driven transitions: {len(driven)}")
    for d in driven:
        omega_MHz = abs(d.rabi_frequency) / (2 * np.pi * 1e6)
        det_MHz   = d.detuning / 1e6  # detuning is in Hz, not rad/s
        print(f"    {d.lower.spectroscopic_name:35s}"
              f" -> {d.upper.spectroscopic_name:35s}"
              f"   Omega/2pi={omega_MHz:7.3f} MHz   delta/2pi={det_MHz:+.3f} MHz")

    offdiag = H.copy()
    np.fill_diagonal(offdiag, 0)
    print(f"  H shape: {H.shape}"
          f"   ||H_offdiag|| = {np.linalg.norm(offdiag):.4e} rad/s")


# ── Summary ───────────────────────────────────────────────────────────────────

print()
print("=" * 60)
print("  Summary")
print("=" * 60)
print(f"  {'Isotope':>8}  {'I':>4}  {'Basis':>6}  {'Transitions':>12}")
for A, d in summary.items():
    atom = d['atom']
    print(f"  Yb-{A:<3d}   {atom.I:>4.1f}  {d['n_basis']:>6}  {len(atom.transition_list):>12}")


# ── Section 3: 399 nm absorption spectrum (steady-state OBE) ──────────────────

print()
print("=" * 60)
print("  SECTION 3 -- 399 nm absorption spectrum (steady-state OBE)")
print(f"  q={POLARIZATION}   B={B_FIELD*1e4:.0f} G   I={INTENSITY_W_m2} W/m^2")
print("=" * 60)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qutip import Qobj, steadystate

# Natural abundances (IUPAC 2021)
ABUNDANCES = {168: 0.0013, 170: 0.0304, 171: 0.1428, 172: 0.2183,
              173: 0.1613, 174: 0.3183, 176: 0.1276}

COLORS = {168: '#888888', 170: '#4477AA', 171: '#EE6677',
          172: '#228833', 173: '#CCBB44', 174: '#AA3377', 176: '#66CCEE'}

# Reference frequency: Yb-174 1P1 bare level
ref_atom = Atom.make_atom('Yb', 174, DB_PATH, laser_list=[], B=B_FIELD)
ref_t399 = find_399_transition(ref_atom)
nu_ref   = (ref_t399.upper.energy_cm - ref_t399.lower.energy_cm) * c * 100

# Frequency scan: -1.2 to +2.5 GHz around Yb-174, 500 points
# (~7 MHz/point -- resolves 30 MHz natural linewidth with ~4 pts per FWHM)
N_POINTS     = 500
freq_offsets = np.linspace(-1.2e9, 2.5e9, N_POINTS)   # Hz
freq_scan    = nu_ref + freq_offsets

# Reference laser used once per isotope to evaluate Rabi couplings.
# Detuning filter in driven_transitions() is +/-1e13 Hz; all 399 nm components
# lie within 3 GHz of nu_ref so all are captured regardless of isotope.
laser_ref = Laser(
    frequency   = nu_ref,
    intensity   = INTENSITY_W_m2,
    linewidth   = LINEWIDTH_HZ,
    polarization= POLARIZATION,
)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1]})

total_signal = np.zeros(N_POINTS)

for A in YB_ISOTOPES:
    atom0 = Atom.make_atom('Yb', A, DB_PATH, laser_list=[], B=B_FIELD)
    zeeman_states, index = atom0.basis()
    n = len(zeeman_states)

    # 1P1 excited states: energy_cm in (20000, 30000) cm^-1
    # (1S0 ~ 0, 3P1 ~ 17993, 1P1 ~ 25068)
    excited_idx = [index[s] for s in zeeman_states if 20000 < s.energy_cm < 30000]

    # Pre-compute per-state rotating-frame energy offsets.
    # For 1P1 states: E_rot = energy_omega - omega_L  (omega_L varies per scan point)
    # For 1S0 states: E_rot = energy_omega            (ground; value is ~0 + tiny Zeeman)
    # For 3P1 states: set to 0                        (decoupled from 399 nm laser)
    state_E_omega = np.array([s.energy_omega for s in zeeman_states])
    is_1P1  = np.array([20000 < s.energy_cm < 30000 for s in zeeman_states])
    is_gnd  = np.array([s.energy_cm < 1000           for s in zeeman_states])
    # states that are neither 1P1 nor ground (i.e. 3P1): diagonal stays 0

    # Collapse operators and Rabi couplings: computed once, laser-freq-independent.
    atom_ref = dataclasses.replace(atom0, laser_list=[laser_ref])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        c_ops  = atom_ref.build_collapse_operators()
        driven = atom_ref.driven_transitions()
        # rabi_frequency already divides by hbar; off-diagonal entry is Omega/2
        couplings = [
            (index[d.transition.lower],
             index[d.transition.upper],
             d.rabi_frequency / 2)
            for d in driven
        ]

    print(f'  Yb-{A}: {n} states, {len(couplings)} driven Zeeman transitions', flush=True)

    # Build the off-diagonal coupling block (constant across scan)
    H_offdiag = np.zeros((n, n), dtype=complex)
    for gi, ei, omega_half in couplings:
        H_offdiag[gi, ei] += omega_half
        H_offdiag[ei, gi] += np.conj(omega_half)

    excited_pop = np.zeros(N_POINTS)

    for k, nu_L in enumerate(freq_scan):
        omega_L = 2 * np.pi * nu_L

        # Rotating-frame diagonal: ground keeps energy_omega (~0),
        # 1P1 gets energy_omega - omega_L (= detuning in rad/s),
        # 3P1 stays 0 (decoupled).
        diag = np.where(is_1P1, state_E_omega - omega_L,
               np.where(is_gnd, state_E_omega, 0.0))

        H_arr = H_offdiag.copy()
        np.fill_diagonal(H_arr, diag)

        rho_ss = steadystate(Qobj(H_arr), c_ops)
        rho_full = rho_ss.full()
        excited_pop[k] = sum(rho_full[i, i].real for i in excited_idx)

        if (k + 1) % 100 == 0:
            print(f'    Yb-{A}: {k+1}/{N_POINTS} points done', flush=True)

    abundance = ABUNDANCES.get(A, 0)
    # Plot absolute excited population scaled by abundance (correct relative heights)
    ax1.plot(freq_offsets / 1e6, abundance * excited_pop,
             color=COLORS[A], linewidth=1.5, alpha=0.85,
             label=f'$^{{{A}}}$Yb  ({abundance*100:.2f}%)')

    total_signal += abundance * excited_pop
    print(f'    peak P_excited = {excited_pop.max():.4f}')

# Natural-abundance-weighted total
norm_tot = total_signal.max() if total_signal.max() > 0 else 1.0
ax2.fill_between(freq_offsets / 1e6, total_signal / norm_tot, color='k', alpha=0.25)
ax2.plot(freq_offsets / 1e6, total_signal / norm_tot, 'k-', linewidth=1.2)

ax1.set_ylabel('Excited state population x abundance', fontsize=10)
ax1.legend(fontsize=9, ncol=2, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_title(
    'Yb 399 nm ($^1S_0 \\to\\ ^1P_1$) absorption spectrum  '
    f'[steady-state OBE,  $q={POLARIZATION}$,  $B={B_FIELD*1e4:.0f}$ G,  '
    f'$I={INTENSITY_W_m2}$ W/m$^2$]',
    fontsize=11,
)

ax2.set_xlabel('Laser frequency offset from $^{174}$Yb [MHz]', fontsize=11)
ax2.set_ylabel('Total signal\n(natural abund.)', fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(bottom=0)

fig.tight_layout()
outpath = str(Path(__file__).parent / 'Yb_399nm_spectrum.png')
fig.savefig(outpath, dpi=150, bbox_inches='tight')
print(f'\n  Spectrum saved -> {outpath}')
plt.close(fig)
