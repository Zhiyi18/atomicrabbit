"""
Yb_polarization_scan.py
Computes the steady-state 399 nm absorption spectrum for all stable Yb isotopes
at three laser polarizations (sigma+, pi, sigma-) and saves the results as CSV.

Parameters: B = 0.7 G, I = 15.92 W/m²

Output columns per CSV:
  freq_offset_MHz   -- laser frequency offset from ^174Yb resonance
  Pe_Yb<A>          -- steady-state excited-state population for each isotope (raw)
  total             -- natural-abundance-weighted sum: sum_A(abundance_A * Pe_A)

Run from the project root:
  python -m atomicrabbit.tests.yb_B07G_I1592.Yb_polarization_scan
"""

import io
import contextlib
import dataclasses
import csv
from pathlib import Path

import numpy as np
from scipy.constants import c
from qutip import Qobj, steadystate

from atomicrabbit.physics.atom import Atom
from atomicrabbit.physics.lasers import Laser

DB_PATH    = str(Path(__file__).parent.parent.parent / "db" / "spectroscopy.db")
OUT_DIR    = Path(__file__).parent
YB_ISOTOPES = [168, 170, 171, 172, 173, 174, 176]

INTENSITY_W_m2 = 15.92
LINEWIDTH_HZ   = 30e6
B_FIELD        = 0.7e-4       # 0.7 Gauss

ABUNDANCES = {168: 0.0013, 170: 0.0304, 171: 0.1428, 172: 0.2183,
              173: 0.1613, 174: 0.3183, 176: 0.1276}

N_POINTS     = 250
FREQ_LOW_HZ  = -1.2e9
FREQ_HIGH_HZ =  2.5e9

POL_NAMES = {-1: 'sigma_minus', 0: 'pi', 1: 'sigma_plus'}


def find_399_transition(atom):
    return next(
        (t for t in atom.transition_list if abs(t.upper.energy_cm - 25068) < 100),
        None,
    )


def compute_spectrum(polarization):
    """Return (freq_offsets_Hz, {A: excited_pop array}) for one polarization."""

    ref_atom = Atom.make_atom('Yb', 174, DB_PATH, laser_list=[], B=B_FIELD)
    ref_t    = find_399_transition(ref_atom)
    nu_ref   = (ref_t.upper.energy_cm - ref_t.lower.energy_cm) * c * 100

    freq_offsets = np.linspace(FREQ_LOW_HZ, FREQ_HIGH_HZ, N_POINTS)
    freq_scan    = nu_ref + freq_offsets

    laser_ref = Laser(
        frequency   = nu_ref,
        intensity   = INTENSITY_W_m2,
        linewidth   = LINEWIDTH_HZ,
        polarization= polarization,
    )

    results = {}

    for A in YB_ISOTOPES:
        atom0    = Atom.make_atom('Yb', A, DB_PATH, laser_list=[], B=B_FIELD)
        atom_ref = dataclasses.replace(atom0, laser_list=[laser_ref])

        zeeman_states, index = atom0.basis()
        n = len(zeeman_states)
        excited_idx = [i for i, s in enumerate(zeeman_states) if 20000 < s.energy_cm < 30000]

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            c_ops    = atom_ref.build_collapse_operators()
            H_ref, _ = atom_ref.build_hamiltonian()

        diag_ref  = np.diag(H_ref).copy()
        H_offdiag = H_ref.copy()
        np.fill_diagonal(H_offdiag, 0)

        is_rf_1P1 = np.array([
            20000 < zeeman_states[i].energy_cm < 30000
            for i in range(n)
        ])

        print(f'  Yb-{A}: {n} states, {is_rf_1P1.sum()} driven 1P1 states', flush=True)

        excited_pop = np.zeros(N_POINTS)

        for k, nu_L in enumerate(freq_scan):
            diag_k = diag_ref.copy()
            diag_k[is_rf_1P1] += 2 * np.pi * (nu_ref - nu_L)
            H_arr = H_offdiag.copy()
            np.fill_diagonal(H_arr, diag_k)

            rho_ss      = steadystate(Qobj(H_arr), c_ops)
            rho_full    = rho_ss.full()
            excited_pop[k] = sum(rho_full[i, i].real for i in excited_idx)

            if (k + 1) % 50 == 0:
                print(f'    Yb-{A}: {k+1}/{N_POINTS} done', flush=True)

        results[A] = excited_pop
        print(f'    peak P_excited = {excited_pop.max():.4f}', flush=True)

    return freq_offsets, results


def save_csv(polarization, freq_offsets, results):
    pol_name = POL_NAMES[polarization]
    outpath  = OUT_DIR / f'Yb_spectrum_{pol_name}.csv'

    total = sum(ABUNDANCES[A] * results[A] for A in YB_ISOTOPES)

    with open(outpath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        f.write(f'# Yb 399 nm absorption spectrum  q={polarization} ({pol_name})\n')
        f.write(f'# B={B_FIELD*1e4:.2f} G   I={INTENSITY_W_m2} W/m^2   N_points={N_POINTS}\n')
        f.write(f'# Pe_Yb<A>: raw steady-state excited population for each isotope\n')
        f.write(f'# total: natural-abundance-weighted sum\n')

        header = (['freq_offset_MHz']
                  + [f'Pe_Yb{A}' for A in YB_ISOTOPES]
                  + ['total'])
        writer.writerow(header)

        for k in range(N_POINTS):
            row = ([f'{freq_offsets[k]/1e6:.4f}']
                   + [f'{results[A][k]:.6f}' for A in YB_ISOTOPES]
                   + [f'{total[k]:.6f}'])
            writer.writerow(row)

    print(f'  Saved -> {outpath}', flush=True)


# ── main ──────────────────────────────────────────────────────────────────────

for q in [1, 0, -1]:
    print()
    print('=' * 60)
    print(f'  Polarization q={q}  ({POL_NAMES[q]})')
    print('=' * 60)
    freq_offsets, results = compute_spectrum(q)
    save_csv(q, freq_offsets, results)

print()
print('All done.')
