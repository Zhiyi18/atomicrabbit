"""
Yb_dynamics_test.py
Population dynamics for Yb-171 and Yb-173 driven on the 399 nm 1S0 -> 1P1 transition.

Two panels per isotope:
  - Ground manifold (1S0): population of each |F, mF> state vs time
  - Excited manifold (1P1): population of each |F, mF> state vs time

Initial state: ground |F, mF=-F> (most negative Zeeman sublevel).
Laser: sigma+ (q=+1), resonant with the bare 1S0->1P1 fine transition.

Run from the project root:
  python -m atomicrabbit.tests.Yb_dynamics_test
"""

import dataclasses
import io
import contextlib
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import qutip as qt
from scipy.constants import c

from atomicrabbit.physics.atom import Atom
from atomicrabbit.physics.lasers import Laser

DB_PATH = str(Path(__file__).parent.parent / "db" / "spectroscopy.db")

B_FIELD        = 0.7e-4   # Tesla (0.7 G)
INTENSITY_W_m2 = 15.92    # W/m²
LINEWIDTH_HZ   = 28.9e6     # Hz
POLARIZATION   = 1        # sigma+

T_END   = 1e-4   # s
N_STEPS = 2000

OUT_DIR = Path(__file__).parent


def find_399_transition(atom):
    return next(
        (t for t in atom.transition_list if abs(t.upper.energy_cm - 25068) < 100),
        None,
    )


def fmt_state(s):
    """Format |F, mF> label, handling half-integer quantum numbers."""
    F2  = round(2 * s.F)
    mF2 = round(2 * s.mF)
    sign = '+' if mF2 >= 0 else ''
    if F2 % 2 == 0:
        return f"F={F2//2}, m$_F$={sign}{mF2//2}"
    else:
        return f"F={F2}/2, m$_F$={sign}{mF2}/2"


def fmt_state_csv(s):
    """Plain-text |F, mF> column name for CSV headers, tagged by manifold."""
    F2  = round(2 * s.F)
    mF2 = round(2 * s.mF)
    sign = 'p' if mF2 >= 0 else 'm'
    F_str  = f"{F2//2}" if F2 % 2 == 0 else f"{F2}-2"
    mF_str = f"{sign}{abs(mF2)//2}" if mF2 % 2 == 0 else f"{sign}{abs(mF2)}-2"
    if s.energy_cm < 1000:
        manifold = "gnd"
    elif 20000 < s.energy_cm < 30000:
        manifold = "exc"
    else:
        manifold = "other"
    return f"{manifold}_F{F_str}_mF{mF_str}"


def save_csv(A, zeeman_states, index, result, tlist):
    ordered_states = sorted(zeeman_states, key=lambda s: (s.energy_cm, s.F, s.mF))
    header = ["time_ns"] + [fmt_state_csv(s) for s in ordered_states]
    data = np.column_stack(
        [tlist * 1e9] + [result.expect[index[s]] for s in ordered_states]
    )
    outpath = OUT_DIR / f"Yb{A}_dynamics.csv"
    np.savetxt(outpath, data, delimiter=",", header=",".join(header),
               comments="", fmt="%.8e")
    print(f"  Saved -> {outpath}")


def run_dynamics(A):
    print(f"\n{'='*60}")
    print(f"  Yb-{A}  B={B_FIELD*1e4:.1f} G  I={INTENSITY_W_m2} W/m^2  q={POLARIZATION}")
    print(f"{'='*60}")

    atom0 = Atom.make_atom('Yb', A, DB_PATH, laser_list=[], B=B_FIELD)
    t399  = find_399_transition(atom0)
    nu399 = (t399.upper.energy_cm - t399.lower.energy_cm) * c * 100

    laser = Laser(
        frequency    = nu399,
        intensity    = INTENSITY_W_m2,
        linewidth    = LINEWIDTH_HZ,
        polarization = POLARIZATION,
    )
    atom = dataclasses.replace(atom0, laser_list=[laser])

    zeeman_states, index = atom.basis()
    n = len(zeeman_states)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        H_arr, _ = atom.build_hamiltonian()
        c_ops    = atom.build_collapse_operators()

    n_driven = int(np.count_nonzero(np.triu(H_arr, 1)))
    print(f"  {n} basis states, {n_driven} driven Zeeman transitions")
    H = qt.Qobj(H_arr)

    L = qt.liouvillian(H, c_ops)
    rates = np.sort(np.abs(L.eigenenergies()))
    rates = rates[rates > 1e-6]
    print(f"  Liouvillian spectrum: |eig| in [{rates.min():.3e}, {rates.max():.3e}] rad/s"
          f"  (span factor {rates.max()/rates.min():.1e})")

    # Initial state: ground level, lowest mF
    gnd_states = [(s, index[s]) for s in zeeman_states if s.energy_cm < 1000]
    psi0_state, psi0_idx = min(gnd_states, key=lambda x: x[0].mF)
    print(f"  Initial state: {fmt_state(psi0_state)}")

    rho0  = qt.ket2dm(qt.basis(n, psi0_idx))
    tlist = np.linspace(0, T_END, N_STEPS)
    e_ops = [qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]

    print(f"  Running mesolve  ({N_STEPS} steps, t_end={T_END*1e9:.0f} ns)...")
    result = qt.mesolve(H, rho0, tlist, c_ops, e_ops=e_ops,
                        options={"method": "diag"})
    print("  Done.")

    save_csv(A, zeeman_states, index, result, tlist)

    return zeeman_states, index, result, tlist


def plot_dynamics(A, zeeman_states, index, result, tlist):
    gnd_states = sorted(
        [s for s in zeeman_states if s.energy_cm < 1000],
        key=lambda s: (s.F, s.mF)
    )
    exc_states = sorted(
        [s for s in zeeman_states if 20000 < s.energy_cm < 30000],
        key=lambda s: (s.F, s.mF)
    )

    panels = [
        (f"$^{{{A}}}$Yb  Ground  $^1S_0$", gnd_states, (0, 1)),
        (f"$^{{{A}}}$Yb  Excited $^1P_1$", exc_states, (0, 0.1)),
    ]

    colors = plt.cm.tab20.colors
    ls_by_F = {}  # assign linestyle per F value
    linestyles = ['-', '--', '-.', ':']

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(
        f"$^{{{A}}}$Yb population dynamics — 399 nm $\\sigma^+$, "
        f"B={B_FIELD*1e4:.1f} G, I={INTENSITY_W_m2} W/m²",
        fontsize=12
    )

    for ax, (title, states, ylim) in zip(axes, panels):
        # Assign a linestyle per F group
        F_vals = sorted(set(s.F for s in states))
        for fi, Fv in enumerate(F_vals):
            ls_by_F[Fv] = linestyles[fi % len(linestyles)]

        # Color cycles within each F group
        color_offset = 0
        for Fv in F_vals:
            group = [s for s in states if s.F == Fv]
            for j, s in enumerate(group):
                idx = index[s]
                ax.plot(tlist * 1e9, result.expect[idx],
                        label=fmt_state(s),
                        color=colors[(color_offset + j) % 20],
                        ls=ls_by_F[Fv], lw=1.2)
            color_offset += len(group)

        ax.set_ylabel("Population", fontsize=10)
        ax.set_title(title, fontsize=11)
        ncol = max(1, len(states) // 6)
        ax.legend(fontsize=7.5, ncol=ncol, loc='upper right')
        ax.set_ylim(*ylim)
        ax.grid(True, which='major', alpha=0.25)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    axes[-1].set_xlabel("Time (ns)", fontsize=11)
    fig.tight_layout()

    outpath = OUT_DIR / f"Yb{A}_dynamics.png"
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"  Saved -> {outpath}")
    plt.show()


# ── main ─────────────────────────────────────────────────────────────────────

for A in [171, 173]:
    zeeman_states, index, result, tlist = run_dynamics(A)

#%%
for A in [171, 173]:
    plot_dynamics(A, zeeman_states, index, result, tlist)
    
print("\nAll done.")
