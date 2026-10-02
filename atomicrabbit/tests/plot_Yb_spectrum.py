"""
plot_Yb_spectrum.py
Plot the Yb 399 nm absorption spectrum for sigma+, pi, sigma- polarizations.

Run from the project root:
  python -m atomicrabbit.tests.plot_Yb_spectrum
"""

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

DIR = Path(__file__).parent

POLARIZATIONS = [
    ("sigma_plus",  r"$\sigma^+$",  "C0"),
    ("pi",          r"$\pi$",       "C1"),
    ("sigma_minus", r"$\sigma^-$",  "C2"),
]

ISOTOPES   = [168, 170, 171, 172, 173, 174, 176]
ABUNDANCES = {168: 0.0013, 170: 0.0304, 171: 0.1428, 172: 0.2183,
              173: 0.1613, 174: 0.3183, 176: 0.1276}

ODD  = {171, 173}
EVEN = {168, 170, 172, 174, 176}


def load_csv(pol_name):
    path = DIR / f"Yb_spectrum_{pol_name}.csv"
    freq, total, per_isotope = [], [], {A: [] for A in ISOTOPES}
    with open(path, newline='') as f:
        for line in f:
            if line.startswith('#'):
                continue
            break
        reader = csv.DictReader([line] + list(f))
        for row in reader:
            freq.append(float(row['freq_offset_MHz']))
            total.append(float(row['total']))
            for A in ISOTOPES:
                per_isotope[A].append(float(row[f'Pe_Yb{A}']))
    return np.array(freq), np.array(total), {A: np.array(v) for A, v in per_isotope.items()}


fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
fig.suptitle("Yb 399 nm absorption spectrum  (B = 1 G,  I = 10 W/m²)", fontsize=13)

for ax, (pol_name, pol_label, color) in zip(axes, POLARIZATIONS):
    freq, total, per_isotope = load_csv(pol_name)

    # natural-abundance-weighted total
    ax.plot(freq, total, color=color, lw=1.8, label=f"total  {pol_label}", zorder=3)

    # individual isotopes — thin lines
    for A in ISOTOPES:
        w = ABUNDANCES[A]
        ls = '--' if A in ODD else '-'
        ax.plot(freq, w * per_isotope[A], lw=0.7, ls=ls, color='0.55', zorder=2)

    # label the dominant even-isotope peaks
    for A in EVEN:
        w = ABUNDANCES[A]
        trace = w * per_isotope[A]
        peak_idx = np.argmax(trace)
        if trace[peak_idx] > 0.002:
            ax.annotate(f"$^{{{A}}}$Yb",
                        xy=(freq[peak_idx], trace[peak_idx]),
                        xytext=(0, 5), textcoords='offset points',
                        ha='center', fontsize=7.5, color='0.3')

    ax.set_ylabel(r"$P_e \times$ abundance", fontsize=9)
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper left', fontsize=9)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.grid(True, which='major', alpha=0.25)

axes[-1].set_xlabel("Laser frequency offset from $^{174}$Yb  (MHz)", fontsize=10)

plt.tight_layout()
out = DIR / "Yb_spectrum_polarizations.png"
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved → {out}")
plt.show()
