"""inspect_db.py — verify NIST energy levels are stored in cm^-1."""
import sqlite3
from pathlib import Path

DB = Path(__file__).parent.parent.parent.parent / "atomicrabbit" / "db" / "spectroscopy.db"
conn = sqlite3.connect(str(DB))
conn.row_factory = sqlite3.Row

print("=== NIST energy levels (ordered by energy) ===")
for r in conn.execute(
    "SELECT level_id, configuration, term_symbol, J, parity, energy_cm1 "
    "FROM energy_levels WHERE source_type='nist' ORDER BY energy_cm1"
):
    print(f"  [{r['level_id']}] {r['configuration']!s:20s} {r['term_symbol']!s:6s} "
          f"J={r['J']} {r['parity']:4s}  {r['energy_cm1']:.3f} cm-1")

print()
print("=== NIST transitions ===")
for r in conn.execute(
    "SELECT t.transition_id, "
    "el.configuration AS lower_config, el.term_symbol AS lower_term, el.J AS lower_J, "
    "eu.configuration AS upper_config, eu.term_symbol AS upper_term, eu.J AS upper_J, "
    "t.wavelength_nm, t.einstein_A_s1 "
    "FROM transitions t "
    "JOIN energy_levels el ON el.level_id = t.lower_level_id "
    "JOIN energy_levels eu ON eu.level_id = t.upper_level_id "
    "WHERE t.source_type = 'nist' "
    "ORDER BY t.wavelength_nm"
):
    print(f"  [{r['transition_id']}] {r['wavelength_nm']:.4f} nm  "
          f"{r['lower_config']} {r['lower_term']} J={r['lower_J']} -> "
          f"{r['upper_config']} {r['upper_term']} J={r['upper_J']}  "
          f"A={r['einstein_A_s1']}")

conn.close()
