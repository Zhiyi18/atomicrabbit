"""
example_usage.py
Demonstrates: create db -> write Rb-87 D2-line data -> read it back.
Covers both the NIST path (direct insert) and the paper path (manual call).
"""

from atomicrabbit.db.db_utils import (
    get_connection, create_database,
    add_atom, add_paper, add_energy_level, add_transition,
    get_transitions_for_atom, get_energy_level, get_paper_for_source,
)
from atomicrabbit.db.nist.load_tables import load_nist_to_db

DB_PATH = "spectroscopy.db"

# ---------------------------------------------------------------------
# 1. CREATE
# ---------------------------------------------------------------------
create_database(DB_PATH)

# ---------------------------------------------------------------------
# 2a. NIST path — direct insert, no JSON intermediary
# ---------------------------------------------------------------------
with get_connection(DB_PATH) as conn:
    load_nist_to_db(conn, element="Rb", mass_number=87,
                    wavelength_min_nm=750, wavelength_max_nm=800,
                    nuclear_spin=1.5)

# ---------------------------------------------------------------------
# 2b. Paper path — manual call mirroring what the paper parser will do
# ---------------------------------------------------------------------
with get_connection(DB_PATH) as conn:
    rb87 = add_atom(conn, element="Rb", mass_number=87,
                    nuclear_spin=1.5, mu_I=-1.8342)

    ground = add_energy_level(
        conn, atom_id=rb87, energy_cm1=0.0,
        configuration="5s", term_symbol="2S", L=0, S=0.5, J=0.5,
        source_type="nist", is_recommended=1,
    )
    excited = add_energy_level(
        conn, atom_id=rb87, energy_cm1=12816.545,
        configuration="5p", term_symbol="2P", L=1, S=0.5, J=1.5,
        hyperfine_A_MHz=84.7, hyperfine_B_MHz=12.5,
        source_type="nist", is_recommended=1,
    )

    paper_id = add_paper(
        conn,
        title="Precision isotope shift measurements in Rb D2 line",
        authors="Doe, J. et al.",
        year=2019,
        doi="10.1234/example.rb87",
        data_types=["isotope_shift"],
    )

    add_transition(
        conn, atom_id=rb87, upper_level_id=excited, lower_level_id=ground,
        wavelength_nm=780.24, frequency_THz=384.23,
        einstein_A_s1=3.8e7,
        isotope_shift_MHz=78.3, isotope_shift_ref_mass=85,
        source_type="paper", source_paper_id=paper_id, is_recommended=0,
    )

# ---------------------------------------------------------------------
# 3. READ — simulation lookup flow
# ---------------------------------------------------------------------
with get_connection(DB_PATH) as conn:
    transitions = get_transitions_for_atom(conn, element="Rb", mass_number=87)

    papers_used = set()
    for t in transitions:
        upper = get_energy_level(conn, t["upper_level_id"])
        lower = get_energy_level(conn, t["lower_level_id"])
        paper_title = get_paper_for_source(conn, t["source_type"], t["source_paper_id"])
        if paper_title:
            papers_used.add(paper_title)
        print(f"{lower['term_symbol']} → {upper['term_symbol']} "
              f"@ {t['wavelength_nm']} nm  [recommended={t['is_recommended']}]")

    print("\nData sourced from papers:")
    for p in papers_used:
        print(" ", p)
