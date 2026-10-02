"""
db_utils.py
Thin helper layer around the spectroscopy SQLite database.

Usage pattern this supports (matches your simulation workflow):
    1. get_transitions_for_atom(...) to find candidate transitions
    2. get_energy_level(...) / get_papers_for_source(...) to pull details
       only when needed
    3. build your Hamiltonian / collapse operators from the results
"""

import sqlite3
from pathlib import Path
from contextlib import contextmanager

SCHEMA_PATH = Path(__file__).parent / "schema.sql"


# ---------------------------------------------------------------------
# Connection / setup
# ---------------------------------------------------------------------
@contextmanager
def get_connection(db_path="spectroscopy.db"):
    """Context-managed connection with foreign keys on and dict-like rows."""
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.row_factory = sqlite3.Row
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def create_database(db_path="spectroscopy.db", schema_path=SCHEMA_PATH):
    """Run schema.sql against db_path. Safe to call repeatedly (uses IF NOT EXISTS)."""
    with open(schema_path, "r") as f:
        schema_sql = f.read()
    with get_connection(db_path) as conn:
        conn.executescript(schema_sql)


# ---------------------------------------------------------------------
# Write helpers
# ---------------------------------------------------------------------
def add_atom(conn, element, mass_number=None, nuclear_spin=None, mu_I=None,
             abundance=None, properties=None):
    cur = conn.execute(
        """INSERT INTO atoms (element, mass_number, nuclear_spin, mu_I, abundance, properties)
           VALUES (?, ?, ?, ?, ?, ?)
           ON CONFLICT(element, mass_number) DO UPDATE SET
               nuclear_spin=excluded.nuclear_spin,
               mu_I=excluded.mu_I,
               abundance=excluded.abundance,
               properties=excluded.properties
           RETURNING atom_id""",
        (element, mass_number, nuclear_spin, mu_I, abundance, properties),
    )
    return cur.fetchone()[0]


def add_paper(conn, title, authors=None, year=None, doi=None, journal=None,
              notes=None, data_types=None, references_titles=None):
    """
    data_types: list like ['transition', 'isotope_shift']
    references_titles: list of paper titles this paper cites (must already exist in DB)
    """
    cur = conn.execute(
        """INSERT INTO papers (title, authors, year, doi, journal, notes)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (title, authors, year, doi, journal, notes),
    )
    paper_id = cur.lastrowid

    for dt in data_types or []:
        conn.execute(
            "INSERT OR IGNORE INTO paper_data_types (paper_id, data_type) VALUES (?, ?)",
            (paper_id, dt),
        )

    for ref_title in references_titles or []:
        row = conn.execute(
            "SELECT paper_id FROM papers WHERE title = ?", (ref_title,)
        ).fetchone()
        if row is None:
            raise ValueError(f"Referenced paper not found yet: {ref_title!r}")
        conn.execute(
            """INSERT OR IGNORE INTO paper_references (citing_paper_id, cited_paper_id)
               VALUES (?, ?)""",
            (paper_id, row["paper_id"]),
        )

    return paper_id


def add_energy_level(conn, atom_id, energy_cm1, configuration=None, term_symbol=None,
                      L=None, S=None, J=None, F=None, parity=None,
                      energy_uncertainty=None, hyperfine_A_MHz=None, hyperfine_B_MHz=None,
                      lifetime_s=None, is_recommended=0, properties=None,
                      source_type="nist", source_paper_id=None):
    cur = conn.execute(
        """INSERT INTO energy_levels
           (atom_id, configuration, term_symbol, L, S, J, F, parity, energy_cm1,
            energy_uncertainty, hyperfine_A_MHz, hyperfine_B_MHz, lifetime_s,
            is_recommended, properties, source_type, source_paper_id)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (atom_id, configuration, term_symbol, L, S, J, F, parity, energy_cm1,
         energy_uncertainty, hyperfine_A_MHz, hyperfine_B_MHz, lifetime_s,
         is_recommended, properties, source_type, source_paper_id),
    )
    return cur.lastrowid


def add_transition(conn, atom_id, upper_level_id, lower_level_id,
                    wavelength_nm=None, frequency_THz=None, einstein_A_s1=None,
                    linewidth_MHz=None, isotope_shift_MHz=None, isotope_shift_ref_mass=None,
                    branching_ratio=None, is_recommended=0, properties=None,
                    source_type="nist", source_paper_id=None):
    cur = conn.execute(
        """INSERT INTO transitions
           (atom_id, upper_level_id, lower_level_id, wavelength_nm, frequency_THz,
            einstein_A_s1, linewidth_MHz, isotope_shift_MHz, isotope_shift_ref_mass,
            branching_ratio, is_recommended, properties, source_type, source_paper_id)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (atom_id, upper_level_id, lower_level_id, wavelength_nm, frequency_THz,
         einstein_A_s1, linewidth_MHz, isotope_shift_MHz, isotope_shift_ref_mass,
         branching_ratio, is_recommended, properties, source_type, source_paper_id),
    )
    return cur.lastrowid


# ---------------------------------------------------------------------
# Read helpers
# ---------------------------------------------------------------------
def get_transitions_for_atom(conn, element, mass_number=None,
                              wavelength_range_nm=None):
    """Find candidate transitions for a given isotope, optionally filtered by
    wavelength range (nm) -- e.g. to find transitions near your laser."""
    query = """
        SELECT t.* FROM transitions t
        JOIN atoms a ON a.atom_id = t.atom_id
        WHERE a.element = ?
    """
    params = [element]
    if mass_number is not None:
        query += " AND a.mass_number = ?"
        params.append(mass_number)
    if wavelength_range_nm is not None:
        lo, hi = wavelength_range_nm
        query += " AND t.wavelength_nm BETWEEN ? AND ?"
        params += [lo, hi]
    return conn.execute(query, params).fetchall()


def get_energy_level(conn, level_id):
    return conn.execute(
        "SELECT * FROM energy_levels WHERE level_id = ?", (level_id,)
    ).fetchone()


def get_paper(conn, paper_id):
    return conn.execute(
        "SELECT * FROM papers WHERE paper_id = ?", (paper_id,)
    ).fetchone()


def get_paper_for_source(conn, source_type, source_paper_id):
    """Given a transition/energy_level row's source fields, return the paper
    title if it came from a paper, else None."""
    if source_type != "paper" or source_paper_id is None:
        return None
    row = get_paper(conn, source_paper_id)
    return row["title"] if row else None


def get_referenced_papers(conn, paper_id):
    """All papers that `paper_id` cites."""
    return conn.execute(
        """SELECT p.* FROM papers p
           JOIN paper_references r ON r.cited_paper_id = p.paper_id
           WHERE r.citing_paper_id = ?""",
        (paper_id,),
    ).fetchall()
