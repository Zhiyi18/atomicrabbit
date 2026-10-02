"""
load_tables.py
Query NIST ASD via astroquery and write results directly to the spectroscopy DB.
No JSON intermediary — NIST data is machine-readable and written inline.
"""

import re
from fractions import Fraction
from astroquery.nist import Nist
import astropy.units as u

L_LETTERS = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5}


def parse_term_symbol(term_str):
    """'2P*_1' or '3D_2' → {'S', 'L', 'J', 'parity', 'term_symbol'}, or None."""
    if not isinstance(term_str, str):
        return None
    parity = 'odd' if '*' in term_str else 'even'
    match = re.match(r'(\d+)([SPDFGH])\*?_([\d./]+)', term_str.strip())
    if not match:
        return None
    multiplicity = int(match.group(1))
    L_letter = match.group(2)
    try:
        J = float(Fraction(match.group(3)))
    except (ValueError, ZeroDivisionError):
        return None
    return {
        'S': (multiplicity - 1) / 2,
        'L': L_LETTERS.get(L_letter),
        'J': J,
        'parity': parity,
        'term_symbol': term_str.strip(),
    }


def parse_level(level_str, energy_cm):
    """
    Parse a NIST level string 'config | term | J' into a DB-ready dict, or None.
    Handles LS-coupling terms (e.g. '3P*') and jj-coupling terms (e.g. '(7/2,5/2)*')
    with L=None, S=None for the jj case.
    """
    if not isinstance(level_str, str):
        return None
    parts = [p.strip() for p in level_str.split('|')]
    if len(parts) != 3:
        return None
    config, term, J_str = parts
    try:
        J = float(Fraction(J_str.strip()))
    except (ValueError, ZeroDivisionError):
        return None
    parity = 'odd' if '*' in term else 'even'
    match = re.match(r'(\d+)([SPDFGH])\*?', term.strip())
    if match:
        L_letter = match.group(2)
        if L_letter not in L_LETTERS:
            return None
        multiplicity = int(match.group(1))
        L = L_LETTERS[L_letter]
        S = (multiplicity - 1) / 2
    else:
        # jj-coupling or non-standard term (e.g. '(7/2,5/2)*') — store without L/S
        if not term:
            return None
        L = None
        S = None
    return {
        'configuration': config,
        'term_symbol': term.strip().rstrip('*'),  # parity stored in parity column; * is redundant
        'S': S,
        'L': L,
        'J': J,
        'parity': parity,
        'energy_cm1': float(energy_cm),
    }


def _find_col(colnames, *keywords):
    """Return the first column name containing any keyword (case-insensitive)."""
    for name in colnames:
        for kw in keywords:
            if kw.lower() in name.lower():
                return name
    return None


def _get_or_insert_level(conn, atom_id, level, is_recommended):
    """Return level_id of an existing NIST level, or insert and return the new id."""
    row = conn.execute(
        """SELECT level_id FROM energy_levels
           WHERE atom_id = ? AND configuration = ? AND term_symbol = ? AND J = ?
             AND source_type = 'nist'""",
        (atom_id, level['configuration'], level['term_symbol'], level['J'])
    ).fetchone()
    if row:
        return row[0]
    cur = conn.execute(
        """INSERT INTO energy_levels
           (atom_id, configuration, term_symbol, L, S, J, parity,
            energy_cm1, source_type, is_recommended)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, 'nist', ?)""",
        (atom_id, level['configuration'], level['term_symbol'],
         level['L'], level['S'], level['J'], level['parity'],
         level['energy_cm1'], is_recommended)
    )
    return cur.lastrowid


def load_nist_to_db(conn, element, mass_number, wavelength_min_nm, wavelength_max_nm,
                    nuclear_spin=None, linename=None, is_recommended=1):
    """
    Query NIST ASD for an element over a wavelength range and write to the DB.

    Parameters
    ----------
    conn : sqlite3.Connection
        Open connection from db_utils.get_connection. Caller owns commit/rollback.
    element : str
        Element symbol, e.g. 'Rb'.
    mass_number : int or None
        Isotope mass number; None for natural / unspecified.
    wavelength_min_nm, wavelength_max_nm : float
        Query range in nanometres.
    nuclear_spin : float or None
        Stored on the atom row if provided.
    linename : str or None
        NIST line name, e.g. 'Rb I'. Defaults to '{element} I'.
    is_recommended : int
        Applied to all inserted levels and transitions. Defaults to 1 so that
        NIST data is immediately usable as the simulation baseline.

    Returns
    -------
    atom_id : int
    """
    if linename is None:
        linename = f'{element} I'

    table = Nist.query(
        wavelength_min_nm * u.nm,
        wavelength_max_nm * u.nm,
        linename=linename,
    )

    ei_col = _find_col(table.colnames, 'Ei')
    wl_col = _find_col(table.colnames, 'obs_wl', 'Observed', 'wl_vac', 'wl_air', 'Ritz')
    a_col  = _find_col(table.colnames, 'Aki', 'A(s')

    if ei_col is None:
        raise ValueError(
            f"Cannot find energy column in NIST table. Columns: {table.colnames}"
        )

    # insert or retrieve atom row
    existing = conn.execute(
        "SELECT atom_id FROM atoms WHERE element = ? AND mass_number IS ?",
        (element, mass_number)
    ).fetchone()
    if existing:
        atom_id = existing[0]
        if nuclear_spin is not None:
            conn.execute(
                "UPDATE atoms SET nuclear_spin = ? WHERE atom_id = ?",
                (nuclear_spin, atom_id)
            )
    else:
        cur = conn.execute(
            "INSERT INTO atoms (element, mass_number, nuclear_spin) VALUES (?, ?, ?)",
            (element, mass_number, nuclear_spin)
        )
        atom_id = cur.lastrowid

    # NIST ASD returns level energies in eV in the Ei-Ek column.
    EV_TO_CM1 = 8065.543937  # 1 eV in cm^-1 (NIST CODATA 2018)

    inserted = skipped = 0

    for row in table:
        # parse energies from the combined 'Ei - Ek' column (values are in eV)
        try:
            ei_str, ek_str = str(row[ei_col]).split('-')
            ei = float(ei_str.strip()) * EV_TO_CM1
            ek = float(ek_str.strip()) * EV_TO_CM1
        except (ValueError, AttributeError):
            skipped += 1
            continue

        lower = parse_level(str(row['Lower level']), ei)
        upper = parse_level(str(row['Upper level']), ek)
        if lower is None or upper is None:
            skipped += 1
            continue

        lower_id = _get_or_insert_level(conn, atom_id, lower, is_recommended)
        upper_id = _get_or_insert_level(conn, atom_id, upper, is_recommended)

        # skip duplicate NIST transitions (NULL != NULL in UNIQUE, so check explicitly)
        if conn.execute(
            """SELECT 1 FROM transitions
               WHERE upper_level_id = ? AND lower_level_id = ? AND source_type = 'nist'""",
            (upper_id, lower_id)
        ).fetchone():
            skipped += 1
            continue

        wl_nm = None
        if wl_col:
            try:
                wl_nm = float(row[wl_col])
            except (ValueError, TypeError):
                pass

        A = None
        if a_col:
            try:
                A = float(row[a_col])
            except (ValueError, TypeError):
                pass

        conn.execute(
            """INSERT INTO transitions
               (atom_id, upper_level_id, lower_level_id,
                wavelength_nm, einstein_A_s1, source_type, is_recommended)
               VALUES (?, ?, ?, ?, ?, 'nist', ?)""",
            (atom_id, upper_id, lower_id, wl_nm, A, is_recommended)
        )
        inserted += 1

    print(f"NIST {linename} [{wavelength_min_nm}-{wavelength_max_nm} nm]: "
          f"{inserted} transitions inserted, {skipped} skipped.")
    return atom_id
