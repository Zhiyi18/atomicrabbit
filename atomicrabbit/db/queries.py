"""
queries.py — read-only SQL helpers for constructing physics objects from spectroscopy.db.
No physics imports; returns plain dicts only.
"""

# Fixed citation for all NIST ASD data
_NIST_ASD_BIBTEX = """\
@misc{NIST_ASD,
  author       = {Kramida, A. and Ralchenko, Yu. and Reader, J. and {NIST ASD Team}},
  title        = {{NIST} Atomic Spectra Database (ver. 5.11)},
  year         = {2024},
  howpublished = {\\url{https://physics.nist.gov/asd}},
  note         = {National Institute of Standards and Technology, Gaithersburg, MD. DOI: 10.18434/T4W30F}
}"""

_NIST_ASD_TEXT = (
    "Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024). "
    "NIST Atomic Spectra Database (ver. 5.11). "
    "National Institute of Standards and Technology, Gaithersburg, MD. "
    "https://doi.org/10.18434/T4W30F"
)


def get_atom_properties(conn, element, mass_number):
    """Return {nuclear_spin, mu_I} for a specific isotope, or {} if not found."""
    row = conn.execute(
        "SELECT nuclear_spin, mu_I FROM atoms WHERE element = ? AND mass_number = ?",
        (element, mass_number),
    ).fetchone()
    return dict(row) if row else {}


def get_fine_states(conn, element, mass_number):
    """
    Return is_recommended energy levels for a specific isotope, ordered by energy.

    Each dict has: configuration, term_symbol, L, S, J, parity, energy_cm1,
                   hyperfine_A_MHz, hyperfine_B_MHz.
    """
    rows = conn.execute(
        """
        SELECT el.configuration, el.term_symbol, el.L, el.S, el.J,
               el.parity, el.energy_cm1, el.hyperfine_A_MHz, el.hyperfine_B_MHz
        FROM energy_levels el
        JOIN atoms a ON el.atom_id = a.atom_id
        WHERE a.element = ? AND a.mass_number = ? AND el.is_recommended = 1
        ORDER BY el.energy_cm1
        """,
        (element, mass_number),
    ).fetchall()
    return [dict(r) for r in rows]


def get_nist_transitions(conn, element):
    """
    Return NIST transitions with their level energies for cross-referencing.

    Each dict has: upper_E (cm⁻¹), upper_J, lower_E (cm⁻¹), lower_J,
                   einstein_A_s1, wavelength_nm.
    Only transitions with a non-null einstein_A_s1 are returned.

    Note: NIST configurations use dot notation with closed shells (e.g. '4f14.6s2')
    while paper configurations use shorthand ('6s2'). Config strings therefore
    cannot be used as join keys — the caller must match by energy proximity.
    """
    rows = conn.execute(
        """
        SELECT eu.energy_cm1 AS upper_E, eu.J AS upper_J,
               el.energy_cm1 AS lower_E, el.J AS lower_J,
               t.einstein_A_s1, t.wavelength_nm
        FROM transitions t
        JOIN atoms a         ON t.atom_id          = a.atom_id
        JOIN energy_levels eu ON t.upper_level_id  = eu.level_id
        JOIN energy_levels el ON t.lower_level_id  = el.level_id
        WHERE a.element = ? AND a.mass_number IS NULL
          AND t.source_type = 'nist' AND t.is_recommended = 1
          AND t.einstein_A_s1 IS NOT NULL
        """,
        (element,),
    ).fetchall()
    return [dict(r) for r in rows]


def get_sources(conn, element, mass_number):
    """
    Return the distinct sources used for is_recommended data for (element, mass_number).

    Covers NIST transitions (element-wide) and paper-sourced energy levels /
    hyperfine constants for the specific isotope.

    Returns a list of dicts.  NIST entry: {'source_type': 'nist'}.
    Paper entry: {'source_type': 'paper', 'paper_id', 'title', 'authors',
                  'year', 'doi', 'journal'}.
    """
    # NIST transitions (mass_number IS NULL in DB, shared across all isotopes)
    nist_exists = conn.execute(
        """
        SELECT 1 FROM transitions t
        JOIN atoms a ON t.atom_id = a.atom_id
        WHERE a.element = ? AND a.mass_number IS NULL
          AND t.source_type = 'nist' AND t.is_recommended = 1
        LIMIT 1
        """,
        (element,),
    ).fetchone()

    # Paper-sourced energy levels for this isotope
    rows = conn.execute(
        """
        SELECT DISTINCT p.paper_id, p.title, p.authors, p.year, p.doi, p.journal
        FROM energy_levels el
        JOIN papers p ON el.source_paper_id = p.paper_id
        JOIN atoms a  ON el.atom_id = a.atom_id
        WHERE a.element = ? AND a.mass_number = ?
          AND el.is_recommended = 1 AND el.source_type = 'paper'
        ORDER BY p.year, p.paper_id
        """,
        (element, mass_number),
    ).fetchall()

    sources = []
    if nist_exists:
        sources.append({'source_type': 'nist'})
    for r in rows:
        sources.append({'source_type': 'paper', **dict(r)})
    return sources


def _bibtex_key(source):
    first_author_last = source['authors'].split(';')[0].split(',')[0].strip()
    return f"{first_author_last}{source['year']}"


def format_sources(sources, fmt='bibtex'):
    """
    Format a list of sources returned by get_sources() as BibTeX or plain text.

    Parameters
    ----------
    sources : list[dict]   Output of get_sources().
    fmt     : str          'bibtex' (default) or 'text'.

    Returns
    -------
    str
    """
    if fmt == 'bibtex':
        entries = []
        for s in sources:
            if s['source_type'] == 'nist':
                entries.append(_NIST_ASD_BIBTEX)
            else:
                key = _bibtex_key(s)
                authors_bibtex = ' and '.join(a.strip() for a in s['authors'].split(';'))
                lines = [
                    f"@article{{{key},",
                    f"  author  = {{{authors_bibtex}}},",
                    f"  title   = {{{{{s['title']}}}}},",
                    f"  journal = {{{s['journal']}}},",
                    f"  year    = {{{s['year']}}},",
                ]
                if s['doi']:
                    lines.append(f"  doi     = {{{s['doi']}}},")
                lines.append("}")
                entries.append('\n'.join(lines))
        return '\n\n'.join(entries)

    elif fmt == 'text':
        lines = []
        for i, s in enumerate(sources, 1):
            if s['source_type'] == 'nist':
                lines.append(f"[{i}] {_NIST_ASD_TEXT}")
            else:
                ref = f"[{i}] {s['authors']} ({s['year']}). {s['title']}. {s['journal']}."
                if s['doi']:
                    ref += f" DOI: {s['doi']}"
                lines.append(ref)
        return '\n'.join(lines)

    else:
        raise ValueError(f"Unknown format '{fmt}'. Use 'bibtex' or 'text'.")
