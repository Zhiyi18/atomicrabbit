"""
parser.py
Read a paper JSON file from db/papers/json/ and insert its records into the DB.

Resolves symbolic FK references (element/mass_number → atom_id, doi → paper_id,
config/term/J → level_id) before inserting. Warnings are issued for incomplete
data; errors per record are collected and returned rather than aborting the file.

Typical usage:
    from atomicrabbit.db.db_utils import get_connection, create_database
    from atomicrabbit.db.papers.parser import parse_file

    create_database('spectroscopy.db')
    with get_connection('spectroscopy.db') as conn:
        counts, errors = parse_file(conn, 'db/papers/json/steck_rb87_2001.json')
"""

import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Insertion order: records are sorted by this before processing
_TABLE_ORDER = {'atoms': 0, 'papers': 1, 'energy_levels': 2, 'transitions': 3}

# Symbolic fields that must be popped from data before building the INSERT;
# they are resolved to FK columns instead.
_SYMBOLIC_FIELDS = {
    'atoms':         set(),
    'papers':        {'data_types'},
    'energy_levels': {'element', 'mass_number', 'source_doi', 'source_title'},
    'transitions':   {
        'element', 'mass_number', 'source_doi', 'source_title',
        'upper_config', 'upper_term', 'upper_J',
        'lower_config', 'lower_term', 'lower_J',
    },
}


# ---------------------------------------------------------------------------
# FK resolution helpers
# ---------------------------------------------------------------------------

def _resolve_atom_id(conn, element, mass_number):
    row = conn.execute(
        "SELECT atom_id FROM atoms WHERE element = ? AND mass_number IS ?",
        (element, mass_number),
    ).fetchone()
    if row is None:
        raise ValueError(
            f"Atom not found: element={element!r}, mass_number={mass_number}. "
            "Insert the atoms record before energy_levels or transitions."
        )
    return row[0]


def _resolve_paper_id(conn, source_doi, source_title):
    """Look up paper_id by DOI; fall back to title for pre-DOI papers."""
    if source_doi is not None:
        row = conn.execute(
            "SELECT paper_id FROM papers WHERE doi = ?", (source_doi,)
        ).fetchone()
        if row:
            return row[0]
    if source_title is not None:
        row = conn.execute(
            "SELECT paper_id FROM papers WHERE title = ?", (source_title,)
        ).fetchone()
        if row:
            return row[0]
    raise ValueError(
        f"Paper not found: doi={source_doi!r}, title={source_title!r}. "
        "Insert the papers record before energy_levels or transitions."
    )


def _resolve_level_id(conn, atom_id, config, term, J, label):
    # ORDER BY is_recommended DESC so the preferred row wins when multiple
    # papers inserted the same (config, term, J) level.
    row = conn.execute(
        """SELECT level_id FROM energy_levels
           WHERE atom_id = ? AND configuration = ? AND term_symbol = ? AND J = ?
           ORDER BY is_recommended DESC""",
        (atom_id, config, term, float(J)),
    ).fetchone()
    if row is None:
        raise ValueError(
            f"{label} level not found: atom_id={atom_id}, "
            f"configuration={config!r}, term_symbol={term!r}, J={J}."
        )
    return row[0]


# ---------------------------------------------------------------------------
# Per-table insert handlers
# ---------------------------------------------------------------------------

def _insert_atom(conn, data):
    for field in ('mass_number', 'nuclear_spin'):
        if data.get(field) is None:
            logger.warning(
                "Atom '%s': %s is null — FK resolution for dependent records "
                "may be ambiguous.", data.get('element'), field
            )

    data = data.copy()
    _serialize_properties(data)
    row = {k: data.get(k) for k in ('element', 'mass_number', 'nuclear_spin', 'mu_I', 'abundance', 'properties')}

    # ON CONFLICT(element, mass_number) won't fire when mass_number is NULL
    # (SQLite treats NULL != NULL in UNIQUE constraints), so check manually.
    if data.get('mass_number') is None:
        existing = conn.execute(
            "SELECT atom_id FROM atoms WHERE element = ? AND mass_number IS NULL",
            (data['element'],),
        ).fetchone()
        if not existing:
            conn.execute(
                """INSERT INTO atoms (element, mass_number, nuclear_spin, mu_I, abundance, properties)
                   VALUES (:element, :mass_number, :nuclear_spin, :mu_I, :abundance, :properties)""",
                row,
            )
    else:
        conn.execute(
            """INSERT INTO atoms (element, mass_number, nuclear_spin, mu_I, abundance, properties)
               VALUES (:element, :mass_number, :nuclear_spin, :mu_I, :abundance, :properties)
               ON CONFLICT(element, mass_number) DO UPDATE SET
                   nuclear_spin = COALESCE(excluded.nuclear_spin, nuclear_spin),
                   mu_I         = COALESCE(excluded.mu_I, mu_I),
                   abundance    = COALESCE(excluded.abundance, abundance),
                   properties   = COALESCE(excluded.properties, properties)""",
            row,
        )


def _insert_paper(conn, data):
    data = data.copy()
    data_types = data.pop('data_types', [])
    row = {k: data.get(k) for k in ('title', 'authors', 'year', 'doi', 'journal', 'notes')}

    if data.get('doi') is None:
        # NULL doi: ON CONFLICT(doi) won't fire — deduplicate by title.
        existing = conn.execute(
            "SELECT paper_id FROM papers WHERE title = ?", (data['title'],)
        ).fetchone()
        if existing:
            paper_id = existing[0]
        else:
            cur = conn.execute(
                """INSERT INTO papers (title, authors, year, doi, journal, notes)
                   VALUES (:title, :authors, :year, :doi, :journal, :notes)""",
                row,
            )
            paper_id = cur.lastrowid
    else:
        cur = conn.execute(
            """INSERT INTO papers (title, authors, year, doi, journal, notes)
               VALUES (:title, :authors, :year, :doi, :journal, :notes)
               ON CONFLICT(doi) DO UPDATE SET
                   title   = excluded.title,
                   authors = excluded.authors,
                   year    = excluded.year,
                   journal = excluded.journal,
                   notes   = excluded.notes
               RETURNING paper_id""",
            row,
        )
        paper_id = cur.fetchone()[0]

    for dt in data_types:
        conn.execute(
            "INSERT OR IGNORE INTO paper_data_types (paper_id, data_type) VALUES (?, ?)",
            (paper_id, dt),
        )


def _insert_energy_level(conn, data):
    data = data.copy()
    _serialize_properties(data)

    element      = data.pop('element')
    mass_number  = data.pop('mass_number')
    source_doi   = data.pop('source_doi', None)
    source_title = data.pop('source_title', None)

    data['atom_id']          = _resolve_atom_id(conn, element, mass_number)
    data['source_paper_id']  = _resolve_paper_id(conn, source_doi, source_title)

    _generic_insert(conn, 'energy_levels', data)


def _insert_transition(conn, data):
    data = data.copy()
    _serialize_properties(data)

    element      = data.pop('element')
    mass_number  = data.pop('mass_number')
    source_doi   = data.pop('source_doi', None)
    source_title = data.pop('source_title', None)
    upper_config = data.pop('upper_config')
    upper_term   = data.pop('upper_term')
    upper_J      = data.pop('upper_J')
    lower_config = data.pop('lower_config')
    lower_term   = data.pop('lower_term')
    lower_J      = data.pop('lower_J')

    atom_id = _resolve_atom_id(conn, element, mass_number)

    data['atom_id']          = atom_id
    data['source_paper_id']  = _resolve_paper_id(conn, source_doi, source_title)
    data['upper_level_id']   = _resolve_level_id(conn, atom_id, upper_config, upper_term, upper_J, 'Upper')
    data['lower_level_id']   = _resolve_level_id(conn, atom_id, lower_config, lower_term, lower_J, 'Lower')

    _generic_insert(conn, 'transitions', data)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _serialize_properties(data):
    """Serialize properties dict to a JSON string in-place, if needed."""
    if isinstance(data.get('properties'), dict):
        data['properties'] = json.dumps(data['properties'])


def _generic_insert(conn, table, data):
    """Build and execute an INSERT from a flat dict of column → value."""
    cols         = ', '.join(data.keys())
    placeholders = ', '.join(f':{k}' for k in data.keys())
    conn.execute(f"INSERT INTO {table} ({cols}) VALUES ({placeholders})", data)


_HANDLERS = {
    'atoms':         _insert_atom,
    'papers':        _insert_paper,
    'energy_levels': _insert_energy_level,
    'transitions':   _insert_transition,
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_file(conn, json_path):
    """
    Read a paper JSON file and insert all records into the DB.

    Records are sorted by insertion order (atoms → papers → energy_levels →
    transitions) regardless of their order in the file. Errors on individual
    records are collected and returned rather than aborting the whole file.

    Note: the connection is NOT committed here — the caller controls that via
    the get_connection context manager. To roll back on any error, check the
    returned errors list before exiting the context.

    Parameters
    ----------
    conn      : sqlite3.Connection from db_utils.get_connection
    json_path : str or Path

    Returns
    -------
    counts : dict  — {'atoms': n, 'papers': n, 'energy_levels': n, 'transitions': n}
    errors : list  — error strings for any records that failed
    """
    json_path = Path(json_path)
    with open(json_path, encoding='utf-8') as f:
        records = json.load(f)

    records.sort(key=lambda r: _TABLE_ORDER.get(r['table'], 99))

    counts = {t: 0 for t in _TABLE_ORDER}
    errors = []

    for record in records:
        table   = record.get('table')
        handler = _HANDLERS.get(table)

        if handler is None:
            logger.warning("[%s] Unknown table %r — skipped.", json_path.name, table)
            continue

        try:
            handler(conn, record['data'])
            counts[table] += 1
        except Exception as e:
            msg = f"[{json_path.name}] {table}: {e}"
            logger.error(msg)
            errors.append(msg)

    logger.info(
        "[%s] inserted %s — %d error(s).",
        json_path.name, counts, len(errors),
    )
    return counts, errors


def parse_directory(conn, json_dir):
    """
    Run parse_file on every .json file in json_dir.

    Returns
    -------
    results : dict mapping filename → (counts, errors)
    """
    json_dir = Path(json_dir)
    results  = {}

    for path in sorted(json_dir.glob('*.json')):
        logger.info("Parsing %s ...", path.name)
        results[path.name] = parse_file(conn, path)

    total_errors = sum(len(e) for _, e in results.values())
    logger.info("parse_directory done — %d file(s), %d total error(s).", len(results), total_errors)
    return results
