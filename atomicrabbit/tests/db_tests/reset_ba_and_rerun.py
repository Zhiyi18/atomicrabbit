"""Reset all Ba data from the DB and re-run the parser. Reports clean error counts."""
import logging
import sqlite3
from pathlib import Path
from atomicrabbit.db.db_utils import get_connection
from atomicrabbit.db.papers.parser import parse_directory

logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')

DB_PATH = "atomicrabbit/db/spectroscopy.db"
JSON_DIR = "atomicrabbit/db/papers/json_barium_ion"

# --- reset: delete all Ba rows so re-run is idempotent ---
with get_connection(DB_PATH) as conn:
    atom_ids = [r[0] for r in conn.execute(
        "SELECT atom_id FROM atoms WHERE element='Ba'"
    ).fetchall()]
    if atom_ids:
        placeholders = ",".join("?" * len(atom_ids))
        conn.execute(f"DELETE FROM transitions WHERE atom_id IN ({placeholders})", atom_ids)
        conn.execute(f"DELETE FROM energy_levels WHERE atom_id IN ({placeholders})", atom_ids)
        conn.execute(f"DELETE FROM atoms WHERE atom_id IN ({placeholders})", atom_ids)
        # delete papers that only exist for Ba (by DOI list from json files)
        ba_dois = [
            "10.1007/BF01413460",
            "10.1103/PhysRevA.33.2117",
            "10.1103/PhysRevLett.119.100501",
            "10.1103/PhysRevA.25.1476",
            "10.1103/PhysRevA.88.012518",
            "10.1103/PhysRevA.86.022512",
            "10.1016/j.adt.2005.04.001",
            "10.1007/s10751-019-1619-x",
            "10.1103/PhysRevA.101.042509",
        ]
        conn.execute(
            "DELETE FROM papers WHERE doi IN ({})".format(",".join("?" * len(ba_dois))),
            ba_dois,
        )
        # null-doi papers by title
        ba_titles = [
            "Hyperfine structure and isotope shifts of Ba+ by fast-beam laser spectroscopy",
            "Hyperfine structure and nuclear magnetic moments of the stable barium isotopes",
        ]
        conn.execute(
            "DELETE FROM papers WHERE title IN ({})".format(",".join("?" * len(ba_titles))),
            ba_titles,
        )
    print("DB reset done.")

# --- re-run parser ---
with get_connection(DB_PATH) as conn:
    results = parse_directory(conn, JSON_DIR)

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
total_errors = 0
for fname, (counts, errors) in sorted(results.items()):
    status = "OK" if not errors else f"{len(errors)} error(s)"
    print(f"\n[{status}] {fname}")
    print(f"  inserted: {counts}")
    for e in errors:
        print(f"  ERROR: {e.split('] ', 1)[-1]}")
    total_errors += len(errors)

print(f"\nTotal errors: {total_errors}")
