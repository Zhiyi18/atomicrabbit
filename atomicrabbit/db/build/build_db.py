"""
build_db.py
Build spectroscopy.db from all element directories under db/json/.

Each element directory that contains _manifest.json is processed:
  1. NIST ASD live query for each specified mass number and wavelength range
  2. Paper JSON files parsed from <element>/paper/
  3. NIST snapshot JSONs parsed from <element>/nist/ (when present)

Usage (from project root):
    python -m atomicrabbit.db.build.build_db
"""
import sys
import importlib.util
import logging
from pathlib import Path

# Locate project root: build/ -> db/ -> atomicrabbit/ -> project root
ROOT     = Path(__file__).resolve().parent.parent.parent.parent
DB_PKG   = ROOT / "atomicrabbit" / "db"
DB_PATH  = DB_PKG / "spectroscopy.db"
JSON_DIR = DB_PKG / "json"


# ---------------------------------------------------------------------------
# Direct module loading — avoids triggering physics/__init__ which needs numpy
# ---------------------------------------------------------------------------

def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_db_utils     = _load("atomicrabbit.db.db_utils",               DB_PKG / "db_utils.py")
_nist_mod     = _load("atomicrabbit.db.nist_tools.load_tables", DB_PKG / "nist_tools" / "load_tables.py")
_paper_parser = _load("atomicrabbit.db.json.paper_parser",      DB_PKG / "json" / "paper_parser.py")

get_connection  = _db_utils.get_connection
create_database = _db_utils.create_database
load_nist_to_db = _nist_mod.load_nist_to_db
parse_file      = _paper_parser.parse_file

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")


# ---------------------------------------------------------------------------
# Per-element loader
# ---------------------------------------------------------------------------

def load_element(conn, element_dir):
    """
    Load one element directory into the DB.

    Reads _manifest.json, queries NIST, then parses paper and nist-snapshot
    JSON files. Returns (n_paper_errors, n_snap_errors).
    """
    import json as _json  # built-in, not local package

    manifest_path = element_dir / "_manifest.json"
    with open(manifest_path, encoding='utf-8') as f:
        manifest = _json.load(f)

    linename    = manifest["linename"]
    element     = linename.split()[0]
    nist_masses = manifest.get("nist_mass_numbers", [])
    nist_ranges = manifest.get("nist_wavelength_ranges_nm", [[100, 500], [500, 1000], [1000, 2000]])
    paper_order = manifest.get("paper_order", [])

    print(f"\n{'='*60}")
    print(f"  {linename}")
    print(f"{'='*60}")

    # -- NIST live queries --------------------------------------------------
    if nist_masses:
        n_queries = len(nist_masses) * len(nist_ranges)
        print(f"  NIST: {len(nist_masses)} mass(es) × {len(nist_ranges)} range(s) = {n_queries} queries")
        for mass in nist_masses:
            for lo, hi in nist_ranges:
                try:
                    load_nist_to_db(conn, element=element, mass_number=mass,
                                    wavelength_min_nm=lo, wavelength_max_nm=hi,
                                    linename=linename)
                except Exception as e:
                    print(f"  WARNING: NIST {linename} mass={mass} [{lo}-{hi} nm] failed: {e}")

    # -- Paper JSON files ---------------------------------------------------
    paper_errors = 0
    paper_dir = element_dir / "papers"
    if paper_dir.exists():
        if paper_order:
            paths = [paper_dir / f for f in paper_order]
            for p in paths:
                if not p.exists():
                    print(f"  MISSING  {p.name}")
        else:
            paths = sorted(paper_dir.glob("*.json"))

        print(f"  Papers:  {len(paths)} file(s)")
        for path in paths:
            counts, errors = parse_file(conn, path)
            status = "OK  " if not errors else "FAIL"
            print(f"    [{status}] {path.name}  {dict(counts)}")
            for e in errors:
                print(f"           ERROR: {e.split('] ', 1)[-1]}")
            paper_errors += len(errors)

    # -- NIST snapshot JSONs (future) ---------------------------------------
    snap_errors = 0
    nist_snap_dir = element_dir / "nist"
    if nist_snap_dir.exists():
        snap_paths = sorted(nist_snap_dir.glob("*.json"))
        if snap_paths:
            print(f"  NIST snapshots: {len(snap_paths)} file(s)")
            for path in snap_paths:
                counts, errors = parse_file(conn, path)
                status = "OK  " if not errors else "FAIL"
                print(f"    [{status}] {path.name}  {dict(counts)}")
                for e in errors:
                    print(f"           ERROR: {e.split('] ', 1)[-1]}")
                snap_errors += len(errors)

    return paper_errors, snap_errors


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    print(f"Database : {DB_PATH}")
    print(f"JSON dir : {JSON_DIR}\n")
    create_database(str(DB_PATH))

    element_dirs = sorted(
        p for p in JSON_DIR.iterdir()
        if p.is_dir()
        and not p.name.startswith(('_', '.'))
        and (p / "_manifest.json").exists()
    )

    if not element_dirs:
        print("No element directories with _manifest.json found.")
        return 0

    print(f"Elements found: {[d.name for d in element_dirs]}")

    total_paper_errors = 0
    total_snap_errors  = 0
    with get_connection(str(DB_PATH)) as conn:
        for element_dir in element_dirs:
            pe, se = load_element(conn, element_dir)
            total_paper_errors += pe
            total_snap_errors  += se

    # Summary
    print(f"\n{'='*60}")
    with get_connection(str(DB_PATH)) as conn:
        def _count(q):
            return conn.execute(q).fetchone()[0]
        n_atoms   = _count("SELECT COUNT(*) FROM atoms")
        n_papers  = _count("SELECT COUNT(*) FROM papers")
        n_levels  = _count("SELECT COUNT(*) FROM energy_levels")
        n_trans   = _count("SELECT COUNT(*) FROM transitions")
        n_nist_l  = _count("SELECT COUNT(*) FROM energy_levels WHERE source_type='nist'")
        n_paper_l = _count("SELECT COUNT(*) FROM energy_levels WHERE source_type='paper'")
        n_nist_t  = _count("SELECT COUNT(*) FROM transitions WHERE source_type='nist'")
        n_paper_t = _count("SELECT COUNT(*) FROM transitions WHERE source_type='paper'")

    print(f"  DB summary:")
    print(f"    atoms          : {n_atoms}")
    print(f"    papers         : {n_papers}")
    print(f"    energy_levels  : {n_levels}  (nist={n_nist_l}, paper={n_paper_l})")
    print(f"    transitions    : {n_trans}  (nist={n_nist_t}, paper={n_paper_t})")
    print(f"  Paper errors     : {total_paper_errors}")
    print(f"  Snapshot errors  : {total_snap_errors}")
    print(f"{'='*60}")

    return total_paper_errors + total_snap_errors


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
