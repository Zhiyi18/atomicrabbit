"""
test_parser.py  —  run all debug + Yb JSON files through parse_file and report results.
"""
import sys, tempfile, logging
from pathlib import Path

ROOT = Path(__file__).parent.parent.parent.parent  # yb_test -> db -> atomicrabbit(pkg) -> project root
DB_PKG = ROOT / "atomicrabbit" / "db"

# Import db modules directly to avoid triggering atomicrabbit/__init__.py
# (which pulls in numpy/astroquery that may not be installed in this env).
import importlib.util

def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

db_utils = _load("atomicrabbit.db.db_utils",    DB_PKG / "db_utils.py")
parser   = _load("atomicrabbit.db.papers.parser", DB_PKG / "papers" / "parser.py")

get_connection  = db_utils.get_connection
create_database = db_utils.create_database
parse_file      = parser.parse_file

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

JSON_DEBUG = ROOT / "atomicrabbit" / "db" / "papers" / "json_debug"
JSON_YB    = ROOT / "atomicrabbit" / "db" / "papers" / "json_ytterbium"

DEBUG_ORDER = [
    "steck_rb87_2001.json",       # creates Rb87 atom + levels
    "heavens_rb_1961.json",        # needs Rb87 levels (is_recommended=0 alternative)
    "hyperfine_update_rb87.json",  # needs Rb87 levels (second paper, same level type)
    "missing_isotope_ba.json",     # null mass_number; parser should warn but insert
    "novel_params_cs.json",        # properties dict on levels + transitions
    "no_doi_old_paper.json",       # null DOI; paper OK; transition EXPECTED FAIL (Rb/null atom)
    "partial_data_sr.json",        # nullable wavelength + A coefficient
]

YB_ORDER = [
    "nudat3_yb_nuclear_moments.json",      # atoms only, no deps
    "kleinert_yb_2016_399nm.json",         # 399nm, 1P1 for all isotopes
    "atkinson_yb_2019_556nm.json",         # 556nm 3P1 -- inserts 174Yb 3P1 (needed by suryanarayana)
    "kawasaki_yb_2024_431nm.json",         # 431nm inner-shell
    "jones_yb171_2023_556nm.json",         # 171Yb 3P1 high-precision
    "kim_yb171_2021_clock.json",           # 171Yb clock 578nm
    "qiao_yb_2023_770nm_mts.json",         # 174Yb 770nm 3P2->3S1
    "qiao_yb171_2023_649nm_770nm.json",    # 171Yb 649/770nm + A constants
    "suryanarayana_yb_2023_581nm.json",    # 174Yb 581nm inner-shell (needs 174Yb 3P1 from atkinson)
    "borisov_yb_1996_581nm.json",          # 581nm isotope shifts all isotopes + 171/173 HFS components
    "kroeze_yb171_2025.json",              # 171Yb compilation (is_recommended=0, goes last)
]

EXPECTED_ERRORS = {
    "no_doi_old_paper.json": 1,  # transition: Rb/null atom not in DB
}


def run_suite(db_path, file_dir, file_order, suite_name):
    print(f"\n{'='*62}")
    print(f"  {suite_name}")
    print(f"{'='*62}")
    create_database(str(db_path))

    total_errors = 0
    unexpected_errors = 0

    with get_connection(str(db_path)) as conn:
        for fname in file_order:
            path = file_dir / fname
            if not path.exists():
                print(f"  MISSING   {fname}")
                continue

            counts, errors = parse_file(conn, path)
            n_inserted = sum(counts.values())
            expected = EXPECTED_ERRORS.get(fname, 0)
            n_unexpected = max(0, len(errors) - expected)

            if not errors:
                status = "OK      "
            elif n_unexpected == 0:
                status = "OK(exp) "   # errors were expected
            else:
                status = "FAIL    "

            print(f"  {status} {fname}")
            print(f"           inserted: {dict(counts)}")
            for e in errors:
                tag = "[expected]" if expected > 0 else "[UNEXPECTED]"
                print(f"           {tag} {e}")
            total_errors   += len(errors)
            unexpected_errors += n_unexpected

    return total_errors, unexpected_errors


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        d_total, d_unexp = run_suite(
            Path(tmp) / "debug.db", JSON_DEBUG, DEBUG_ORDER, "DEBUG FILE SUITE"
        )
        y_total, y_unexp = run_suite(
            Path(tmp) / "yb.db", JSON_YB, YB_ORDER, "YTTERBIUM FILE SUITE"
        )

    print(f"\n{'='*62}")
    print(f"  Total errors : {d_total + y_total}  "
          f"(unexpected: {d_unexp + y_unexp})")
    print(f"{'='*62}\n")
    sys.exit(0 if (d_unexp + y_unexp) == 0 else 1)
