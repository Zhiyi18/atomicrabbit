"""
build_yb_db.py  —  create spectroscopy.db and populate it with Yb I data.

Step 1: NIST ASD  (natural Yb, source_type='nist')
Step 2: Paper JSON files from json_ytterbium/ (isotope-specific, source_type='paper')
"""
import sys, importlib.util, logging
from pathlib import Path

ROOT    = Path(__file__).parent.parent.parent.parent  # yb_test -> db -> atomicrabbit(pkg) -> project root
DB_PKG  = ROOT / "atomicrabbit" / "db"
DB_PATH = ROOT / "atomicrabbit" / "db" / "spectroscopy.db"

# ── direct module loading (avoids triggering physics/__init__ imports) ─────────
def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

db_utils  = _load("atomicrabbit.db.db_utils",         DB_PKG / "db_utils.py")
nist_mod  = _load("atomicrabbit.db.nist.load_tables", DB_PKG / "nist" / "load_tables.py")
parser    = _load("atomicrabbit.db.papers.parser",     DB_PKG / "papers" / "parser.py")

get_connection  = db_utils.get_connection
create_database = db_utils.create_database
load_nist_to_db = nist_mod.load_nist_to_db
parse_file      = parser.parse_file

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

# ── Yb I wavelength ranges to query from NIST ─────────────────────────────────
# Split into ranges to avoid overly large responses; covers all AMO-relevant lines.
NIST_RANGES = [
    (300,  500),   # 399nm 1S0->1P1, 431nm inner-shell
    (500,  700),   # 556nm 1S0->3P1, 578nm clock, 649nm 3P0->3S1
    (700,  900),   # 770nm 3P2->3S1
    (900, 2000),   # 1389nm 3P0->3D1 and other IR repumps
]

# ── Paper files, in dependency order ──────────────────────────────────────────
JSON_YB = ROOT / "atomicrabbit" / "db" / "papers" / "json_ytterbium"

YB_PAPER_ORDER = [
    "nudat3_yb_nuclear_moments.json",      # 173Yb mu_I; atoms only
    "kleinert_yb_2016_399nm.json",         # 399nm 1S0->1P1 isotope shifts + A(1P1)
    "atkinson_yb_2019_556nm.json",         # 556nm 1S0->3P1 isotope shifts + A(3P1)
    "kawasaki_yb_2024_431nm.json",         # 431nm inner-shell isotope shifts
    "jones_yb171_2023_556nm.json",         # 171Yb 3P1 high-precision A constant
    "kim_yb171_2021_clock.json",           # 171Yb clock transition absolute frequency
    "qiao_yb_2023_770nm_mts.json",         # 174Yb 770nm absolute frequency
    "qiao_yb171_2023_649nm_770nm.json",    # 171Yb 649/770nm + A(3S1, 3P2)
    "suryanarayana_yb_2023_581nm.json",    # 174Yb 581nm inner-shell scheme
    "borisov_yb_1996_581nm.json",          # 581nm isotope shifts all isotopes + 171/173 HFS components
    "kroeze_yb171_2025.json",              # 171Yb reference compilation (is_recommended=0)
]


def main():
    print(f"Database: {DB_PATH}\n")
    create_database(str(DB_PATH))

    # ── Step 1: NIST ──────────────────────────────────────────────────────────
    print("=" * 60)
    print("  STEP 1 - NIST ASD (Yb I, natural isotope mix)")
    print("=" * 60)
    with get_connection(str(DB_PATH)) as conn:
        for lo, hi in NIST_RANGES:
            try:
                load_nist_to_db(conn, element="Yb", mass_number=None,
                                wavelength_min_nm=lo, wavelength_max_nm=hi)
            except Exception as e:
                print(f"  WARNING: NIST query {lo}-{hi} nm failed: {e}")

    # ── Step 2: Papers ────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("  STEP 2 - Paper JSON files (isotope-specific)")
    print("=" * 60)
    total_errors = 0
    with get_connection(str(DB_PATH)) as conn:
        for fname in YB_PAPER_ORDER:
            path = JSON_YB / fname
            if not path.exists():
                print(f"  MISSING  {fname}")
                continue
            counts, errors = parse_file(conn, path)
            status = "OK  " if not errors else "FAIL"
            print(f"  [{status}] {fname}")
            print(f"         inserted: {dict(counts)}")
            for e in errors:
                print(f"         ERROR: {e}")
            total_errors += len(errors)

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    with get_connection(str(DB_PATH)) as conn:
        n_atoms   = conn.execute("SELECT COUNT(*) FROM atoms").fetchone()[0]
        n_papers  = conn.execute("SELECT COUNT(*) FROM papers").fetchone()[0]
        n_levels  = conn.execute("SELECT COUNT(*) FROM energy_levels").fetchone()[0]
        n_trans   = conn.execute("SELECT COUNT(*) FROM transitions").fetchone()[0]
        n_nist_t  = conn.execute(
            "SELECT COUNT(*) FROM transitions WHERE source_type='nist'").fetchone()[0]
        n_paper_t = conn.execute(
            "SELECT COUNT(*) FROM transitions WHERE source_type='paper'").fetchone()[0]
        n_nist_l  = conn.execute(
            "SELECT COUNT(*) FROM energy_levels WHERE source_type='nist'").fetchone()[0]
        n_paper_l = conn.execute(
            "SELECT COUNT(*) FROM energy_levels WHERE source_type='paper'").fetchone()[0]

    print(f"  DB summary:")
    print(f"    atoms          : {n_atoms}")
    print(f"    papers         : {n_papers}")
    print(f"    energy_levels  : {n_levels}  (nist={n_nist_l}, paper={n_paper_l})")
    print(f"    transitions    : {n_trans}  (nist={n_nist_t}, paper={n_paper_t})")
    print(f"  Paper errors     : {total_errors}")
    print("=" * 60)

    return total_errors


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
