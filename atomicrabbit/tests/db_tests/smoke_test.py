import sys, importlib.util, json
from pathlib import Path

ROOT   = Path("C:/Users/USER-01/Documents/Claude/atomicrabbit")
DB_PKG = ROOT / "atomicrabbit" / "db"

def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

db_utils     = _load("atomicrabbit.db.db_utils",               DB_PKG / "db_utils.py")
nist_mod     = _load("atomicrabbit.db.nist_tools.load_tables", DB_PKG / "nist_tools" / "load_tables.py")
paper_parser = _load("atomicrabbit.db.json.paper_parser",      DB_PKG / "json" / "paper_parser.py")

print("db_utils  OK:", hasattr(db_utils, "get_connection"))
print("nist_mod  OK:", hasattr(nist_mod, "load_nist_to_db"))
print("parser    OK:", hasattr(paper_parser, "parse_file"))

for el in ["Ba II", "Yb I"]:
    m = DB_PKG / "json" / el / "_manifest.json"
    data = json.load(open(m))
    print(f"{el}: linename={data['linename']}, masses={data['nist_mass_numbers']}")

print("All OK")
