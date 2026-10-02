"""Run the paper parser on all Ba ion JSON files and print a summary of errors."""
import logging
from pathlib import Path
from atomicrabbit.db.db_utils import get_connection
from atomicrabbit.db.papers.parser import parse_directory

logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')

DB_PATH = "atomicrabbit/db/spectroscopy.db"
JSON_DIR = "atomicrabbit/db/papers/json_barium_ion"

with get_connection(DB_PATH) as conn:
    results = parse_directory(conn, JSON_DIR)

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
total_errors = 0
for fname, (counts, errors) in sorted(results.items()):
    print(f"\n{fname}")
    print(f"  inserted: {counts}")
    if errors:
        for e in errors:
            short = e.split('] ', 1)[-1]
            print(f"  ERROR: {short}")
    total_errors += len(errors)

print(f"\nTotal errors: {total_errors}")
