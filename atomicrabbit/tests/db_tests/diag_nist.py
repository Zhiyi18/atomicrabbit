"""Diagnose NIST query column names and row format for Ba II and Yb I."""
from astroquery.nist import Nist
import astropy.units as u

for linename, lo, hi in [("Ba II", 450, 500), ("Yb I", 390, 410)]:
    print(f"\n{'='*60}")
    print(f"{linename}  [{lo}-{hi} nm]")
    print(f"{'='*60}")
    try:
        t = Nist.query(lo * u.nm, hi * u.nm, linename=linename)
        print(f"Columns: {t.colnames}")
        print(f"Rows: {len(t)}")
        for i, row in enumerate(t[:3]):
            print(f"\n--- row {i} ---")
            for col in t.colnames:
                print(f"  {col!r:30} = {str(row[col])!r}")
    except Exception as e:
        print(f"ERROR: {e}")
